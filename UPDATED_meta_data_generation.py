import os
import json
from tqdm import tqdm
from Bio import Entrez
from openai import OpenAI
import re
import time
from dotenv import load_dotenv
from xml.etree import ElementTree as ET

# Initialize OpenAI client
load_dotenv()
client = OpenAI(api_key=os.getenv("OPEN_AI_API"))

# Always set your email for Entrez
Entrez.email = "priyankp@wisdomsquare.net"


# ---- Function to fetch title + abstract ----
def fetch_extract_and_abstract(pmid):
    pmid = str(pmid)
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="abstract", retmode="xml")
    records = Entrez.read(handle)
    handle.close()

    article = records["PubmedArticle"][0]["MedlineCitation"]["Article"]

    # Title
    title = article.get("ArticleTitle", "")

# Abstract
    abstract_text = ""
    if "Abstract" in article:
        abstract_parts = article["Abstract"].get("AbstractText", [])
        abstract_text = " ".join(str(p) for p in abstract_parts)

    if abstract_text:
        return {"pmid": pmid, "title": title, "abstract": abstract_text}

    # --- Fallback: use PMC full text as "abstract" ---
    def _strip_ns(root):
        for el in root.iter():
            if "}" in el.tag:
                el.tag = el.tag.split("}", 1)[1]

    def _get_text(el):
        return " ".join("".join(el.itertext()).split())

    try:
        time.sleep(0.34)
        lnk = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid, linkname="pubmed_pmc")
        linkset = Entrez.read(lnk)
        lnk.close()

        pmc_numeric = None
        for lsdb in linkset[0].get("LinkSetDb", []):
            for link in lsdb.get("Link", []):
                pmc_numeric = link.get("Id"); break
            if pmc_numeric: break

        if not pmc_numeric:
            return {"pmid": pmid, "title": title, "abstract": ""}

        pmcid = f"PMC{pmc_numeric}"
        time.sleep(0.34)
        h2 = Entrez.efetch(db="pmc", id=pmcid, rettype="full", retmode="xml")
        xml_bytes = h2.read()
        h2.close()

        root = ET.fromstring(xml_bytes)
        _strip_ns(root)

        # Prefer body-only text (skips references/back matter)
        body = root.find(".//body")
        if body is not None:
            paras = [_get_text(p) for p in body.findall(".//p")]
            full_text = " ".join(p for p in paras if p) or _get_text(body)
            return {"pmid": pmid, "title": title, "abstract": full_text}

        # Fallback: whole doc text if body missing
        return {"pmid": pmid, "title": title, "abstract": _get_text(root)}

    except Exception:
        # If anything fails, just return what we have
        return {"pmid": pmid, "title": title, "abstract": ""}
    

# ---- Function to Check None type Focus ----
    
def normalize_text(text):
    """Converts text to a consistent format for matching."""
    # 1. Convert to lowercase
    # 2. Replace all hyphens with spaces
    # 3. Collapse multiple spaces into a single space
    normalized = text.lower().replace('-', ' ')
    return re.sub(r'\s+', ' ', normalized).strip()

def find_synonyms_in_text(synonyms_data, text):
    """
    Finds synonyms in text, handling hyphen and space variations.

    Args:
        synonyms_data (set): A set with a semicolon-separated string of synonyms.
        text (str): The text abstract to search within.

    Returns:
        set: A set of the original synonyms that were found.
    """
    # 1. Parse the synonym string into a clean list
    synonym_list = [s.strip() for s in synonyms_data.split(';')]
    print(synonym_list)

    # 2. Create a mapping from a normalized synonym to its original form
    # This allows us to find a match and know which original term it came from.
    # Example: {'acetyl l carnitine': 'Acetyl L Carnitine', 
    #           'β acetyl l carnitine': 'β-acetyl-L-carnitine'}
    normalized_map = {normalize_text(s): s for s in synonym_list}

    # 3. Normalize the entire abstract for a consistent search
    normalized_abstract = normalize_text(text)

    # 4. Find matches
    found_synonyms = set()
    for norm_syn, original_syn in normalized_map.items():
        # Use word boundaries (\b) to ensure we match whole phrases
        pattern = r'\b' + re.escape(norm_syn) + r'\b'
        if re.search(pattern, normalized_abstract):
            found_synonyms.add(original_syn) # Add the original form to the results

    #print(">>>>>Found Synonyms:",found_synonyms)
            
    return found_synonyms


# ---- Function to extract metadata ----
def extract_metadata(synonym, title, abstract, pubmed_type, focus_status):
    metadata_schema = {
        "type": "object",
        "properties": {
            "duration_days": {
                "type": ["integer", "null"],
                "description": "Convert weeks/months/years to days, null if not reported."
            },
            "sample_size": {
                "type": ["integer", "null"],
                "description": "Number of participants, null if not reported."
            },
            "sample_gender": {
                "type": "array",
                "items": {
                    "type": "string",
                    "enum": ["male", "female", "other", "unspecified"]

                },
                "description": "Participant genders included."
            },
            "species": {
                "type": "array",
                "items": {
                    "type": "string",
                    "enum": ["humans","animals","other"]
                },
                    "description": "Species studied in the article.\
                                    Categories include:\
                                    - 'humans': Human participants of any age group (adults, children, adolescents).\
                                    - 'animals': animals such as rodents, dogs, cows, poultry, fish, insects, etc.\
                                    - 'other': Includes eveything else other than humans and animals."
},
            "population": {
                "type": "string",
                "description": "Short descriptive text of the study population."
            },
            "study_type": {
                "type": "array",
                "items": {
                    "type": "string",
                    "enum": [
                        "randomized controlled trial",
                        "non-randomized trial",
                        "open-label trial",
                        "single-blind trial",
                        "double-blind trial",
                        "triple-blind trial",
                        "crossover trial",
                        "parallel group trial",
                        "factorial design",
                        "cluster randomized trial",
                        "unspecified"
                    ]
                },
                "description": f"Type of study or trial design used in the research. \
                The aticle type is {pubmed_type}, which provides context for interpreting the study type. \
                Categories include: \
                - 'randomized controlled trial (rct)': Participants are randomly assigned to groups (e.g., treatment vs placebo). \
                - 'non-randomized trial': Groups assigned by investigator judgment or participant choice. \
                - 'open-label trial': Both participants and researchers know which treatment is given. \
                - 'single-blind trial': Either participants or investigators are unaware of treatment assignment. \
                - 'double-blind trial': Both participants and investigators are unaware of treatment assignment. \
                - 'triple-blind trial': Participants, investigators, and data analysts/statisticians are blinded. \
                - 'crossover trial': Participants receive multiple interventions sequentially (each subject acts as their own control). \
                - 'parallel group trial': Each group receives a different intervention simultaneously. \
                - 'factorial design': Tests two or more interventions in combination (e.g., A, B, A+B). \
                - 'cluster randomized trial': Groups (not individuals) are randomized (e.g., schools, clinics). \
                - 'unspecified': Trial aticle type is not clinical or cannot be determined."
            },
            "focus": {
                "type": "array",
                "items": {
                    "type": "string",
                    "enum": ["primary", "secondary"]
                },
                "description": f"""Determine the focus of the study into one of the below category.
                    - primary: If the {synonym} is the main focus of the study. The study investigates its composition, effects, mechanisms, applications, or biological activity.
                    - secondary: If the {synonym} is mentioned along with the other ingredients but not main subject of the study. The study’s main conclusions or experiments are not about this {synonym}."""
            },
            "benefits": {
                "type": "array",
                "items": {
                    "type": "string"
                },
                "description": f"One concise line with rich context describing benefits of {synonym}. Use the same wordings mentioned in the article as much as possible without using acronyms or short forms. For example, It can help in improving sleep quality."
            },
            "synergies_interactions_positive": {
                "type": "array",
                "items": {
                    "type": "string"
                },
                "description": f"One concise line with rich context describing positive synergies or beneficial interactions of {synonym}  with other compounds, ingredients, or treatments. Use the same wordings mentioned in the article as much as possible without using acronyms or short forms. For example, 'It has enhanced calming effect when combined with lavender oil."
            },
            "synergies_interactions_negative": {
                "type": "array",
                "items": {
                    "type": "string"
                },
                "description": f"One concise line with rich context describing negative synergies interactions of {synonym} with other compounds, ingredients, or treatments, including reduced efficacy or harmful effects. Use the same wordings mentioned in the article as much as possible without using acronyms or short forms. For example, 'It may increase drowsiness when combined with alcohol."
            },
            "safety_side_effects": {
                "type": "array",
                "items": {
                    "type": "string"
                },
                "description": f"One concise line with rich context describing safety concerns, risks, toxicity, or side effects associated with the {synonym}. Use the same wordings mentioned in the article as much as possible without using acronyms or short forms. For example, It may cause mild stomach upset if taken in high doses."
            },
            "interventions": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                        "ingredient": {
                            "type": "string",
                            "description": "Name of the intervention substance (e.g.,Ashwagandha)."
                        },
                        "daily_dosage": {
                    "type": ["number", "null"],
                    "description": "Total daily dosage normalized according to the unit reported in the paper. Example: '600 mg twice daily' → 1200, '0.5 g/day' → 0.5 (if unit is g), '10 ml/day' → 10. If normalization is not possible, set null."
                },
                "units": {
                    "type": "string",
                    "description": "Unit of the dosage as reported in the paper. Examples: 'mg', 'g', 'ml', 'IU', 'mcg', 'microg/ml' etc."
                },
                "original_text": {
                    "type": "string",
                    "description": "Exact dosage text as reported in the paper."
                }
                },
                "required": ["ingredient", "dosage_value", "units", "original_text"]

                    }
                    },
            "usage": {
                    "type": "array",
                    "items": {
                        "type": "string",
                        "enum": ["ingestion", "inhalation", "topical", "injection", "unspecified"]
                    },
                    "description": f"Route of administration of {synonym} or the treatment containing it. \
                Look for mentions of how the substance was administered to participants in the methods or intervention description. All specific routes are grouped into one of five categories: \
                - 'ingestion': Includes 'oral', 'capsule', 'pill', 'tablet', 'drink', 'solution', 'sublingual', 'chew', 'swallowed, or any route involving swallowing or absorption through the digestive or mucosal system. \
                - 'inhalation': Includes inhaled', 'breathed', 'vapor', 'aromatherapy', 'smoke', 'diffused', 'olfactory exposure', 'nasal spray', 'nebulizer' or any route where the compound is inhaled through the respiratory system. \
                - 'topical': Includes 'applied', 'ointment', 'cream', 'serum', 'oil', 'massage', 'on skin', 'on hair', 'rubbed', 'transdermal' or ocular/ophthalmic routes where the compound is applied externally or absorbed through the skin or eyes. \
                - 'injection': Includes 'IV', 'intramuscular', 'subcutaneous', 'injected', 'parenteral', 'infused', 'intravenous' and all parenteral routes where the substance is administered via needle or infusion. \
                - 'unspecified': When the route of administration is not mentioned explicitly.",
                "minItems": 1,
                },
            "conditions": {
                "type": "array",
                    "items": {"type": "string"},
                "description": "Medical conditions studied.",
                "default": []
        },
            "biomarkers": {
                "type": "array",
                    "items": {"type": "string"},
                "description": "Physiological measures studied.",
                "default": []
        },
            "functions": {
                "type": "array",
                    "items": {"type": "string"},
            "description": "Functional domains studied (e.g., cognition, sleep).",
            "default": []
        },
            "purpose": {
                "type": "string",
                "description": "Short of description main purpose or aim of the study."
            },
            "conclusion":{
                "type": "string",
                "description": f"The main conclusion of the study focusing on {synonym}, summarized in one to two concise sentences with very rich context. If no explicit conclusion is provided, summarize the key finding instead. Do not include background, methods, or introduction details. Do not be generic; Use the same wordings mentioned in the article as much as possible without using acronyms or short forms"
            },
            "outcomes": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "name": {
                            "type": "string",
                            "description": "Canonicalized outcome name."
                        },
                        "domain": {
                            "type": "string",
                            "enum": ["condition", "function", "biomarker"],
                            "description": "Outcome classification domain."
                        },
                        "type": {
                            "type": "string",
                            "enum": ["primary", "secondary"],
                            "description": "Whether the outcome was primary or secondary."
                        },
                        "result": {
                            "type": "string",
                            "enum": ["improved", "worsened", "no_effect", "mixed", "not_reported"],
                            "description": "Reported result direction."
                        }
                    },
                    "description":"Must include all the values in conditions, biomarkers and functions with name, domain, type, result.",
                    "required": ["name", "domain", "type", "result"]
                }
            },
            "diseases": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Standardized disease names."
            },
            "symptoms": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Standardized symptoms studied."
            },
            "keywords": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Keywords or important terms."
            },
            "location": {
                "type": ["string", "null"],
                "description": "Geographic country name of the study with normalization."
            },
            "mechanism": {
                "type": "array",
                "items": {"type": "string"},
                "description": "Proposed biological mechanisms in short terms only."
            }
        },
        "required": ["interventions","outcomes", "conditions", "biomarkers", "functions","usage"]
    }

    

    if focus_status:
        metadata_schema["properties"]["focus"] = {
                                            "type": "string",
                                            #"description": "Focus determination skipped since status=False.",
                                            "default": "none"
                                        }

    print(">>>>>Schema:",metadata_schema)

    system_prompt = """
    You are a scientific research metadata extractor. 

    Extract the following metadata from the research paper text provided. 
    Return it in JSON format exactly as specified below. 

    Normalization Rules for all the key value pairs:
    - Convert all text to lowercase.
    - Replace Greek letters with their English names (e.g., α → alpha, β → beta, γ → gamma).
    - Transliterate all other Unicode characters to ASCII equivalents (e.g., µ → u, ± → +/-)
    - Remove trademark and copyright symbols (® ™ © ℠).
    - Keep meaningful scientific units and math symbols (+, -, %, mg, ml, um, etc.).
    """

    user_prompt = f"""
    Extract structured metadata from the following PubMed paper of ingredient {synonym}:

    ---
    Title: {title}
    Abstract: {abstract}
    ---
    """

    
    response = client.chat.completions.create(
        model="gpt-4o",
        messages=[
            {"role": "system", "content": system_prompt},
            {"role": "user", "content": user_prompt}
        ],
        temperature=0,
        response_format={
            "type": "json_schema",
            "json_schema": {
                "name": "metadata_schema",
                "schema": json.loads(json.dumps(metadata_schema)),
            }
        }
    )

    result = json.loads(response.choices[0].message.content)

    if focus_status:
        # Force focus to None regardless of GPT output
        result["focus"] = "none"

    return result


def process_pmids(root_names, search_terms, synonyms, pmids, pubmed_types,output_file="pubmed_metadata.json"):
    results = []

    # Load existing results if file already exists
    if os.path.exists(output_file):
        with open(output_file, "r") as f:
            try:
                results = json.load(f)
            except json.JSONDecodeError:
                results = []

    # Build lookup of already processed PMIDs
    done_pmids = {(r.get("PMID"),r.get("root_name")) for r in results}

    # Process remaining
    for root_name, search_term ,synonym, pmid, pubmed_type in tqdm(zip(root_names, search_terms, synonyms, pmids,pubmed_types), total=len(pmids), desc="Processing PMIDs", unit="pmid"):
        if (pmid,root_name) in done_pmids:
            continue  # Skip already processed

        try:
            paper = fetch_extract_and_abstract(pmid)
            print(pmid)
            # checking for no abstracts
            if paper['abstract'].strip():
                print("Synonyms",synonym)
                text = paper["title"] + paper["abstract"]
                # Set focus_status to TRUE for No Match Abstracts
                match = find_synonyms_in_text(synonym, text)
                focus_status = (match == {} or match == {''})

                print(">>>>>>>>>>Focus",focus_status)
                metadata_json = extract_metadata(synonym, paper["title"], paper["abstract"],pubmed_type,focus_status)
                result = {"root_name": root_name, "search_term": search_term,"synonyms" : synonym, "PMID": pmid, "pubmed_type": pubmed_type, "metadata": metadata_json}
            else:
                continue
        except Exception as e:
            result = {"root_name": root_name, "search_term": search_term, "synonyms" : synonym, "PMID": pmid, "pubmed_type": pubmed_type,"Error": str(e)}

        results.append(result)

        # Save incrementally after each PMID
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2)

    return results
