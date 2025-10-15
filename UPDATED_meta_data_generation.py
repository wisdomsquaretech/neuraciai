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
    medline_citation = records["PubmedArticle"][0]["MedlineCitation"]

    # Title
    title = article.get("ArticleTitle", "")

# Abstract
    abstract_text = ""
    if "Abstract" in article:
        abstract_parts = article["Abstract"].get("AbstractText", [])
        abstract_text = " ".join(str(p) for p in abstract_parts)
    
    keywords = []
    if "KeywordList" in medline_citation:
        for keyword_list in medline_citation["KeywordList"]:
            for kw in keyword_list:
                keywords.append(str(kw))

    if abstract_text:
        return {"pmid": pmid, "title": title, "abstract": abstract_text, "keywords":keywords}

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
        return {"pmid": pmid, "title": title, "abstract": _get_text(root), "keywords":keywords}

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
    print(">>>>>> Synonyms List",synonym_list)

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

    print(">>>>>Found Synonyms:",found_synonyms)
            
    return found_synonyms


# ---- Function to extract metadata ----
def extract_metadata(synonym, title, abstract, pubmed_type, focus_status):

    metadata_schema = {
        "type": "object",
        "properties": {
            "duration_days": {
                    "oneOf": [
                        {
                        "type": "integer",
                        "description": "Study or treatment duration expressed in days. Convert weeks/months/years to total days."
                        },
                        {
                        "const": "not mentioned",
                        "description": "Duration information was not reported or cannot be determined from the text."
                        },
                    ],
                    "description": "Duration of the study or experimental intervention, normalized to days."
                    },
            "sample_size": {
                    "oneOf": [
                        {
                        "type": "integer",
                        "description": "Number of participants included in the study."
                        },
                        {
                        "const": "not mentioned",
                        "description": "The sample size was not reported or cannot be determined from the study."
                        },
                        ],
                    "description": "Total number of participants or samples in the study."
                    },
            "sample_gender": {
                "type": "array",
                "items": {
                    "oneOf": [
                    {
                        "const": "male",
                        "description": "Participant identified as male."
                    },
                    {
                        "const": "female",
                        "description": "Participant identified as female."
                    },
                    {
                        "const": "other",
                        "description": "Participant gender is not male or female (includes nonbinary or diverse identities)."
                    },
                    {
                        "const": "not mentioned",
                        "description": "Gender information was not provided or not stated in the source."
                    }
                    ]
                },
                "description": "Participant genders included in the study.",
                "minItems": 1,
                "uniqueItems": True,
                "default": ["not mentioned"]
                },
            "species": {
                "type": "array",
                "items": {
                "oneOf": [
                {
                "const": "humans",
                "description": "Studies involving human participants, volunteers, or patients."
                },
                {
                "const": "animals",
                "description": "Studies involving non-human animals such as mice, rats, rabbits, or other vertebrates/invertebrates."
                },
                {
                "const": "cell lines",
                "description": "Immortalized or continuously cultured cells derived from humans or animals, used for in vitro experiments."
                },
                {
                "const": "primary cells",
                "description": "Cells directly isolated from tissues of humans, animals, or plants, cultured briefly for experiments."
                },
                {
                "const": "tissues",
                "description": "Ex vivo samples of tissues from humans, animals, or plants used in experiments or analysis."
                },
                {
                "const": "organoids",
                "description": "3D miniaturized and simplified versions of organs grown in vitro from stem cells to model organ function or disease."
                },
                {
                "const": "microorganisms",
                "description": "Microscopic organisms including bacteria, fungi, protozoa, algae, or viruses."
                },
                {
                "const": "plants",
                "description": "Whole plants or plant parts (leaves, roots, seeds, etc.) used in experiments or studies."
                },
                {
                "const": "computer models",
                "description": "Computational or in silico models, simulations, or algorithmic experiments."
                },
                {
                "const": "other",
                "description": "Any experimental system or species that does not fall into the predefined categories."
                },
                {
                "const": "not mentioned",
                "description": "Species or biological system was not specified in the article."
                }
                ]
                },
                "description": "Species or biological systems studied in the article.",
                "minItems": 1,
                "uniqueItems": True,
                "default": ["not mentioned"]
                },
            "experimental_model": {
                "type": "array",
                "items": {
                    "oneOf": [
                    {
                        "const": "in vitro",
                        "description": "Experiments conducted outside living organisms, typically in test tubes, culture plates, or petri dishes."
                    },
                    {
                        "const": "ex vivo",
                        "description": "Experiments conducted on tissues taken from living organisms but outside the organism."
                    },
                    {
                        "const": "in vivo",
                        "description": "Experiments conducted within living organisms (humans, animals, or plants)."
                    },
                    {
                        "const": "in silico",
                        "description": "Studies conducted using computer simulations, computational modeling, or digital analyses."
                    },
                    {
                        "const": "other",
                        "description": "Any experimental approach or environment not covered by the above categories like ex vivo, in vivo etc."
                    },
                    {
                        "const": "not mentioned",
                        "description": "The experimental model type was not reported."
                    }
                    ]
                },
                "description": "Type(s) of experimental models used in the study.",
                "minItems": 1,
                "uniqueItems": True,
                "default": ["not mentioned"]
                },
            "population": {
                "oneOf": [
                    {
                    "type": "string",
                    "description": "Short descriptive text of the study population, including condition/diagnosis and any defining traits or subgroups (e.g., early-treatment vs therapeutic-treatment), but excluding participant numbers, age, or sex unless reported as defining characteristics."
                    },
                ],
                "description": "Text description of the study population.",
                },
            "study_type": {
                "type": "array",
                "items": {
                    "oneOf": [
                    {
                        "const": "randomized controlled trial",
                        "description": "Participants are randomly assigned to different intervention or control groups to minimize bias and allow direct comparison of outcomes. Synonyms: RCT, placebo-controlled, controlled trial, multicenter randomized."
                    },
                    {
                        "const": "non-randomized trial",
                        "description": "Participants are assigned to groups without randomization, based on investigator judgment or predefined criteria. Synonyms: non-randomized, nonrandomised, investigator assigned."
                    },
                    {
                        "const": "open-label trial",
                        "description": "Both participants and investigators know which intervention is being administered. Synonyms: open-label, open trial, unblinded."
                    },
                    {
                        "const": "single-blind trial",
                        "description": "Either participants or investigators (usually participants) are unaware of which intervention has been assigned. Synonyms: single-blind, single-masked."
                    },
                    {
                        "const": "double-blind trial",
                        "description": "Both participants and investigators are unaware of treatment assignments to minimize bias. Synonyms: double-blind, double-masked, placebo-controlled double-blind."
                    },
                    {
                        "const": "triple-blind trial",
                        "description": "Participants, investigators, and data analysts/statisticians are all blinded to treatment assignments. Synonyms: triple-blind, triple-masked."
                    },
                    {
                        "const": "crossover trial",
                        "description": "Each participant receives multiple interventions in sequence, acting as their own control. Synonyms: crossover design, cross-over, cross over study."
                    },
                    {
                        "const": "parallel group trial",
                        "description": "Participants are assigned to one intervention group and remain there; groups are compared at the end. Synonyms: parallel group, parallel-arm, parallel design, parallel assignment."
                    },
                    {
                        "const": "factorial design",
                        "description": "Two or more interventions are tested simultaneously to evaluate individual and combined effects. Synonyms: factorial trial, 2x2 design, multifactorial."
                    },
                    {
                        "const": "cluster randomized trial",
                        "description": "Groups or clusters are randomized instead of individuals. Synonyms: cluster randomized, group randomized, cluster trial, community randomized."
                    },
                    {
                        "const": "adaptive trial",
                        "description": "Trial design allows planned modifications based on interim analyses without undermining validity. Synonyms: adaptive design, adaptive randomization, adaptive study."
                    },
                    {
                        "const": "pragmatic clinical trial",
                        "description": "Designed to evaluate interventions in routine clinical practice, emphasizing real-world effectiveness. Synonyms: pragmatic trial, effectiveness trial, practical clinical trial."
                    },
                    {
                        "const": "other",
                        "description": "Trial design cannot be classified into any of the predefined types such as 'randomized controlled trial', 'single-blind trial', 'double-blind trial', etc."
                    },
                    {
                        "const": "not mentioned",
                        "description": "Study type was not reported from the article."
                    }
                    ]
                },
                "description": f"Type(s) of study or trial design used in the research. Given that article type is {pubmed_type}",
                "minItems": 1,
                "uniqueItems": True,
                "default": ["not mentioned"]
                },
            "focus": {
            "type": "array",
            "items": {
                "oneOf": [
                {
                    "const": "primary",
                    "description": f"The {synonym} is the **main focus** of the study. The study is primarily designed to investigate this {synonym}’s composition, structure, biological activity, health effects, therapeutic applications, dosage, or safety. It typically appears as the main intervention, central subject of analysis, or key variable in study objectives and conclusions."
                },
                {
                    "const": "secondary",
                    "description": f"The {synonym} is mentioned but not the main subject of the study. This includes situations where the {synonym} is part of a multi-ingredient formula, appears only in background context or literature references, is used as a comparator or control, or the study’s main conclusions concern something else."
                },
                ]
            },
            "description": "Determines the focus of the study for the ingredient(s).",
            "minItems": 1,
            "uniqueItems": True,
            },
            "benefits": {
            "type": "array",
            "items": {
                "oneOf": [
                {
                    "type": "string",
                    "description": f"A concise functional or clinical benefit of the {synonym} (1-3 words), extracted directly from the study. Examples: 'blood pressure reduction', 'anti-inflammatory effect', 'improved memory'."
                },
                {
                    "const": "not mentioned",
                    "description": f"No functional or clinical benefits of the {synonym} were reported or could be extracted from the study."
                }
                ]
            },
            "description": "Functional or clinical benefits of the ingredient reported in the study. Return only measurable outcomes or effects; do not include qualitative descriptors like 'safe' or 'promising'.",
            "minItems": 1,
            "uniqueItems": True,
            "default": ["not mentioned"]
            },

          "synergies_interactions_positive": {
            "type": "array",
            "items": {
                "oneOf": [
                {
                    "type": "string",
                    "description": f"A concise statement describing a positive synergy, beneficial interaction, or enhanced effect of the {synonym} with other compounds, ingredients, or treatments. Use wording from the study; avoid acronyms or short forms. Only include functional or therapeutic effects, not general qualitative descriptors unless they describe measurable outcomes. Example: 'Enhanced calming effect when combined with lavender oil.'"
                },
                {
                    "const": "not mentioned",
                    "description": f"No positive synergies or beneficial interactions of {synonym} were reported or could be determined from the study."
                }
                ]
            },
            "description": "List of positive synergies or beneficial interactions involving the ingredient.",
            "minItems": 1,
            "uniqueItems": True,
            "default": ["not mentioned"]
            },
            "synergies_interactions_negative": {
                "type": "array",
                "items": {
                    "oneOf": [
                    {
                        "type": "string",
                        "description": f"A concise statement describing a negative interaction, reduced efficacy, or adverse effect of the {synonym} with other compounds, ingredients, or treatments. Use wording from the study; avoid acronyms or short forms. Only include measurable functional or clinical outcomes, not general qualitative descriptors unless tied to a functional outcome. Example: 'Co-administration with certain antidepressants increased risk of insomnia.'"
                    },
                    {
                        "const": "not mentioned",
                        "description": f"No negative interactions or adverse effects of {synonym} were reported or could be determined from the study."
                    }
                    ]
                },
                "description": "List of negative synergies, adverse interactions, or reduced effects involving the ingredient.",
                "minItems": 1,
                "uniqueItems": True,
                "default": ["not mentioned"]
                },
            "safety_side_effects": {
                "type": "array",
                "items": {
                    "oneOf": [
                    {
                        "type": "string",
                        "description": f"A concise functional or clinical adverse effect caused by the {synonym}, either alone or in combination with other compounds or treatments. Examples: 'nausea', 'liver toxicity', 'headache'. Only include measurable negative outcomes, not general descriptors like 'safe' or 'well-tolerated'."
                    },
                    {
                        "const": "not mentioned",
                        "description": f"No adverse effects of {synonym} were reported or could be determined from the study."
                    }
                    ]
                },
                "description": "Functional or clinical adverse effects of the ingredient reported in the study.",
                "minItems": 1,
                "uniqueItems": True,
                "default": ["not mentioned"]
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
                    "oneOf": [
                    {
                        "const": "ingestion",
                        "description": "Oral, sublingual, buccal, rectal, vaginal, oral-transmucosal, or any route involving swallowing or absorption through the digestive or mucosal system."
                    },
                    {
                        "const": "inhalation",
                        "description": "Aromatherapy, vapor, aerosol, nasal spray, nebulizer, or any route where the compound is inhaled through the respiratory system."
                    },
                    {
                        "const": "topical",
                        "description": "Dermal application (cream, serum, oil), transdermal patches, or ocular/ophthalmic routes where the compound is applied externally or absorbed through the skin or eyes."
                    },
                    {
                        "const": "injection",
                        "description": "All parenteral routes (intravenous, intramuscular, subcutaneous, etc.) where the substance is administered via needle or infusion."
                    },
                    {
                        "const": "other",
                        "description": "Route of administration does not fall into the predefined categories."
                    },
                    {
                        "const": "not mentioned",
                        "description": "The route of administration was not reported or could not be determined from the study."
                    }
                    ]
                },
                "description": f"Route(s) of administration of the {synonym} or treatment.",
                "minItems": 1,
                "uniqueItems": True,
                "default": ["not mentioned"]
                },
            "conditions": {
                "type": "array",
                    "items": {"type": "string"},
                "description": "Medical conditions, diseases, or health statuses that are the focus of the study population. Include all conditions explicitly mentioned in the title, abstract, or methods section that the study is investigating. If multiple conditions are included, list each one separately.",
                "default": [] 
        },
            "biomarkers": {
                "type": "array",
                    "items": {"type": "string"},
                "description": "Physiological, biochemical, or molecular measures that are studied or reported in the trial. Include all relevant markers explicitly mentioned in the title, abstract, or methods section. Examples include blood pressure, heart rate, hormone levels, inflammatory markers, or neurochemical measures.",
                "default": []
        },
            "functions": {
                "type": "array",
                    "items": {"type": "string"},
            "description": "Functional or behavioral domains that are studied or assessed in the trial. Include all domains explicitly mentioned in the title, abstract, or methods section, such as cognition, sleep, mood, memory, attention, pain, or motor function. Return concise keywords for each domain.",
            "default": []
        },
            "purpose": {
                "type": "string",
                "description": "A concise statement describing the main objective, aim, or purpose of the study, as explicitly reported in the title, abstract, or introduction. Do not infer outcomes or results; focus only on the study’s stated intent."

            },
            "conclusion":{
                "type": "string",
                "description": f"""The **main conclusion** or key finding of the study specifically related to the {synonym}, summarized in one to two concise sentences with rich context. If no explicit conclusion is provided, extract the most important result from the results section. Do not include background, methods, or introductory details. Avoid generic statements; use wording directly from the article whenever possible, without acronyms or short forms."""

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
            "description": "RAG-optimized keywords extracted from the PubMed abstract, capturing the main biomedical concepts, diseases, genes, proteins, or processes. These keywords are designed to improve semantic retrieval, context-aware search, and accurate document linking in retrieval-augmented generation workflows.",
            "minItems": 5,
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


def process_pmids(root_names, synonyms, pmids, pubmed_types,output_file="pubmed_metadata.json"):
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
    for root_name ,synonym, pmid, pubmed_type in tqdm(zip(root_names, synonyms, pmids,pubmed_types), total=len(pmids), desc="Processing PMIDs", unit="pmid"):
        if (pmid,root_name) in done_pmids:
            continue  # Skip already processed

        try:
            paper = fetch_extract_and_abstract(pmid)
            #print(pmid)
            # checking for no abstracts
            if paper['abstract'].strip():
                print("Synonyms",synonym)
                abstract_with_keywords = paper["abstract"] + "\nKeywords: " + ",".join(paper["keywords"])
                text = paper["title"] + "\n" +abstract_with_keywords
                print(text)
                # Set focus_status to TRUE for No Match Abstracts
                match = find_synonyms_in_text(synonym, text)
                focus_status = (match == {} or match == {''} or match == set())
                print(">>>>>>>>>>Focus",focus_status)
                metadata_json = extract_metadata(synonym, paper["title"], abstract_with_keywords, pubmed_type,focus_status)
                result = {"root_name": root_name, "synonyms" : synonym, "PMID": pmid, "pubmed_type": pubmed_type, "metadata": metadata_json}
            else:
                print("No abstract:",pmid)
                continue
        except Exception as e:
            result = {"root_name": root_name,  "synonyms" : synonym, "PMID": pmid, "pubmed_type": pubmed_type,"Error": str(e)}

        results.append(result)

        # Save incrementally after each PMID
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2)

    return results
