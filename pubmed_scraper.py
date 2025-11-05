"""
PubMed Ingredient Scraper

This file:
  1) Reads an Excel file of ingredients (default: "76 Ingredient List.xlsx") from *input_dir*.
  2) Builds an OR query per ingredient across its synonyms (Common, Scientific, Other).
  3) Queries PubMed (Entrez) restricted to allowed study types.
  4) Appends per-hit rows to a raw CSV in *progress_dir*, deduping by PMID across runs.
  5) Writes an audit JSONL in *progress_dir* with per-root run stats:
     {"root","ts_utc","found_in_query","new_appended","total_known_after"}.
  6) Produces a de-duplicated aggregated Excel in *output_dir*.
"""

# =========================
# = GLOBAL CONFIGURATION  =
# =========================

import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from dotenv import load_dotenv
import pandas as pd
from Bio import Entrez
from tqdm.auto import tqdm
import os
load_dotenv()
# ---------- Entrez Credentials ----------
ENTREZ_EMAIL =  os.environ.get("ENTREZ_EMAIL", "")
ENTREZ_API_KEY = os.environ.get("ENTREZ_API_KEY", None)

# ---------- Filenames (relative to the respective dirs) ----------
INGREDIENT_XLSX_NAME = "Ingredients.xlsx"
RAW_CSV_NAME         = "pubmed_hits_raw.csv"
PROGRESS_JL_NAME     = "pubmed_progress.jl"   # JSON Lines progress
FINAL_XLSX_NAME      = "pubmed_ingredients.xlsx"

# ---------- Query / Fetch Parameters ----------
BATCH_SIZE       = 200
SLEEP_SECS       = 0.5
SLICE_THRESHOLD  = 10_000

# ---------- Allowed Study Types ----------
SYSTEMATIC_META_ANALYSIS = [
    "Meta-Analysis",
    "Network Meta-Analysis",
    "Systematic Review",
]

CLINICAL_TRIAL = [
    "Adaptive Clinical Trial",
    "Clinical Study",
    "Clinical Trial",
    "Clinical Trial Protocol",
    "Clinical Trial, Phase I",
    "Clinical Trial, Phase II",
    "Clinical Trial, Phase III",
    "Clinical Trial, Phase IV",
    "Clinical Trial, Veterinary",
    "Controlled Clinical Trial",
    "Equivalence Trial",
    "Multicenter Study",
    "Pragmatic Clinical Trial",
    "Randomized Controlled Trial",
    "Randomized Controlled Trial, Veterinary",
]

OBSERVATIONAL_STUDIES = [
    "Case Reports",
    "Comparative Study",
    "Evaluation Study",
    "Observational Study",
    "Observational Study, Veterinary",
    "Twin Study",
    "Validation Study",
    "Review",
    "Scientific Integrity Review",
    "Scoping Review",
]

ALLOWED_STUDY_TYPES = (
    SYSTEMATIC_META_ANALYSIS + CLINICAL_TRIAL + OBSERVATIONAL_STUDIES
)

# =========================
# = Helper Functions      =
# =========================

def _require_entrez_credentials():
    """
    Ensure Entrez credentials are available and configure the Bio.
    Entrez client accordingly. Raises a RuntimeError if ENTREZ_EMAIL is missing.
    """
    if not ENTREZ_EMAIL:
        raise RuntimeError(
            "ENTREZ_EMAIL is not set. Set env var ENTREZ_EMAIL or the module constant."
        )
    Entrez.email = ENTREZ_EMAIL
    if ENTREZ_API_KEY:
        Entrez.api_key = ENTREZ_API_KEY

def _normalize_year(year_value):
    """
    Return a 4-digit publication year extracted from a variety of PubMed 
    date fields, or 'NA' if not found. Accepts strings or numbers.
    """
    if year_value is None:
        return "NA"
    s = str(year_value).strip()
    if not s:
        return "NA"
    m = re.search(r"\b(\d{4})\b", s)
    return m.group(1) if m else "NA"

def _split_general(s):
    """
    Split common/other search terms by commas (outside parentheses), semicolons, 
    newlines, slashes, or pipes. Returns a cleaned list of unique-looking tokens.
    """
    if not isinstance(s, str) or not s.strip():
        return []
    parts = re.split(r'(?:,(?![^()]*\))|;|\n|/|\|)', s)
    return [p.strip() for p in parts if p and p.strip()]

def _split_semicolon_only(s):
    """
    Split scientific names only on semicolons or newlines, preserving 
    commas within names. Returns a list of trimmed parts.
    """
    if not isinstance(s, str) or not s.strip():
        return []
    parts = re.split(r'(?:;|\n)', s)
    return [p.strip() for p in parts if p and p.strip()]

def _build_type_clause(allowed_types):
    """
    Build the PubMed publication-type filter clause 
    (e.g., '"Clinical Trial"[ptyp] OR "Review"[ptyp]').
    """
    return " OR ".join([f'"{t}"[ptyp]' for t in allowed_types if t])

def _build_term_clause(term_or_list):
    """
    Build the search term clause, accepting either a single term or 
    a list; lists are OR-joined. Terms are searched in [tiab].
    """
    if isinstance(term_or_list, (list, tuple)):
        return " OR ".join([f'"{t}"[tiab]' for t in term_or_list if t])
    return f'"{term_or_list}"[tiab]'

def build_query(term_or_list, allowed_types):
    """
    Compose the final PubMed query by AND-ing the 
    term clause with the publication-type clause.
    """
    return f"({ _build_term_clause(term_or_list) }) AND ({ _build_type_clause(allowed_types) })"

def _count_hits(query: str) -> int:
    """Return PubMed esearch 'Count' for an arbitrary query string."""
    h = Entrez.esearch(db="pubmed", term=query, retmax=0)
    r = Entrez.read(h)
    return int(r.get("Count", 0))

def iterate_pubmed_in_slices(query, batch_size = BATCH_SIZE, years_step = 5):
    """
    Generator that yields lists of PubMed records by iterating over year slices 
    to avoid huge result sets. Handles simple retry and pacing between requests.
    """
    now_year = datetime.now().year
    slices = [(y, min(y + years_step - 1, now_year)) for y in range(1940, now_year + 1, years_step)]
    for y0, y1 in slices:
        sh = Entrez.esearch(db="pubmed", term=query, usehistory="y", retmax=0,
                            datetype="pdat", mindate=str(y0), maxdate=str(y1))
        sr = Entrez.read(sh)
        count = int(sr.get("Count", 0))
        if count == 0:
            continue
        webenv, qk = sr["WebEnv"], sr["QueryKey"]
        for start in range(0, count, batch_size):
            last_err = None
            for _ in range(3):
                try:
                    fh = Entrez.efetch(db="pubmed", retmode="xml", retstart=start, retmax=batch_size, webenv=webenv, query_key=qk)
                    rec = Entrez.read(fh)
                    yield rec.get("PubmedArticle", []) + rec.get("PubmedBookArticle", [])
                    break
                except Exception as e:
                    last_err = e
                    time.sleep(1.5)
            else:
                raise last_err
            time.sleep(SLEEP_SECS)

# =========================
# = Progress (JSONL)      =
# =========================

def _append_audit_log(progress_path, root, found_in_query, new_appended, total_known_after):
    """Append an audit record to JSONL showing per-root run stats (timestamped UTC)."""
    rec = {
        "root": root,
        "ts_utc": datetime.now(timezone.utc).isoformat(),
        "found_in_query": int(found_in_query),
        "new_appended": int(new_appended),
        "total_known_after": int(total_known_after),
    }
    with progress_path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")

def _load_progress_jl(path):
    """
    Read a JSONL progress file and return a set of completed root ingredient names. 
    Skips empty or malformed lines without failing the run.
    """
    done = set()
    if not path.exists():
        return done
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
                root = str(rec.get("root", "")).strip()
                if root:
                    done.add(root)
            except Exception:
                # Skip malformed lines but continue
                continue
    return done

def _append_progress_jl(path, root):
    """
    Append a single JSONL record like {"root": "<ingredient>"} to mark a completed root. 
    Creates the file if it does not exist.
    """
    rec = {"root": root}
    with path.open("a", encoding="utf-8") as f:
        f.write(json.dumps(rec, ensure_ascii=False) + "\n")

# =========================
# = Core Logic            =
# =========================

def _build_synonym_map(tbl):
    """
    Construct a mapping {root -> [root + unique synonyms]} using the 
    input DataFrame columns. Deduplicates while preserving order.
    """
    out = {}
    for _, row in tbl.iterrows():
        root = str(row.get("Ingredient", "")).strip()
        if not root:
            continue
        terms = [root]
        terms += _split_general(row.get("Common Name(s)", ""))
        terms += _split_general(row.get("Scientific Name(s)", ""))
        terms += _split_general(row.get("Other Search Words", ""))

        seen = set()
        uniq = []
        for t in terms:
            if t and t not in seen:
                seen.add(t)
                uniq.append(t)
        out[root] = uniq
    return out

def _ensure_raw_csv(path):
    """
    Create the append-safe raw CSV with headers if it does not already exist. 
    This avoids header duplication on subsequent appends.
    """
    if not path.exists():
        cols = ["Root Name","Search Term","PMID","Title","Article Type","Publication Year","PubMed URL"]
        pd.DataFrame(columns=cols).to_csv(path, index=False)

def _process_papers(papers, root, term_repr):
    """
    Transform a batch of PubMed XML records into normalized row dicts ready for CSV append. 
    Extracts PMID, title, types, year, and URL fields.
    """
    rows = []
    for paper in papers:
        try:
            mc = paper["MedlineCitation"]
            pmid = str(mc["PMID"])
            art  = mc.get("Article", {})
            title = str(art.get("ArticleTitle", ""))

            pt_list = art.get("PublicationTypeList", [])
            pub_types = [str(x) for x in pt_list] if pt_list else ["NA"]

            pub_year = None
            try:
                ad = art.get("ArticleDate", [])
                if ad:
                    pub_year = ad[0].get("Year")
                if not pub_year:
                    jd = art["Journal"]["JournalIssue"]["PubDate"]
                    pub_year = jd.get("Year") or jd.get("MedlineDate")
            except Exception:
                pub_year = "NA"
            pub_year = _normalize_year(pub_year)

            rows.append({
                "Root Name": root,
                "Search Term": term_repr,
                "PMID": pmid,
                "Title": title,
                "Article Type": "; ".join(pub_types),
                "Publication Year": pub_year,
                "PubMed URL": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            })
        except Exception:
            continue
    return rows

def _process_combined_query(root, terms):
    """
    Run one combined OR PubMed query for all terms of a given root and return 
    per-paper result rows and total new hits. Uses year slicing for very large result sets.
    """
    query = build_query(list(terms), ALLOWED_STUDY_TYPES)
    print("🔵Query",query)
    probe = Entrez.esearch(db="pubmed", term=query, retmax=0, usehistory="y")
    pr = Entrez.read(probe)
    total_hits = int(pr.get("Count", 0))
    if total_hits == 0:
        return [], 0

    out_rows = []
    term_repr = "; ".join(terms)

    if total_hits <= SLICE_THRESHOLD:
        webenv, qk = pr["WebEnv"], pr["QueryKey"]
        with tqdm(total=total_hits, desc=f"{root} | combined", unit="papers", leave=False) as pbar:
            for start in range(0, total_hits, BATCH_SIZE):
                last_err = None
                for _ in range(3):
                    try:
                        fh = Entrez.efetch(db="pubmed", retmode="xml", retstart=start, retmax=BATCH_SIZE, webenv=webenv, query_key=qk)
                        rec = Entrez.read(fh)
                        papers = rec.get("PubmedArticle", []) + rec.get("PubmedBookArticle", [])
                        out_rows.extend(_process_papers(papers, root, term_repr))
                        break
                    except Exception as e:
                        last_err = e
                        time.sleep(1.5)
                else:
                    raise last_err
                pbar.update(len(papers))
                time.sleep(SLEEP_SECS)
    else:
        with tqdm(desc=f"{root} | combined (sliced)", unit="papers", leave=False) as pbar:
            for papers in iterate_pubmed_in_slices(query, batch_size=BATCH_SIZE, years_step=5):
                got = _process_papers(papers, root, term_repr)
                out_rows.extend(got)
                pbar.update(len(got))
                time.sleep(SLEEP_SECS)

    return out_rows, total_hits

def _aggregate_final(raw_csv):
    """
    Aggregate raw per-hit rows into a de-duplicated per-(Root, PMID) table 
    with merged publication types and search terms. Returns a DataFrame.
    """
    df = pd.read_csv(raw_csv, dtype=str, encoding="utf-8-sig", keep_default_na=False)
    if df.empty:
        return df
    agg = (
        df.groupby(["Root Name","PMID"], as_index=False)
          .agg({
              "Article Type": lambda s: "; ".join(sorted(set("; ".join(s).split("; ")))) if len(s) else "",
              "Publication Year": "first",
              "PubMed URL": "first",
          })
    )
    return agg

def _load_seen_pmids_by_root(raw_csv):
    """Fast-load existing (Root, PMID) map from the raw CSV to avoid duplicates on reruns."""
    seen = {}
    if not raw_csv.exists():
        return seen
    try:
        usecols = ["Root Name", "PMID"]
        df = pd.read_csv(raw_csv, usecols=usecols, dtype=str, keep_default_na=False)
        for root, group in df.groupby("Root Name"):
            seen[root] = set(group["PMID"].astype(str))
    except Exception:
        # Fallback to empty map on read issues
        seen = {}
    return seen

# =========================
# = Entry Function        =
# =========================

def pubmed_ingredient_search(input_dir, output_dir, progress_dir):
    """
    End-to-end runner: reads inputs, queries PubMed per root, appends raw rows, 
    tracks JSONL progress, and writes the final aggregated Excel output.
    """
    print(f"🚀 Beginning PubMed scraper...")
    
    _require_entrez_credentials()

    input_dir = Path(input_dir).expanduser().resolve()
    output_dir = Path(output_dir).expanduser().resolve()
    progress_dir = Path(progress_dir).expanduser().resolve()

    output_dir.mkdir(parents=True, exist_ok=True)
    progress_dir.mkdir(parents=True, exist_ok=True)

    ingredient_path = input_dir / INGREDIENT_XLSX_NAME
    raw_csv_path    = progress_dir / RAW_CSV_NAME
    progress_jl     = progress_dir / PROGRESS_JL_NAME
    final_xlsx_path = output_dir / FINAL_XLSX_NAME

    if not ingredient_path.exists():
        raise FileNotFoundError(f"Missing input Excel: {ingredient_path}")
    tbl = pd.read_excel(ingredient_path, dtype=str, keep_default_na=False)

    if "Ingredient" not in tbl.columns:
        raise ValueError("Input Excel must contain an 'Ingredient' column.")
    ingredients_list = tbl["Ingredient"].dropna().astype(str).str.strip().tolist()

    syn_map = _build_synonym_map(tbl)
    roots = [r for r in ingredients_list if r in syn_map]

    _ensure_raw_csv(raw_csv_path)
    
    # Load already-seen PMIDs by root to support dedup + "new this run" counts
    seen_by_root = _load_seen_pmids_by_root(raw_csv_path)

    total_known_global = sum(len(s) for s in seen_by_root.values())
    if total_known_global > 0:
        print(f"⏯️ Resuming: {total_known_global} known PMIDs across {len(seen_by_root)} roots from previous runs.")
    else:
        print("🆕 Starting fresh: no prior records found in raw CSV.")

    total_new_all_roots = 0

    for idx, root in enumerate(roots, 1):
        pre_seen = seen_by_root.get(root, set())
    
        terms = syn_map[root]
        print(f"[{idx}/{len(roots)}] ▶️ {root}: known_before={len(pre_seen)}")

        query_unfiltered = f"({_build_term_clause(terms)})"
        query_filtered   = build_query(terms, ALLOWED_STUDY_TYPES)

        # Count hits for each and print
        total_hits_unfiltered = _count_hits(query_unfiltered)
        total_hits_filtered   = _count_hits(query_filtered)

        if total_hits_filtered > SLICE_THRESHOLD:
            print(f"[{idx}/{len(roots)}] 🚩🚩🚩 {root}: filtered hits {total_hits_filtered:,} exceed threshold ({SLICE_THRESHOLD:,}); year-sliced fetching will be used.")

        print(f"[{idx}/{len(roots)}] {root}: hits_unfiltered={total_hits_unfiltered}, hits_filtered={total_hits_filtered}")

        try:
            rows, total_hits_est = _process_combined_query(root, syn_map[root])
            
            processed_count = len(rows)
            print(f"[{idx}/{len(roots)}] {root}: processed={processed_count} record(s) from API")

        except Exception as e:
            print(f"[{idx}/{len(roots)}] {root}: ERROR during query -> {e}")
            # Log a failed attempt with zeros (optional)
            _append_audit_log(progress_jl, root=root, found_in_query=0, new_appended=0, total_known_after=len(pre_seen))
            continue

        # Filter to only *new* PMIDs for append
        new_rows = []
        new_pmids = []
        for r in rows:
            pmid = r.get("PMID")
            if pmid and pmid not in pre_seen:
                new_rows.append(r)
                new_pmids.append(pmid)

        # Append only new rows
        if new_rows:
            pd.DataFrame.from_records(new_rows).to_csv(raw_csv_path, mode="a", header=False, index=False)
            # Update in-memory set
            pre_seen.update(new_pmids)
            seen_by_root[root] = pre_seen

        new_count = len(new_rows)
        total_known_after = len(pre_seen)
        total_new_all_roots += new_count

        # Print concise per-root summary
        print(f"[{idx}/{len(roots)}] ✅ {root}: hits_unfiltered={total_hits_unfiltered}, hits_filtered={total_hits_filtered}, new_this_run={len(new_rows)}, total_known={total_known_after}")
        
        # Write audit log entry
        _append_audit_log(progress_jl, root=root, found_in_query=total_hits_est, new_appended=new_count, total_known_after=total_known_after)

    # Aggregate final output
    agg = _aggregate_final(raw_csv_path)
    if not agg.empty:
        agg.to_excel(final_xlsx_path, index=False)

    # Final summary
    print("\n✅ Periodic run complete.")
    print(f"New rows appended this run (all roots): {total_new_all_roots}")
    print(f"Raw per-hit file: {raw_csv_path}")
    print(f"Audit log (JSONL): {progress_jl}")
    
    if final_xlsx_path.exists():
        print(f"Final Excel: {final_xlsx_path} | rows: {len(agg)}")
    else:
        print("Final Excel: (no rows) -> not written")
    
    print(f"🏁 Completed PubMed scraper.")


# pubmed_ingredient_search(
#     "ingredients_scraper/Data",
#     "ingredients_scraper/Progress",
#     "ingredients_scraper/Output"
# )
