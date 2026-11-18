"""NCBI access: lineage/group lookup and host / isolation-source evidence."""
import json
import os
import re
import time

import requests

BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
GENOME_API = BASE + "/genome/taxon/{taxid}/dataset_report"
VIRUS_API = BASE + "/virus/taxon/{taxid}/dataset_report"
TAXONOMY_API = BASE + "/taxonomy/taxon/{taxid}/dataset_report"
DEFAULT_CACHE = os.path.join(os.path.expanduser("~"), ".cache", "pathoannot", "cache_v3.json")
MIN_ANNOTATED = 5

LINEAGE_GROUPS = [(10239, "Virus"), (2, "Bacteria"), (2157, "Archaea"), (4751, "Fungi"),
                  (2759, "Eukaryota (non-fungal)")]

HUMAN_PAT = re.compile(r"homo sapiens|\bhuman\b|\bpatient\b", re.I)
CLINICAL_PAT = re.compile(
    r"blood|csf|cerebrospinal|sputum|bronch|urine|wound|abscess|pus|swab|"
    r"stool|fec|rectal|throat|nasal|nasopharyn|clinical|tissue|biopsy|"
    r"respiratory|lung|skin|vagin|cervical|ascit|pleural|synovial|catheter", re.I)
ENV_PAT = re.compile(
    r"soil|sediment|seawater|marine|river|lake|freshwater|wastewater|sewage|"
    r"compost|rhizosph|plant|leaf|root|hot spring|ocean|air\b|dust|surface", re.I)
NON_HUMAN_HOST = re.compile(
    r"\b(bos|sus|gallus|mus|rattus|canis|felis|equus|ovis|capra|apis|bird|"
    r"cattle|pig|chicken|dog|cat|horse|sheep|goat|fish|mosquito|tick|bat)\b", re.I)


# ------------------------------------------------------------------ cache
def load_cache(path):
    if os.path.exists(path):
        try:
            with open(path) as fh:
                return json.load(fh)
        except (OSError, ValueError):
            return {}
    return {}


def save_cache(cache, path):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        json.dump(cache, fh)


def cached(cache, key, fn, cache_path=None):
    """Return cache[key], computing it with fn() if missing. Failures are not cached."""
    if key in cache and cache[key] is not None:
        return cache[key]
    try:
        value = fn()
    except Exception:
        return None
    cache[key] = value
    if cache_path:
        save_cache(cache, cache_path)
    return value


# ------------------------------------------------------------------ http
def _get_json(url, params, api_key):
    headers = {"Accept": "application/json"}
    if api_key:
        headers["api-key"] = api_key
    r = requests.get(url, params=params, headers=headers, timeout=60)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()


def get_lineage(taxid, api_key=None):
    """Return (group, [lineage TaxIDs]) for a taxon."""
    data = _get_json(TAXONOMY_API.format(taxid=taxid), {}, api_key)
    reports = (data or {}).get("reports") or []
    if not reports:
        return "Unknown", []
    tax = reports[0].get("taxonomy") or {}
    ids = set()
    for key in ("parents", "lineage"):
        for x in tax.get(key, []) or []:
            try:
                ids.add(int(x))
            except (TypeError, ValueError):
                pass
    try:
        ids.add(int(tax.get("tax_id", taxid)))
    except (TypeError, ValueError):
        pass
    group = "Unknown"
    for tid, name in LINEAGE_GROUPS:
        if tid in ids:
            group = name
            break
    return group, sorted(ids)


def get_group(taxid, api_key=None):
    return get_lineage(taxid, api_key)[0]


# ------------------------------------------------------------------ NCBI records
def parse_genome_report(rep):
    bs = (rep.get("assembly_info") or {}).get("biosample") or {}
    host = bs.get("host") or ""
    iso = bs.get("isolation_source") or ""
    for att in bs.get("attributes", []) or []:
        n = (att.get("name") or "").lower()
        if not host and n == "host":
            host = att.get("value", "")
        if not iso and n in ("isolation_source", "isolation source"):
            iso = att.get("value", "")
    return host, iso


def parse_virus_report(rep):
    h = rep.get("host") or {}
    host = h.get("sci_name") or h.get("organism_name") or h.get("name") or ""
    iso_d = rep.get("isolate") or {}
    iso = iso_d.get("source") or iso_d.get("isolation_source") or ""
    return host, iso


def fetch_records(url, parser, max_pages=2, api_key=None):
    rows, token, pages, total = [], None, 0, None
    while pages < max_pages:
        params = {"page_size": 500}
        if token:
            params["page_token"] = token
        data = _get_json(url, params, api_key)
        if data is None:
            break
        total = data.get("total_count", total)
        rows.extend(parser(rep) for rep in data.get("reports", []))
        token = data.get("next_page_token")
        pages += 1
        if not token:
            break
        time.sleep(0.4)
    return rows, (total if total is not None else len(rows))


def fetch_genome_reports(taxid, max_pages=2, api_key=None):
    return fetch_records(GENOME_API.format(taxid=taxid), parse_genome_report, max_pages, api_key)


def fetch_virus_reports(taxid, max_pages=2, api_key=None):
    return fetch_records(VIRUS_API.format(taxid=taxid), parse_virus_report, max_pages, api_key)


def summarise(rows):
    n = len(rows)
    if n == 0:
        return dict(n_genomes=0, n_with_host_info=0, pct_human_host=None,
                    pct_clinical_source=None, pct_environmental_source=None,
                    pct_other_animal_host=None)
    with_info = [(h, i) for h, i in rows if h or i]
    human = sum(1 for h, i in rows if HUMAN_PAT.search(h or ""))
    clin = sum(1 for h, i in rows if CLINICAL_PAT.search(i or ""))
    env = sum(1 for h, i in rows if ENV_PAT.search(i or "") and not HUMAN_PAT.search(h or ""))
    animal = sum(1 for h, i in rows if NON_HUMAN_HOST.search(h or "") and not HUMAN_PAT.search(h or ""))
    d = max(len(with_info), 1)
    return dict(n_genomes=n, n_with_host_info=len(with_info),
                pct_human_host=round(100 * human / d, 1),
                pct_clinical_source=round(100 * clin / d, 1),
                pct_environmental_source=round(100 * env / d, 1),
                pct_other_animal_host=round(100 * animal / d, 1))


def classify(s, min_annotated=MIN_ANNOTATED):
    """Heuristic NCBI-based label. Evidence only, not a diagnosis."""
    if s["n_genomes"] == 0:
        return "No genome records"
    if s["n_with_host_info"] < min_annotated:
        return "Insufficient metadata"
    h, c, e = s["pct_human_host"], s["pct_clinical_source"], s["pct_environmental_source"]
    if h >= 50 or c >= 40:
        return "Strong human/clinical association"
    if h >= 15 or c >= 15:
        return "Moderate human association"
    if e >= 50 and h < 5:
        return "Mostly environmental"
    return "Weak/unclear association"


def lookup_ncbi(taxid, group="Unknown", api_key=None):
    """Host/isolation evidence, using NCBI Virus for viruses and genome records otherwise."""
    if group == "Virus":
        order = ["virus"]
    elif group != "Unknown":
        order = ["genome"]
    else:
        order = ["genome", "virus"]
    best = None
    for src in order:
        fetch = fetch_virus_reports if src == "virus" else fetch_genome_reports
        rows, total = fetch(taxid, api_key=api_key)
        cand = dict(record_source=src, n_records_total=total, **summarise(rows))
        if best is None or cand["n_with_host_info"] > best["n_with_host_info"]:
            best = cand
        if best["n_with_host_info"] >= MIN_ANNOTATED:
            break
    return best
