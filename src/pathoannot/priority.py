"""Priority-pathogen lists (WHO built in; ICMR or any other list supplied by the user)."""
import re
from importlib import resources

import pandas as pd

BUILTIN_FILES = ["who_bppl_2024.csv", "who_fppl_2022.csv"]
RANK = {"critical": 4, "high": 3, "medium": 2}


def norm(name):
    """Lower-case, strip [brackets] used by NCBI for unplaced genera, collapse spaces."""
    return re.sub(r"\s+", " ", re.sub(r"[\[\]]", "", str(name))).strip().lower()


def _rows(df):
    out = []
    for _, r in df.iterrows():
        if pd.isna(r.get("taxon")) or not str(r["taxon"]).strip():
            continue
        names = [str(r["taxon"])]
        if pd.notna(r.get("aliases")):
            names += [a for a in str(r["aliases"]).split(";") if a.strip()]
        tid = r.get("taxid")
        out.append(dict(
            list=str(r.get("list", "user list")) if pd.notna(r.get("list")) else "user list",
            priority=str(r.get("priority", "listed")) if pd.notna(r.get("priority")) else "listed",
            names=[norm(n) for n in names],
            taxid=int(tid) if pd.notna(tid) and str(tid).strip() else None,
            note=str(r["note"]) if pd.notna(r.get("note")) else "",
        ))
    return out


def load_builtin():
    entries = []
    for fn in BUILTIN_FILES:
        with resources.files("pathoannot").joinpath(f"data/priority/{fn}").open("r", encoding="utf-8") as fh:
            entries += _rows(pd.read_csv(fh))
    return entries


def load_user(path):
    df = pd.read_csv(path)
    if "taxon" not in df.columns and "taxid" not in df.columns:
        raise ValueError(f"{path}: needs a 'taxon' and/or 'taxid' column")
    if "taxon" not in df.columns:
        df["taxon"] = ""
    entries = []
    for e in _rows(df.assign(taxon=df["taxon"].replace("", pd.NA).fillna("_taxid_only_"))):
        entries.append(e)
    return entries


def match(species, taxid, lineage_ids, entries):
    """Return the list entries that apply to this species (by name, genus prefix, or TaxID/clade)."""
    sp = norm(species)
    lin = set(lineage_ids or [])
    if taxid is not None:
        lin.add(int(taxid))
    hits = []
    for e in entries:
        by_name = any(n != "_taxid_only_" and (sp == n or sp.startswith(n + " ")) for n in e["names"])
        by_id = e["taxid"] is not None and e["taxid"] in lin
        if by_name or by_id:
            hits.append(e)
    return hits


def summarise_hits(hits):
    if not hits:
        return "", ""
    best = max(hits, key=lambda e: RANK.get(e["priority"].lower(), 1))
    text = "; ".join(f"{e['list']} [{e['priority']}]" + (f" ({e['note']})" if e["note"] else "") for e in hits)
    return best["priority"], text
