"""Step 1 (annotate): combine priority lists, literature and NCBI evidence per species."""
import os
import sys

import pandas as pd

from . import contam, core, literature, priority

ALL_SOURCES = ("priority", "literature", "ncbi")
CLASS_OF_CALL = {
    "Likely pathogen/opportunist": "Pathogen/opportunist",
    "Possible pathogen/opportunist": "Pathogen/opportunist",
    "Likely commensal / human-associated": "Commensal",
    "Likely environmental": "Environmental",
}


def read_bracken(path):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"taxonomy_id": "TaxID", "name": "species"})
    missing = {"species", "TaxID", "new_est_reads", "fraction_total_reads"} - set(df.columns)
    if missing:
        raise ValueError(f"{path}: not a Bracken output file (missing {sorted(missing)})")
    return df[["species", "TaxID", "new_est_reads", "fraction_total_reads"]]


def sample_name(path):
    s = os.path.basename(path)
    for ext in (".bracken", ".txt", ".tsv"):
        if s.endswith(ext):
            s = s[: -len(ext)]
    return s


def final_call(support, ncbi_label, lit_label, has_any):
    """Combine independent lines of evidence into one cautious call."""
    if support >= 2:
        return "Likely pathogen/opportunist", ""
    if support >= 1:
        return "Possible pathogen/opportunist", "REVIEW"
    if ncbi_label == "Moderate human association" or lit_label == "Commensal/microbiome-dominated literature":
        return "Likely commensal / human-associated", "REVIEW"
    if ncbi_label == "Mostly environmental" or lit_label == "Environmental-dominated literature":
        return "Likely environmental", ""
    if not has_any:
        return "Unresolved (insufficient evidence)", "REVIEW"
    return "Low pathogen evidence", "REVIEW"


def annotate_sample(path, min_reads=10, sources=ALL_SOURCES, cache=None, cache_path=core.DEFAULT_CACHE,
                    api_key=None, priority_entries=None, min_annotated=core.MIN_ANNOTATED, log=sys.stderr,
                    ntc_paths=(), contaminants=None, ntc_ratio=contam.DEFAULT_RATIO):
    cache = cache if cache is not None else {}
    priority_entries = priority_entries if priority_entries is not None else priority.load_builtin()
    full = read_bracken(path)
    total_reads = float(full["new_est_reads"].sum()) or 1.0
    df = full[full["new_est_reads"] >= min_reads].copy()
    df["TaxID"] = df["TaxID"].astype("Int64")
    name = sample_name(path)
    print(f"[{name}] annotating {len(df)} species (sources: {', '.join(sources)})", file=log)

    rows = []
    for _, r in df.iterrows():
        tid = int(r["TaxID"]) if pd.notna(r["TaxID"]) else None
        rec = dict(sample=name, species=r["species"], TaxID=tid, new_est_reads=int(r["new_est_reads"]),
                   fraction_total_reads=r["fraction_total_reads"],
                   reads_per_million=round(1e6 * r["new_est_reads"] / total_reads, 1))
        lin = core.cached(cache, f"tax:{tid}", lambda: core.get_lineage(tid, api_key), cache_path) if tid else None
        group, lineage = (lin if lin else ("Unknown", []))
        rec["group"] = group

        # 1. priority lists
        if "priority" in sources:
            lvl, text = priority.summarise_hits(priority.match(r["species"], tid, lineage, priority_entries))
            rec["priority_level"], rec["priority_hits"] = lvl, text
        # 2. NCBI host / isolation source
        ncbi_label = ""
        if "ncbi" in sources and tid:
            s = core.cached(cache, f"ncbi:{tid}", lambda: core.lookup_ncbi(tid, group, api_key), cache_path)
            if s:
                ncbi_label = core.classify(s, min_annotated)
                rec.update(ncbi_label=ncbi_label, **{k: s[k] for k in (
                    "record_source", "n_records_total", "pct_human_host", "pct_clinical_source",
                    "pct_environmental_source")})
            else:
                rec["ncbi_label"] = "Query failed"
        # 3. literature
        lit_lab = ""
        if "literature" in sources:
            ls = core.cached(cache, f"lit:{r['species']}", lambda: literature.lit_summary(r["species"]), cache_path)
            if ls:
                lit_lab = literature.lit_label(ls)
                rec.update(ls, lit_label=lit_lab)
            else:
                rec["lit_label"] = "Query failed"

        support, lines = 0.0, []
        if rec.get("priority_level"):
            support += 1
            lines.append("priority list")
        if ncbi_label == "Strong human/clinical association":
            support += 1
            lines.append("NCBI host data")
        elif ncbi_label == "Moderate human association":
            support += 0.5
        if lit_lab == "Disease-dominated literature":
            support += 1
            lines.append("literature")
        has_any = bool(rec.get("priority_level")) or ncbi_label not in ("", "Query failed", "Insufficient metadata", "No genome records") \
            or lit_lab not in ("", "Query failed", "Insufficient literature")
        rec["pathogen_support"] = support
        rec["evidence_lines"] = "; ".join(lines)
        rec["final_call"], rec["review_flag"] = final_call(support, ncbi_label, lit_lab, has_any)
        rows.append(rec)

    out = pd.DataFrame(rows)
    out = contam.apply(out, ntc_paths, contaminants, ntc_ratio)
    out["pathogen_class"] = out["final_call"].map(CLASS_OF_CALL).fillna("Unresolved")
    return out.sort_values("new_est_reads", ascending=False).reset_index(drop=True)
