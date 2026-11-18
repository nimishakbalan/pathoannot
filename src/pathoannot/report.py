"""Reporting tiers and the HTML report."""
import datetime
import html

import pandas as pd

from . import __version__

DEFAULT_MIN_RPM = 10.0
DEFAULT_MIN_READS = 50


def _num(v, default=0.0):
    try:
        return default if pd.isna(v) else float(v)
    except (TypeError, ValueError):
        return default


def assign_report_tier(df, min_rpm=DEFAULT_MIN_RPM, min_reads=DEFAULT_MIN_READS):
    """Reportable / Review / Excluded / Not reported, with a plain-language reason."""
    tiers, reasons = [], []
    for _, r in df.iterrows():
        flag = str(r.get("contamination_flag", ""))
        cls = r.get("pathogen_class", "Unresolved")
        reads = _num(r.get("new_est_reads"))
        rpm = _num(r.get("rpm", r.get("reads_per_million")))
        tier = str(r.get("presence_tier", "Not assessed"))
        if tier == "nan":
            tier = "Not assessed"
        why = []
        if flag.startswith("Likely contaminant"):
            tiers.append("Excluded: likely contaminant")
            reasons.append(f"{flag}; RPM ratio vs NTC {r.get('rpm_ratio_vs_ntc')}")
            continue
        if cls == "Pathogen/opportunist":
            ok_abund = reads >= min_reads and rpm >= min_rpm
            if not ok_abund:
                why.append(f"below abundance thresholds (reads >= {min_reads:g}, RPM >= {min_rpm:g})")
            if tier not in ("Confirmed", "Probable"):
                why.append(f"alignment check: {tier}")
            if "contaminant" in flag.lower():
                why.append(flag)
            if why:
                tiers.append("Review")
                reasons.append("; ".join(why))
            else:
                tiers.append("Reportable")
                reasons.append("pathogen evidence + abundance + alignment confirmation + passes contamination checks")
        elif cls == "Unresolved" and reads >= min_reads and tier == "Confirmed":
            tiers.append("Review")
            reasons.append("class unresolved but organism confirmed present")
        else:
            tiers.append("Not reported")
            reasons.append(f"class: {cls}")
    return df.assign(report_tier=tiers, report_reason=reasons)


# ------------------------------------------------------------------ HTML
CSS = """
:root{--bg:#fff;--fg:#1d2733;--muted:#5b6876;--line:#d8dee6;--head:#eef2f6;--ok:#1b7a3d;--warn:#a56300;--bad:#b3261e}
@media (prefers-color-scheme:dark){:root{--bg:#161b22;--fg:#e6edf3;--muted:#9aa7b4;--line:#30363d;--head:#1f262e;--ok:#56d364;--warn:#e3b341;--bad:#ff7b72}}
body{font:15px/1.5 system-ui,sans-serif;background:var(--bg);color:var(--fg);max-width:1100px;margin:2rem auto;padding:0 1rem}
h1{font-size:1.6rem;margin-bottom:.2rem}h2{font-size:1.15rem;margin-top:2rem;border-bottom:1px solid var(--line);padding-bottom:.3rem}
.muted{color:var(--muted)}table{border-collapse:collapse;width:100%;margin:.6rem 0;font-size:13px}
th,td{border:1px solid var(--line);padding:.35rem .5rem;text-align:left;vertical-align:top}th{background:var(--head)}
.wrap{overflow-x:auto}.ok{color:var(--ok);font-weight:600}.warn{color:var(--warn);font-weight:600}.bad{color:var(--bad);font-weight:600}
.box{border:1px solid var(--line);padding:.7rem 1rem;border-radius:6px;background:var(--head)}
.sign td{height:2.2rem}@media print{body{max-width:none;margin:0}}
"""

COLS_MAIN = [("species", "Organism"), ("group", "Group"), ("pathogen_class", "Class"), ("new_est_reads", "Reads"),
             ("rpm", "RPM"), ("mean_depth", "Depth"), ("breadth_pct", "Breadth %"), ("breadth_ratio", "Breadth ratio"),
             ("presence_tier", "Alignment check"), ("deepac_mean_score", "DeePaC mean"),
             ("deepac_pct_above_threshold", "DeePaC % reads >=0.5"), ("priority_hits", "Priority lists"),
             ("lit_label", "Literature"), ("ncbi_label", "NCBI host data"), ("contamination_flag", "Contamination"),
             ("amr_genes", "Resistance genes"), ("report_reason", "Why")]


def _fmt(v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    if isinstance(v, float):
        return f"{v:.3g}" if abs(v) < 1000 else f"{v:,.0f}"
    return html.escape(str(v))


def _table(df, cols):
    cols = [(c, h) for c, h in cols if c in df.columns]
    if df.empty or not cols:
        return '<p class="muted">None.</p>'
    head = "".join(f"<th>{html.escape(h)}</th>" for _, h in cols)
    body = "".join("<tr>" + "".join(f"<td>{_fmt(r[c])}</td>" for c, _ in cols) + "</tr>" for _, r in df.iterrows())
    return f'<div class="wrap"><table><tr>{head}</tr>{body}</table></div>'


def build_html(sample, df, amr_df=None, meta=None, params=None):
    meta, params = meta or {}, params or {}
    tier = df["report_tier"]
    rep, rev = df[tier == "Reportable"], df[tier == "Review"]
    exc, other = df[tier.str.startswith("Excluded")], df[tier == "Not reported"]
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    meta_rows = "".join(f"<tr><th>{html.escape(k)}</th><td>{_fmt(v)}</td></tr>" for k, v in meta.items())
    prm_rows = "".join(f"<tr><th>{html.escape(k)}</th><td>{_fmt(v)}</td></tr>" for k, v in params.items())
    amr_html = ('<p class="muted">Resistance analysis not run, or no confident hits.</p>' if amr_df is None or amr_df.empty
                else _table(amr_df, [("species", "Attributed species"), ("gene", "Gene"), ("class", "Class"),
                                     ("subclass", "Subclass"), ("identity_pct", "Identity %"),
                                     ("coverage_pct", "Coverage %"), ("contig", "Contig")]))
    counts = (f'<span class="ok">{len(rep)} reportable</span> &middot; <span class="warn">{len(rev)} need review</span>'
              f' &middot; <span class="bad">{len(exc)} excluded as likely contaminants</span> &middot; {len(other)} not reported')
    return f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>pathoannot report: {html.escape(sample)}</title>
<style>{CSS}</style></head><body>
<h1>Metagenomic pathogen report: {html.escape(sample)}</h1>
<p class="muted">Generated {now} by pathoannot {__version__}. Research use: not a validated diagnostic report.</p>
<div class="box"><b>Summary:</b> {counts}</div>
<h2>Run information</h2><table>{meta_rows or '<tr><td class="muted">No metadata supplied.</td></tr>'}</table>
<h2>Organisms meeting reporting criteria</h2>{_table(rep, COLS_MAIN)}
<h2>Organisms needing expert review</h2>{_table(rev, COLS_MAIN)}
<h2>Resistance genes (AMRFinderPlus on assembled contigs)</h2>{amr_html}
<p class="muted">Genes are attributed to a species from the Kraken2 classification of their contig. Plasmid-borne genes
can sit on contigs that are unassigned or assigned to a different species, so attribution is indicative only.</p>
<h2>Excluded as likely contaminants</h2>{_table(exc, COLS_MAIN)}
<h2>Other species (commensal, environmental or low evidence; not reported)</h2>{_table(other.head(40), COLS_MAIN)}
<h2>Methods</h2>
<p>Reads were classified with Kraken2 and abundances re-estimated with Bracken. Each species was assigned a class
(pathogen/opportunist, commensal, environmental, unresolved) from independent evidence: WHO priority pathogen lists
and any user-supplied list (for example ICMR), the Europe PMC literature, and NCBI host and isolation-source metadata.
Contamination was assessed against a negative control (RPM ratio) and a list of common kit-contaminant genera.
Presence was checked by extracting each species' reads and aligning them to a reference genome (minimap2); depth,
breadth and the ratio of observed to expected breadth were computed with samtools. Read-level pathogenic potential
was scored with DeePaC (bacteria). Resistance genes were found with AMRFinderPlus on assembled contigs.</p>
<h2>Thresholds and parameters</h2><table>{prm_rows}</table>
<h2>Limitations</h2><ul>
<li>Thresholds are defaults and have not been validated; tune them with known-positive and negative samples.</li>
<li>WHO lists describe drug-resistant pathogens; a list match shows concern about the species, not resistance in this sample.</li>
<li>Literature and NCBI metadata are biased toward disease isolates; they indicate association, not causation.</li>
<li>DeePaC scores are read-level and noisy; they support, and never replace, the other evidence.</li>
<li>Detection of a pathogen does not prove infection. Interpret with clinical findings and negative controls.</li>
<li>Clinical use requires full validation (for example CLIA or ISO 15189) and sign-off by a qualified person.</li></ul>
<h2>Review and sign-off</h2><table class="sign"><tr><th>Reviewed by</th><td></td><th>Date</th><td></td></tr>
<tr><th>Interpretation</th><td colspan="3"></td></tr><tr><th>Signature</th><td colspan="3"></td></tr></table>
</body></html>"""
