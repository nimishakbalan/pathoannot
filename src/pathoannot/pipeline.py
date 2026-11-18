"""Orchestration of all steps for one sample."""
import os
import sys

import pandas as pd

from . import amr, evidence, reads, report, score
from . import contam as contam_mod

ALL_STEPS = ("score", "amr", "report")


def _log(msg):
    print(msg, file=sys.stderr)


def process_sample(spec, o, cache, entries, contaminants):
    sample = spec["sample"]
    steps = set(o.steps.split(","))
    ntc = tuple(o.ntc_bracken or ())

    # ---- 1. evidence + contamination (needs only the Bracken file)
    annot = evidence.annotate_sample(
        spec["bracken"], min_reads=o.min_reads, sources=tuple(o.sources.split(",")), cache=cache,
        cache_path=o.cache, api_key=o.api_key, priority_entries=entries, min_annotated=o.min_genomes,
        ntc_paths=ntc, contaminants=contaminants, ntc_ratio=o.ntc_ratio)

    kreport = spec.get("kraken_report")
    kmap, total_reads = None, None
    if kreport and os.path.exists(kreport):
        kmap = reads.parse_kreport(kreport)
        total_reads = reads.kreport_total(kreport)
        if total_reads:
            annot["rpm"] = (1e6 * annot["new_est_reads"] / total_reads).round(1)
        stats = reads.parse_kreport_stats(kreport)
        extra = pd.DataFrame([dict(TaxID=t, **{k: v for k, v in d.items() if k != "kraken_clade_reads"})
                              for t, d in stats.items()]) if stats else pd.DataFrame()
        if not extra.empty and len(extra.columns) > 1:
            annot = annot.merge(extra, on="TaxID", how="left")
    if "rpm" not in annot.columns:
        annot["rpm"] = annot["reads_per_million"]

    final = annot

    # ---- 2. read-level scoring + alignment
    can_score = all(spec.get(k) for k in ("kraken_output", "kraken_report", "r1"))
    if "score" in steps and can_score:
        cand = annot[annot["new_est_reads"] >= o.score_min_reads].head(o.score_top)
        sel = cand["TaxID"].astype(int).tolist()
        _log(f"[{sample}] scoring {len(sel)} species")
        ids = reads.collect_read_ids(spec["kraken_output"], sel, kmap)
        scored = score.score_sample(
            cand, spec, o.outdir, o.refs_dir or os.path.join(o.outdir, "refs"), kmap, ids,
            groups_for_deepac=tuple(o.deepac_groups.split(",")), deepac_mode=o.deepac_mode,
            deepac_cap=o.deepac_max_reads, threads=o.threads, api_key=o.api_key,
            do_deepac=not o.no_deepac, do_coverage=not o.no_coverage, keep_work=o.keep_work,
            total_reads=total_reads or float(annot["new_est_reads"].sum()), log=_log)
        new_cols = [c for c in scored.columns if c not in annot.columns]
        final = annot.merge(scored[["TaxID"] + new_cols], on="TaxID", how="left")
    elif "score" in steps:
        _log(f"[{sample}] skipping read scoring (needs kraken_output, kraken_report and r1)")

    # ---- 3. assembly + resistance genes
    amr_df = None
    if "amr" in steps:
        contigs = spec.get("contigs")
        work = os.path.join(o.outdir, "amr", sample)
        try:
            if not contigs and o.assemble and spec.get("r1"):
                _log(f"[{sample}] assembling with MEGAHIT")
                contigs = amr.assemble(spec["r1"], spec.get("r2"), os.path.join(work, "megahit"), o.threads)
            if contigs:
                _log(f"[{sample}] running AMRFinderPlus on {contigs}")
                tsv = amr.run_amrfinder(contigs, work, o.threads)
                amr_df = amr.parse_amrfinder(tsv, o.amr_min_identity, o.amr_min_coverage)
                ctaxa = {}
                if spec.get("contig_kraken"):
                    ctaxa = amr.read_contig_taxa(spec["contig_kraken"])
                elif o.kraken_db:
                    ctaxa = amr.read_contig_taxa(amr.classify_contigs(contigs, o.kraken_db, work, o.threads))
                inv = {t: sp for sp, ts in (kmap or {}).items() for t in ts}
                names = dict(zip(annot["TaxID"], annot["species"]))
                amr_df = amr.attribute(amr_df, ctaxa, inv, names)
                summary = amr.per_species(amr_df)
                if summary:
                    sm = pd.DataFrame([dict(TaxID=t, **d) for t, d in summary.items()])
                    final = final.merge(sm, on="TaxID", how="left")
                amr_df.to_csv(os.path.join(o.outdir, f"{sample}_amr.csv"), index=False)
            else:
                _log(f"[{sample}] skipping resistance genes (give contigs, or --assemble with reads)")
        except Exception as exc:
            _log(f"[{sample}] resistance step failed: {exc}")

    # ---- 4. reporting tiers + report
    final = report.assign_report_tier(final, o.report_min_rpm, o.report_min_reads)
    final.to_csv(os.path.join(o.outdir, f"{sample}_final.csv"), index=False)
    if "report" in steps:
        meta = {"Sample": sample, "Bracken file": os.path.basename(spec["bracken"]),
                "Total reads (Kraken2 report)": total_reads,
                "Host reads removed (%)": spec.get("host_pct"),
                "Negative control(s)": ", ".join(os.path.basename(p) for p in ntc) or "none supplied",
                "Species annotated": len(final)}
        params = {"Min reads to report": o.report_min_reads, "Min RPM to report": o.report_min_rpm,
                  "NTC RPM-ratio cut-off": o.ntc_ratio, "Evidence sources": o.sources,
                  "AMR min identity / coverage %": f"{o.amr_min_identity} / {o.amr_min_coverage}",
                  "DeePaC mode": o.deepac_mode}
        with open(os.path.join(o.outdir, f"{sample}_report.html"), "w", encoding="utf-8") as fh:
            fh.write(report.build_html(sample, final, amr_df, meta, params))
    return final
