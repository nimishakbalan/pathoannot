"""pathoannot command line: annotate | run | report | check."""
import argparse
import os
import sys

import pandas as pd

from . import __version__, amr, contam, core, evidence, literature, pipeline, priority, report, score


def _evidence_args(ap):
    ap.add_argument("-o", "--outdir", default="pathoannot_out", help="output directory (default pathoannot_out)")
    ap.add_argument("-m", "--min-reads", type=int, default=10, help="ignore species below this Bracken read count (default 10)")
    ap.add_argument("--sources", default="priority,literature,ncbi",
                    help="class evidence: any of priority,literature,ncbi (default all)")
    ap.add_argument("--priority-list", action="append", default=[], metavar="CSV",
                    help="extra priority list, e.g. ICMR (columns: taxon[,aliases,taxid,list,priority,note]); repeatable")
    ap.add_argument("--ntc-bracken", nargs="*", default=[], metavar="FILE",
                    help="Bracken file(s) of negative controls (NTC) for the contamination ratio")
    ap.add_argument("--ntc-ratio", type=float, default=contam.DEFAULT_RATIO,
                    help="sample/NTC RPM ratio below which a species is a likely contaminant (default 10)")
    ap.add_argument("--contaminant-list", help="CSV with a 'genus' column to add to the built-in kit-contaminant genera")
    ap.add_argument("--min-genomes", type=int, default=core.MIN_ANNOTATED, help="min NCBI records with host info (default 5)")
    ap.add_argument("--cache", default=core.DEFAULT_CACHE, help=f"lookup cache (default {core.DEFAULT_CACHE})")
    ap.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"), help="NCBI API key (or env NCBI_API_KEY)")


def _report_args(ap):
    ap.add_argument("--report-min-reads", type=float, default=report.DEFAULT_MIN_READS, help="min Bracken reads to report (default 50)")
    ap.add_argument("--report-min-rpm", type=float, default=report.DEFAULT_MIN_RPM, help="min reads per million to report (default 10)")


def _run_args(ap):
    ap.add_argument("--sample-sheet", help="TSV/CSV; columns: sample,bracken[,kraken_report,kraken_output,r1,r2,contigs,contig_kraken,host_pct]")
    ap.add_argument("--sample"), ap.add_argument("--bracken")
    ap.add_argument("--kraken-report"), ap.add_argument("--kraken-output")
    ap.add_argument("--r1", help="host-depleted FASTQ read 1 (.gz ok)"), ap.add_argument("--r2")
    ap.add_argument("--contigs", help="existing assembly (e.g. MEGAHIT final.contigs.fa)")
    ap.add_argument("--contig-kraken", help="existing Kraken2 output on the contigs (to attribute resistance genes to species)")
    ap.add_argument("--host-pct", type=float, help="percent of reads removed as host (shown in the report)")
    ap.add_argument("--steps", default="score,amr,report", help="steps after annotation: score,amr,report (default all; skipped if inputs missing)")
    ap.add_argument("--refs-dir", help="reference genomes named <taxid>.fna (default <outdir>/refs; missing ones are downloaded)")
    ap.add_argument("--score-min-reads", type=int, default=100, help="score species with >= this many reads (default 100)")
    ap.add_argument("--score-top", type=int, default=50, help="score at most this many species (default 50)")
    ap.add_argument("--deepac-groups", default="Bacteria"), ap.add_argument("--deepac-mode", choices=["rapid", "sensitive"], default="rapid")
    ap.add_argument("--deepac-max-reads", type=int, default=50000)
    ap.add_argument("--no-deepac", action="store_true"), ap.add_argument("--no-coverage", action="store_true")
    ap.add_argument("--assemble", action="store_true", help="assemble r1/r2 with MEGAHIT when no --contigs is given")
    ap.add_argument("--kraken-db", help="Kraken2 database, to classify contigs when --contig-kraken is not given")
    ap.add_argument("--amr-min-identity", type=float, default=90.0), ap.add_argument("--amr-min-coverage", type=float, default=80.0)
    ap.add_argument("--threads", type=int, default=4), ap.add_argument("--keep-work", action="store_true")


def build_parser():
    ap = argparse.ArgumentParser(prog="pathoannot", description="From Kraken2/Bracken results to a pathogen report.")
    ap.add_argument("--version", action="version", version=f"pathoannot {__version__}")
    sub = ap.add_subparsers(dest="cmd", required=True)

    a = sub.add_parser("annotate", help="class + contamination evidence from Bracken only (fast)")
    a.add_argument("-b", "--bracken", nargs="+", required=True)
    _evidence_args(a)

    r = sub.add_parser("run", help="full pipeline: evidence, contamination, alignment, DeePaC, AMR, report")
    _evidence_args(r), _run_args(r), _report_args(r)

    p = sub.add_parser("report", help="rebuild the HTML report from existing CSV files")
    p.add_argument("--final", required=True, help="<sample>_final.csv or _annotated.csv")
    p.add_argument("--amr", help="<sample>_amr.csv")
    p.add_argument("--sample", help="sample name for the title")
    p.add_argument("-o", "--outdir", default="pathoannot_out")
    _report_args(p)

    c = sub.add_parser("check", help="test tools, NCBI, Europe PMC and lists")
    c.add_argument("--taxid", type=int, nargs="*", default=[])
    c.add_argument("--name", nargs="*", default=[])
    c.add_argument("--priority-list", action="append", default=[])
    c.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    return ap


def _entries(args):
    ent = priority.load_builtin()
    for p in getattr(args, "priority_list", []):
        ent += priority.load_user(p)
    return ent


def _specs(args):
    if args.sample_sheet:
        sheet = pd.read_csv(args.sample_sheet, sep=None, engine="python")
        if not {"sample", "bracken"} <= set(sheet.columns):
            sys.exit("sample sheet needs at least the columns: sample, bracken")
        return [{k: (None if pd.isna(v) else v) for k, v in row.items()} for row in sheet.to_dict("records")]
    if not (args.sample and args.bracken):
        sys.exit("give --sample-sheet, or at least --sample and --bracken")
    return [dict(sample=args.sample, bracken=args.bracken, kraken_report=args.kraken_report,
                 kraken_output=args.kraken_output, r1=args.r1, r2=args.r2, contigs=args.contigs,
                 contig_kraken=args.contig_kraken, host_pct=args.host_pct)]


def _check(args):
    print("External tools:")
    tools = score.check_tools()
    import shutil
    for t in ("megahit", "amrfinder", "kraken2"):
        tools[t] = shutil.which(t)
    for t, p in tools.items():
        print(f"  {t}: {p or 'NOT FOUND'}")
    print(f"Priority-list entries: {len(_entries(args))}; kit-contaminant genera: {len(contam.load_contaminants())}")
    for tid in args.taxid:
        try:
            group, _ = core.get_lineage(tid, args.api_key)
            s = core.lookup_ncbi(tid, group, args.api_key)
            print(f"TaxID {tid}: group={group}; source={s['record_source']}; records={s['n_records_total']}; "
                  f"label={core.classify(s)}; human%={s['pct_human_host']}; clinical%={s['pct_clinical_source']}")
        except Exception as exc:
            print(f"TaxID {tid}: FAILED ({exc})")
    for nm in args.name:
        try:
            ls = literature.lit_summary(nm)
            print(f"{nm}: {ls} -> {literature.lit_label(ls)}")
        except Exception as exc:
            print(f"{nm}: Europe PMC FAILED ({exc})")
    return 0


def _summary(df):
    cols = [c for c in ("sample", "species", "pathogen_class", "new_est_reads", "rpm", "presence_tier",
                        "contamination_flag", "report_tier") if c in df.columns]
    if "report_tier" in df.columns:
        hot = df[df["report_tier"].isin(["Reportable", "Review"])]
    else:
        hot = df[df["pathogen_class"] == "Pathogen/opportunist"]
    print("\nSpecies to look at (Reportable and Review):")
    print(hot[cols].head(40).to_string(index=False) if len(hot) else "  none")


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.cmd == "check":
        return _check(args)

    os.makedirs(args.outdir, exist_ok=True)
    if args.cmd == "report":
        df = pd.read_csv(args.final)
        amr_df = pd.read_csv(args.amr) if args.amr else None
        df = report.assign_report_tier(df, args.report_min_rpm, args.report_min_reads)
        name = args.sample or (df["sample"].iloc[0] if "sample" in df.columns and len(df) else "sample")
        path = os.path.join(args.outdir, f"{name}_report.html")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(report.build_html(name, df, amr_df, {"Sample": name}, {"Min reads": args.report_min_reads,
                                                                              "Min RPM": args.report_min_rpm}))
        print(f"Wrote {path}")
        return 0

    cache = core.load_cache(args.cache)
    entries = _entries(args)
    contaminants = contam.load_contaminants(args.contaminant_list)

    if args.cmd == "annotate":
        results = []
        for path in args.bracken:
            res = evidence.annotate_sample(path, min_reads=args.min_reads, sources=tuple(args.sources.split(",")),
                                           cache=cache, cache_path=args.cache, api_key=args.api_key,
                                           priority_entries=entries, min_annotated=args.min_genomes,
                                           ntc_paths=tuple(args.ntc_bracken), contaminants=contaminants,
                                           ntc_ratio=args.ntc_ratio)
            res.to_csv(os.path.join(args.outdir, f"{evidence.sample_name(path)}_annotated.csv"), index=False)
            results.append(res)
        allres = pd.concat(results, ignore_index=True)
        allres.to_csv(os.path.join(args.outdir, "ALL_samples_annotated.csv"), index=False)
        _summary(allres)
        return 0

    finals = [pipeline.process_sample(spec, args, cache, entries, contaminants) for spec in _specs(args)]
    allres = pd.concat(finals, ignore_index=True)
    allres.to_csv(os.path.join(args.outdir, "ALL_samples_final.csv"), index=False)
    _summary(allres)
    print(f"\nResults written to {args.outdir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
