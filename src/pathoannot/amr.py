"""Assembly and resistance genes: MEGAHIT, contig classification (Kraken2), AMRFinderPlus."""
import os
import shutil
import subprocess

import pandas as pd


def _run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def assemble(r1, r2, outdir, threads=4, binary="megahit"):
    """Assemble (host-depleted) reads with MEGAHIT; returns the final.contigs.fa path."""
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    cmd = [binary, "-t", str(threads), "-o", outdir] + (["-1", r1, "-2", r2] if r2 else ["-r", r1])
    _run(cmd)
    return os.path.join(outdir, "final.contigs.fa")


def classify_contigs(contigs, db, outdir, threads=4, binary="kraken2"):
    """Run Kraken2 on contigs so each resistance gene can be attributed to a species."""
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "contigs.kraken")
    _run([binary, "--db", db, "--threads", str(threads), "--output", out,
          "--report", os.path.join(outdir, "contigs.kreport"), contigs])
    return out


def read_contig_taxa(kraken_output):
    taxa = {}
    with open(kraken_output) as fh:
        for line in fh:
            p = line.split("\t", 4)
            if len(p) >= 3 and p[0] == "C":
                try:
                    taxa[p[1]] = int(p[2])
                except ValueError:
                    pass
    return taxa


def run_amrfinder(contigs, outdir, threads=4, binary="amrfinder", plus=False):
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "amrfinder.tsv")
    cmd = [binary, "-n", contigs, "--threads", str(threads), "-o", out] + (["--plus"] if plus else [])
    _run(cmd)
    return out


def parse_amrfinder(path, min_identity=90.0, min_coverage=80.0, types=("AMR",)):
    """Read an AMRFinderPlus table, keep confident AMR hits."""
    df = pd.read_csv(path, sep="\t")
    pick = {"Contig id": "contig", "Element symbol": "gene", "Element name": "gene_name", "Type": "type",
            "Subtype": "subtype", "Class": "class", "Subclass": "subclass",
            "% Identity to reference sequence": "identity_pct", "% Coverage of reference sequence": "coverage_pct",
            "Method": "method"}
    df = df.rename(columns=pick)
    for c in pick.values():
        if c not in df.columns:
            df[c] = None
    df = df[list(pick.values())].copy()
    if types:
        df = df[df["type"].isin(types)]
    df = df[(pd.to_numeric(df["identity_pct"], errors="coerce") >= min_identity) &
            (pd.to_numeric(df["coverage_pct"], errors="coerce") >= min_coverage)]
    df["contig"] = df["contig"].astype(str).str.split().str[0]
    return df.reset_index(drop=True)


def attribute(amr_df, contig_taxa, strain_to_species, species_names):
    """Add species TaxID / name to each resistance hit using the contig's Kraken2 assignment."""
    tids, names = [], []
    for c in amr_df["contig"]:
        t = contig_taxa.get(c)
        sp = strain_to_species.get(t, t) if t is not None else None
        tids.append(sp)
        names.append(species_names.get(sp, f"taxid {sp}") if sp is not None else "unassigned")
    return amr_df.assign(species_taxid=tids, species=names)


def per_species(amr_df):
    out = {}
    if amr_df is None or amr_df.empty or "species_taxid" not in amr_df.columns:
        return out
    for tid, g in amr_df.dropna(subset=["species_taxid"]).groupby("species_taxid"):
        out[int(tid)] = dict(n_amr_genes=g["gene"].nunique(),
                             amr_genes=", ".join(sorted(g["gene"].dropna().unique())[:12]),
                             amr_classes=", ".join(sorted(g["class"].dropna().unique())))
    return out
