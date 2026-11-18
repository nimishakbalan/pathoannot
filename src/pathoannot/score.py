"""Step 2 (score): DeePaC pathogenic-potential scores, depth / breadth of coverage, abundance."""
import glob
import io
import math
import os
import shutil
import subprocess
import zipfile

import requests

from . import core
from .reads import fastq_to_fasta

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
CONFIRMED_RATIO, PROBABLE_RATIO = 0.7, 0.3
CONFIRMED_READS, PROBABLE_READS, MIN_READS = 50, 20, 10


# ------------------------------------------------------------------ tools
def check_tools():
    return {t: shutil.which(t) for t in ("minimap2", "samtools", "deepac")}


def _run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


# ------------------------------------------------------------------ DeePaC
def run_deepac(fasta, mode="rapid", binary="deepac"):
    """Run `deepac predict` and return a flat list of per-read pathogenic potentials."""
    import numpy as np
    flag = "-s" if mode == "sensitive" else "-r"
    _run([binary, "predict", flag, fasta])
    pred = os.path.splitext(fasta)[0] + "_predictions.npy"
    arr = np.load(pred)
    return list(arr.reshape(len(arr), -1)[:, -1])


def summarise_scores(scores, threshold=0.5):
    if not scores:
        return dict(deepac_n_reads=0, deepac_mean_score=None, deepac_pct_above_threshold=None)
    return dict(deepac_n_reads=len(scores),
                deepac_mean_score=round(sum(scores) / len(scores), 4),
                deepac_pct_above_threshold=round(100 * sum(s >= threshold for s in scores) / len(scores), 1))


# ------------------------------------------------------------------ references
def _find_local_ref(taxid, refs_dir):
    for ext in ("fa", "fasta", "fna", "fa.gz", "fasta.gz", "fna.gz"):
        p = os.path.join(refs_dir, f"{taxid}.{ext}")
        if os.path.exists(p):
            return p
    return None


def _download_genome_ref(taxid, dest, api_key=None):
    data = core._get_json(core.GENOME_API.format(taxid=taxid), {"page_size": 100}, api_key) or {}
    reports = data.get("reports") or []
    if not reports:
        return None
    pick = next((r for r in reports if (r.get("assembly_info") or {}).get("refseq_category", "")
                 in ("reference genome", "representative genome")), reports[0])
    acc = pick["accession"]
    r = requests.get(f"{core.BASE}/genome/accession/{acc}/download",
                     params={"include_annotation_type": "GENOME_FASTA", "filename": f"{acc}.zip"},
                     headers={"api-key": api_key} if api_key else None, timeout=300)
    r.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(r.content)) as z:
        names = [n for n in z.namelist() if n.endswith((".fna", ".fa", ".fasta"))]
        if not names:
            return None
        with z.open(names[0]) as src, open(dest, "wb") as out:
            shutil.copyfileobj(src, out)
    return dest


def _download_virus_ref(taxid, dest, api_key=None):
    data = core._get_json(core.VIRUS_API.format(taxid=taxid), {"page_size": 500}, api_key) or {}
    reports = data.get("reports") or []
    if not reports:
        return None
    accs = [r.get("accession") for r in reports if r.get("accession")]
    acc = next((a for a in accs if a.startswith("NC_")), accs[0] if accs else None)
    if not acc:
        return None
    r = requests.get(EUTILS, params={"db": "nuccore", "id": acc, "rettype": "fasta", "retmode": "text"}, timeout=120)
    r.raise_for_status()
    with open(dest, "w") as out:
        out.write(r.text)
    return dest


def get_reference(taxid, group, refs_dir, api_key=None):
    """Local <refs_dir>/<taxid>.fa(sta|na)[.gz] if present, else best-effort NCBI download."""
    os.makedirs(refs_dir, exist_ok=True)
    local = _find_local_ref(taxid, refs_dir)
    if local:
        return local
    dest = os.path.join(refs_dir, f"{taxid}.fna")
    try:
        fn = _download_virus_ref if group == "Virus" else _download_genome_ref
        return fn(taxid, dest, api_key)
    except Exception:
        return None


# ------------------------------------------------------------------ mapping / coverage
def parse_coverage(text):
    """Aggregate `samtools coverage` output over all reference sequences."""
    length = cov = reads = 0
    depth_x_len = 0.0
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        start, end, nreads, covbases, meandepth = int(f[1]), int(f[2]), int(f[3]), int(f[4]), float(f[6])
        seq_len = end - start + 1
        length += seq_len
        cov += covbases
        reads += nreads
        depth_x_len += meandepth * seq_len
    if length == 0:
        return None
    depth = depth_x_len / length
    breadth = 100 * cov / length
    expected = 100 * (1 - math.exp(-depth))
    return dict(ref_length_bp=length, mapped_reads=reads, mean_depth=round(depth, 4),
                breadth_pct=round(breadth, 3), expected_breadth_pct=round(expected, 3),
                breadth_ratio=round(breadth / expected, 3) if expected > 0 else None)


def map_and_cover(ref, fastqs, workdir, threads=4):
    """minimap2 (short reads) -> sorted BAM -> samtools coverage; returns parse_coverage() dict."""
    os.makedirs(workdir, exist_ok=True)
    bam = os.path.join(workdir, "aln.bam")
    mm = subprocess.Popen(["minimap2", "-ax", "sr", "-t", str(threads), ref, *fastqs],
                          stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run(["samtools", "sort", "-@", str(threads), "-o", bam, "-"],
                   stdin=mm.stdout, check=True, stderr=subprocess.DEVNULL)
    mm.stdout.close()
    if mm.wait() != 0:
        raise RuntimeError("minimap2 failed")
    cov = _run(["samtools", "coverage", "-q", "1", bam]).stdout
    return parse_coverage(cov)


def presence_tier(cov):
    """Confirmed / Probable / Doubtful from how evenly reads cover the reference.

    Compares observed breadth with the breadth expected for the observed depth if reads were
    spread randomly (1 - exp(-depth)). Reads piled on a few regions give a low ratio."""
    if not cov:
        return "Not assessed"
    r, n = cov["breadth_ratio"], cov["mapped_reads"]
    if n < MIN_READS or r is None:
        return "Doubtful"
    if r >= CONFIRMED_RATIO and n >= CONFIRMED_READS:
        return "Confirmed"
    if r >= PROBABLE_RATIO and n >= PROBABLE_READS:
        return "Probable"
    return "Doubtful"


# ------------------------------------------------------------------ per-sample scoring
def score_sample(annot, spec, outdir, refs_dir, kreport_map, ids_by_species, groups_for_deepac=("Bacteria",),
                 deepac_mode="rapid", deepac_cap=50000, threads=4, api_key=None, do_deepac=True,
                 do_coverage=True, keep_work=False, total_reads=None, log=print):
    """Add DeePaC, depth/breadth and abundance columns to the annotate table for one sample."""
    from .reads import extract_reads
    sample = spec["sample"]
    work = os.path.join(outdir, "work", sample)
    fastqs = [spec["r1"]] + ([spec["r2"]] if spec.get("r2") else [])
    extracted = extract_reads(fastqs, ids_by_species, work)
    extra = []
    for _, row in annot.iterrows():
        tid, rec = row["TaxID"], {}
        if tid in extracted:
            files = extracted[tid]
            rec["reads_extracted"] = sum(1 for _ in open(files[0])) // 4
            if do_deepac and row["group"] in groups_for_deepac and rec["reads_extracted"] > 0:
                try:
                    fa = os.path.join(work, f"{tid}.fasta")
                    fastq_to_fasta(files, fa, deepac_cap)
                    rec.update(summarise_scores(run_deepac(fa, deepac_mode)))
                except Exception as exc:
                    log(f"  [{sample}] DeePaC failed for {tid}: {exc}")
            cov = None
            if do_coverage and rec["reads_extracted"] > 0:
                ref = get_reference(int(tid), row["group"], refs_dir, api_key)
                if ref:
                    try:
                        cov = map_and_cover(ref, files, os.path.join(work, str(tid)), threads)
                    except Exception as exc:
                        log(f"  [{sample}] mapping failed for {tid}: {exc}")
                else:
                    log(f"  [{sample}] no reference genome for {tid}; put <taxid>.fna in {refs_dir}")
            if cov:
                rec.update(cov)
                if total_reads:
                    rec["rpkm"] = round(cov["mapped_reads"] / (cov["ref_length_bp"] / 1e3) / (total_reads / 1e6), 3)
            rec["presence_tier"] = presence_tier(cov)
        extra.append(rec)
    import pandas as pd
    out = pd.concat([annot.reset_index(drop=True), pd.DataFrame(extra)], axis=1)
    if not keep_work:
        shutil.rmtree(work, ignore_errors=True)
    return out
