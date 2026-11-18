"""Kraken2 helpers: species -> descendant TaxIDs, per-read assignments, FASTQ extraction."""
import gzip
import os
import re


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def parse_kreport(path):
    """Return {species_taxid: set(taxids of the species and everything below it)}."""
    species, stack = {}, []  # stack of (depth, rank, taxid)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            rank, taxid, name = f[-3].strip(), f[-2].strip(), f[-1]
            try:
                taxid = int(taxid)
            except ValueError:
                continue
            depth = (len(name) - len(name.lstrip(" "))) // 2
            while stack and stack[-1][0] >= depth:
                stack.pop()
            if rank == "S":
                species[taxid] = {taxid}
            elif rank.startswith("S"):
                for _, prank, ptid in reversed(stack):
                    if prank == "S":
                        species[ptid].add(taxid)
                        break
            stack.append((depth, rank, taxid))
    return species


def collect_read_ids(kraken_output, wanted_species, kreport_map, max_per_species=None):
    """Stream a Kraken2 per-read output; return {species_taxid: set(read_ids)} for wanted species."""
    tax_to_species = {}
    for sp in wanted_species:
        for t in kreport_map.get(sp, {sp}):
            tax_to_species[t] = sp
    ids = {sp: set() for sp in wanted_species}
    pat = re.compile(r"taxid (\d+)")
    with open(kraken_output) as fh:
        for line in fh:
            p = line.split("\t", 4)
            if len(p) < 3 or p[0] != "C":
                continue
            m = pat.search(p[2])
            try:
                tax = int(m.group(1)) if m else int(p[2])
            except ValueError:
                continue
            sp = tax_to_species.get(tax)
            if sp is not None and (max_per_species is None or len(ids[sp]) < max_per_species):
                ids[sp].add(p[1])
    return ids


def extract_reads(fastq_paths, ids_by_species, outdir):
    """Write per-species FASTQ files for each mate. Returns {taxid: [paths]}."""
    os.makedirs(outdir, exist_ok=True)
    id_to_sp = {rid: sp for sp, s in ids_by_species.items() for rid in s}
    outputs = {sp: [] for sp in ids_by_species}
    for mate, fq in enumerate(fastq_paths, start=1):
        handles = {sp: open(os.path.join(outdir, f"{sp}_R{mate}.fastq"), "w") for sp in ids_by_species}
        try:
            with _open(fq) as fh:
                while True:
                    h = fh.readline()
                    if not h:
                        break
                    seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
                    rid = h[1:].split()[0]
                    if rid.endswith("/1") or rid.endswith("/2"):
                        rid = rid[:-2]
                    sp = id_to_sp.get(rid)
                    if sp is not None:
                        handles[sp].write(h + seq + plus + qual)
        finally:
            for sp, hd in handles.items():
                hd.close()
                outputs[sp].append(os.path.join(outdir, f"{sp}_R{mate}.fastq"))
    return outputs


def fastq_to_fasta(fastq_paths, fasta_out, max_reads=None):
    n = 0
    with open(fasta_out, "w") as out:
        for fq in fastq_paths:
            with open(fq) as fh:
                while True:
                    h = fh.readline()
                    if not h:
                        break
                    seq = fh.readline()
                    fh.readline()
                    fh.readline()
                    out.write(f">{h[1:].split()[0]}_{n}\n{seq}")
                    n += 1
                    if max_reads and n >= max_reads:
                        return n
    return n


def parse_kreport_stats(path):
    """Per-species clade reads and (if the report has them) minimizer counts, a genome-breadth proxy.
    Minimizer columns exist when Kraken2 was run with --report-minimizer-data."""
    stats = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6 or f[-3].strip() != "S":
                continue
            try:
                tid = int(f[-2])
                d = {"kraken_clade_reads": int(float(f[1]))}
                if len(f) >= 8:
                    d["minimizers_total"], d["minimizers_distinct"] = int(f[3]), int(f[4])
                    d["minimizer_dup_ratio"] = round(d["minimizers_total"] / d["minimizers_distinct"], 2) if d["minimizers_distinct"] else None
            except ValueError:
                continue
            stats[tid] = d
    return stats


def kreport_total(path):
    """Total reads in the sample (unclassified + root clade) from a Kraken2 report."""
    total = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6 and f[-3].strip() in ("U", "R") and f[-2].strip() in ("0", "1"):
                total += int(float(f[1]))
    return total
