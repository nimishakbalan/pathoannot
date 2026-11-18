"""Contamination control: negative-control (NTC) comparison and a kit-contaminant genus list."""
from importlib import resources

import pandas as pd

from .priority import norm

DEFAULT_RATIO = 10.0


def load_contaminants(extra_path=None):
    with resources.files("pathoannot").joinpath("data/contaminants.csv").open("r", encoding="utf-8") as fh:
        genera = set(norm(g) for g in pd.read_csv(fh)["genus"])
    if extra_path:
        genera |= set(norm(g) for g in pd.read_csv(extra_path)["genus"])
    return genera


def _read(path):
    df = pd.read_csv(path, sep="\t").rename(columns={"taxonomy_id": "TaxID"})
    return df[["TaxID", "new_est_reads"]]


def ntc_table(paths):
    """{TaxID: (reads, rpm)}, taking the highest value seen across the negative controls."""
    best = {}
    for p in paths:
        df = _read(p)
        total = float(df["new_est_reads"].sum()) or 1.0
        for tid, n in zip(df["TaxID"], df["new_est_reads"]):
            rpm = 1e6 * n / total
            if int(tid) not in best or rpm > best[int(tid)][1]:
                best[int(tid)] = (int(n), rpm)
    return best


def flag(ratio, in_ntc, ntc_given, known, cut):
    if ntc_given:
        if in_ntc:
            if ratio < cut:
                return "Likely contaminant (NTC)"
            if ratio < 3 * cut:
                return "Possible contaminant (NTC)"
            return "Passes NTC"
        return "Passes NTC; kit-contaminant genus" if known else "Passes NTC (absent in control)"
    return "Possible contaminant (kit-contaminant genus)" if known else "No NTC provided"


def apply(df, ntc_paths=(), contaminants=None, ratio_cut=DEFAULT_RATIO):
    """Add NTC read counts, RPM ratio and a contamination flag to an annotate table."""
    contaminants = contaminants if contaminants is not None else load_contaminants()
    ntc = ntc_table(ntc_paths) if ntc_paths else {}
    given = bool(ntc_paths)
    cols = dict(ntc_reads=[], ntc_rpm=[], rpm_ratio_vs_ntc=[], known_contaminant_genus=[], contamination_flag=[])
    for _, r in df.iterrows():
        tid = int(r["TaxID"]) if pd.notna(r["TaxID"]) else -1
        genus = norm(str(r["species"]).split(" ")[0].replace("[", "").replace("]", ""))
        known = genus in contaminants
        in_ntc = tid in ntc
        reads, rpm = ntc.get(tid, (0, 0.0))
        ratio = (r["reads_per_million"] / rpm) if in_ntc and rpm > 0 else None
        cols["ntc_reads"].append(reads if given else None)
        cols["ntc_rpm"].append(round(rpm, 2) if given else None)
        cols["rpm_ratio_vs_ntc"].append(round(ratio, 2) if ratio is not None else None)
        cols["known_contaminant_genus"].append(known)
        cols["contamination_flag"].append(flag(ratio, in_ntc, given, known, ratio_cut))
    return df.assign(**cols)
