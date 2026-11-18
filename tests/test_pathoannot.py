import os

import pandas as pd
import pytest

from pathoannot import amr, contam, core, cli, evidence, literature, priority, reads, report, score


# ------------------------------------------------------------------ fakes
class FakeResp:
    def __init__(self, payload, status=200):
        self._p, self.status_code = payload, status

    def json(self):
        return self._p

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(self.status_code)


PARENTS = {573: [1, 2, 91347], 10298: [1, 10239], 5476: [1, 2759, 4751], 999999: [1, 2],
           562: [1, 2, 91347], 12345: [1, 2]}


def fake_get(url, params=None, headers=None, timeout=None):
    if "europepmc" in url:
        q = params["query"]
        if "Klebsiella pneumoniae" in q:
            n = 1000 if "infection" in q else (30 if "commensal" in q else (20 if "soil" in q else 2000))
        else:
            n = 100 if "infection" not in q and "commensal" not in q and "soil" not in q else 5
            if "soil" in q:
                n = 70
        return FakeResp({"hitCount": n})
    if "/taxonomy/taxon/" in url:
        tid = int(url.split("/taxonomy/taxon/")[1].split("/")[0])
        return FakeResp({"reports": [{"taxonomy": {"tax_id": tid, "parents": PARENTS.get(tid, [1])}}]})
    if "/virus/taxon/" in url:
        return FakeResp({"total_count": 800, "reports": [
            {"accession": "NC_001806.2", "host": {"sci_name": "Homo sapiens"}, "isolate": {"source": "throat swab"}}] * 20})
    if "/genome/taxon/" in url:
        tid = int(url.split("/genome/taxon/")[1].split("/")[0])
        if tid in (573, 562):
            reps = [{"assembly_info": {"biosample": {"host": "Homo sapiens", "isolation_source": "blood"}}}] * 20
        else:
            reps = [{"assembly_info": {"biosample": {"isolation_source": "soil"}}}] * 20
        return FakeResp({"total_count": len(reps), "reports": reps})
    return FakeResp({}, 404)


@pytest.fixture(autouse=True)
def patch_net(monkeypatch):
    monkeypatch.setattr(core.requests, "get", fake_get)
    monkeypatch.setattr(core.time, "sleep", lambda s: None)
    monkeypatch.setattr(literature.time, "sleep", lambda s: None)


def write_bracken(path):
    pd.DataFrame({
        "name": ["Klebsiella pneumoniae", "Soil bug", "Human alphaherpesvirus 1", "Escherichia coli"],
        "taxonomy_id": [573, 999999, 10298, 562],
        "taxonomy_lvl": ["S"] * 4,
        "kraken_assigned_reads": [500, 50, 40, 30],
        "added_reads": [10, 5, 3, 2],
        "new_est_reads": [510, 55, 43, 32],
        "fraction_total_reads": [0.8, 0.08, 0.07, 0.05],
    }).to_csv(path, sep="\t", index=False)


# ------------------------------------------------------------------ priority lists
def test_builtin_priority_matching():
    ent = priority.load_builtin()
    assert len(ent) == 31
    hits = priority.match("Klebsiella pneumoniae", 573, [1, 2, 91347], ent)
    assert any(e["list"] == "WHO BPPL 2024" and e["priority"] == "critical" for e in hits)      # via order clade
    assert priority.match("Salmonella enterica", 28901, [], ent)[0]["taxon" if False else "priority"] == "high"
    assert priority.match("[Candida] auris", 498019, [], ent)
    assert priority.match("Nakaseomyces glabratus", 5478, [], ent)
    assert priority.match("Fusarium oxysporum", 5507, [], ent)                                   # genus entry
    assert not priority.match("Streptococcus mitis", 28037, [], ent)


def test_user_priority_list(tmp_path):
    f = tmp_path / "icmr.csv"
    f.write_text("list,priority,taxon,aliases\nICMR,priority,Kingella kingae,K. kingae\n")
    ent = priority.load_user(str(f))
    assert priority.match("Kingella kingae", 504, [], ent)


# ------------------------------------------------------------------ literature
def test_literature_labels():
    s = literature.lit_summary("Klebsiella pneumoniae")
    assert literature.lit_label(s) == "Disease-dominated literature"
    assert literature.lit_label({"lit_total": 5}) == "Insufficient literature"
    env = dict(lit_total=100, lit_disease_frac=0.05, lit_commensal_frac=0.02, lit_environment_frac=0.7)
    assert literature.lit_label(env) == "Environmental-dominated literature"


# ------------------------------------------------------------------ annotate
def test_annotate_end_to_end(tmp_path):
    br = tmp_path / "s1.bracken"
    write_bracken(br)
    res = evidence.annotate_sample(str(br), cache={}, cache_path=str(tmp_path / "c.json"), log=open(os.devnull, "w"))
    call = dict(zip(res["species"], res["final_call"]))
    assert call["Klebsiella pneumoniae"] == "Likely pathogen/opportunist"   # priority + NCBI + literature
    assert call["Human alphaherpesvirus 1"] in ("Possible pathogen/opportunist", "Likely pathogen/opportunist")
    assert call["Soil bug"] == "Likely environmental"
    row = res[res["species"] == "Klebsiella pneumoniae"].iloc[0]
    assert "WHO BPPL 2024" in row["priority_hits"] and row["group"] == "Bacteria"
    assert row["reads_per_million"] > 0


def test_sources_can_be_disabled(tmp_path):
    br = tmp_path / "s1.bracken"
    write_bracken(br)
    res = evidence.annotate_sample(str(br), sources=("priority",), cache={}, cache_path=None, log=open(os.devnull, "w"))
    assert "ncbi_label" not in res.columns and "lit_label" not in res.columns


# ------------------------------------------------------------------ kraken helpers
KREPORT = (
    " 10.00\t10\t10\tU\t0\tunclassified\n"
    " 90.00\t90\t0\tR\t1\troot\n"
    " 80.00\t80\t0\tD\t2\t  Bacteria\n"
    " 60.00\t60\t5\tS\t573\t    Klebsiella pneumoniae\n"
    " 40.00\t40\t40\tS1\t1125630\t      Klebsiella pneumoniae subsp. x\n"
    " 20.00\t20\t20\tS\t562\t    Escherichia coli\n"
)


def test_parse_kreport(tmp_path):
    f = tmp_path / "s.kreport"
    f.write_text(KREPORT)
    m = reads.parse_kreport(str(f))
    assert m[573] == {573, 1125630} and m[562] == {562}


def test_collect_and_extract(tmp_path):
    kr = tmp_path / "s.kraken"
    kr.write_text("C\tr1\t573\t100\t573:5\nC\tr2\t1125630\t100\t1:1\nC\tr3\t562\t100\t562:1\nU\tr4\t0\t100\t0:1\n")
    ids = reads.collect_read_ids(str(kr), [573, 562], {573: {573, 1125630}, 562: {562}})
    assert ids[573] == {"r1", "r2"} and ids[562] == {"r3"}
    fq = tmp_path / "R1.fastq"
    fq.write_text("".join(f"@{r}/1\nACGT\n+\nIIII\n" for r in ["r1", "r2", "r3", "r4"]))
    out = reads.extract_reads([str(fq)], ids, str(tmp_path / "w"))
    assert sum(1 for _ in open(out[573][0])) == 8
    n = reads.fastq_to_fasta(out[573], str(tmp_path / "x.fa"))
    assert n == 2


# ------------------------------------------------------------------ coverage / tiers
COV = ("#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq\n"
       "chr\t1\t1000000\t1000\t500000\t50\t0.1500\t30\t60\n"
       "plasmid\t1\t100000\t100\t10000\t10\t0.0150\t30\t60\n")


def test_parse_coverage_and_tiers():
    cov = score.parse_coverage(COV)
    assert cov["ref_length_bp"] == 1100000 and cov["mapped_reads"] == 1100
    assert round(cov["breadth_pct"], 2) == round(100 * 510000 / 1100000, 2)
    assert score.presence_tier(cov) in ("Confirmed", "Probable", "Doubtful")
    even = dict(mapped_reads=500, breadth_ratio=0.95)
    piled = dict(mapped_reads=500, breadth_ratio=0.05)
    assert score.presence_tier(even) == "Confirmed"
    assert score.presence_tier(piled) == "Doubtful"
    assert score.presence_tier(None) == "Not assessed"


def test_summarise_scores():
    s = score.summarise_scores([0.9, 0.8, 0.2, 0.1])
    assert s["deepac_n_reads"] == 4 and s["deepac_pct_above_threshold"] == 50.0


# ------------------------------------------------------------------ full run (tools mocked)
def test_run_end_to_end(tmp_path, monkeypatch):
    br = tmp_path / "s1.bracken"
    write_bracken(br)
    kr = tmp_path / "s1.kraken"
    kr.write_text("".join(f"C\tk{i}\t573\t100\t573:5\n" for i in range(150)) + "".join(f"C\te{i}\t562\t100\t562:5\n" for i in range(120)))
    krep = tmp_path / "s1.kreport"
    krep.write_text(KREPORT)
    fq = tmp_path / "R1.fastq"
    fq.write_text("".join(f"@{p}{i}\nACGT\n+\nIIII\n" for p in ("k", "e") for i in range(150 if p == "k" else 120)))

    monkeypatch.setattr(score, "run_deepac", lambda fa, mode="rapid", binary="deepac": [0.9, 0.7, 0.2, 0.6])
    monkeypatch.setattr(score, "get_reference", lambda taxid, group, refs_dir, api_key=None: "ref.fna")
    monkeypatch.setattr(score, "map_and_cover", lambda ref, fq, wd, threads=4: score.parse_coverage(COV))

    out = tmp_path / "out"
    rc = cli.main(["run", "-o", str(out), "--sample", "s1", "--bracken", str(br), "--kraken-output", str(kr),
                   "--kraken-report", str(krep), "--r1", str(fq), "--cache", str(tmp_path / "c.json"),
                   "--score-min-reads", "100", "--steps", "score"])
    assert rc == 0
    df = pd.read_csv(out / "s1_final.csv")
    k = df[df["TaxID"] == 573].iloc[0]
    assert k["reads_extracted"] == 150 and k["deepac_n_reads"] == 4
    assert k["deepac_pct_above_threshold"] == 75.0
    assert k["mean_depth"] > 0 and k["breadth_pct"] > 0 and k["presence_tier"] in ("Confirmed", "Probable", "Doubtful")
    assert (out / "ALL_samples_final.csv").exists()


def test_check_command(capsys):
    assert cli.main(["check", "--taxid", "10298", "--name", "Klebsiella pneumoniae"]) == 0
    out = capsys.readouterr().out
    assert "Priority-list entries: 31" in out and "group=Virus" in out


# ------------------------------------------------------------------ contamination
def test_ntc_ratio_and_contaminant_list(tmp_path):
    br = tmp_path / "s1.bracken"
    write_bracken(br)
    ntc = tmp_path / "ntc.bracken"
    pd.DataFrame({"name": ["Klebsiella pneumoniae", "Ralstonia pickettii"], "taxonomy_id": [573, 329],
                  "taxonomy_lvl": ["S", "S"], "kraken_assigned_reads": [5, 80], "added_reads": [0, 0],
                  "new_est_reads": [5, 80], "fraction_total_reads": [0.06, 0.94]}).to_csv(ntc, sep="\t", index=False)
    res = evidence.annotate_sample(str(br), sources=("priority",), cache={}, cache_path=None,
                                   log=open(os.devnull, "w"), ntc_paths=(str(ntc),))
    kp = res[res["species"] == "Klebsiella pneumoniae"].iloc[0]
    assert kp["ntc_reads"] == 5 and kp["rpm_ratio_vs_ntc"] > 0
    assert kp["contamination_flag"] in ("Likely contaminant (NTC)", "Possible contaminant (NTC)", "Passes NTC")
    other = res[res["species"] == "Soil bug"].iloc[0]
    assert other["contamination_flag"] == "Passes NTC (absent in control)"
    no_ntc = evidence.annotate_sample(str(br), sources=("priority",), cache={}, cache_path=None,
                                      log=open(os.devnull, "w"))
    assert (no_ntc["contamination_flag"] == "No NTC provided").all()
    assert contam.flag(2.0, True, True, False, 10) == "Likely contaminant (NTC)"
    assert contam.flag(25.0, True, True, False, 10) == "Possible contaminant (NTC)"
    assert contam.flag(500.0, True, True, False, 10) == "Passes NTC"
    assert contam.flag(None, False, False, True, 10).startswith("Possible contaminant")


# ------------------------------------------------------------------ minimizers
def test_kreport_stats(tmp_path):
    f = tmp_path / "m.kreport"
    f.write_text(" 10.00\t10\t10\t0\t0\tU\t0\tunclassified\n 90.00\t90\t0\t0\t0\tR\t1\troot\n"
                 " 60.00\t60\t60\t1200\t400\tS\t573\t    Klebsiella pneumoniae\n")
    st = reads.parse_kreport_stats(str(f))
    assert st[573]["minimizers_distinct"] == 400 and st[573]["minimizer_dup_ratio"] == 3.0
    assert reads.kreport_total(str(f)) == 100


# ------------------------------------------------------------------ AMR
AMR_TSV = ("Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\tSubclass\tMethod\t"
           "% Coverage of reference sequence\t% Identity to reference sequence\n"
           "p1\tk141_1 flag=1\t1\t900\t+\tblaKPC-2\tKPC carbapenemase\tcore\tAMR\tAMR\tBETA-LACTAM\tCARBAPENEM\tEXACTX\t100\t99.8\n"
           "p2\tk141_2\t1\t900\t+\toqxA\tefflux\tcore\tAMR\tAMR\tQUINOLONE\tQUINOLONE\tPARTIALX\t60\t99\n"
           "p3\tk141_3\t1\t900\t+\tvirX\tvirulence\tplus\tVIRULENCE\tVIRULENCE\tNA\tNA\tEXACTX\t100\t100\n")


def test_amr_parse_and_attribute(tmp_path):
    f = tmp_path / "amr.tsv"
    f.write_text(AMR_TSV)
    df = amr.parse_amrfinder(str(f))
    assert list(df["gene"]) == ["blaKPC-2"]                     # low-coverage and non-AMR hits dropped
    df = amr.attribute(df, {"k141_1": 1125630}, {1125630: 573}, {573: "Klebsiella pneumoniae"})
    assert df["species"].iloc[0] == "Klebsiella pneumoniae"
    assert amr.per_species(df)[573]["amr_genes"] == "blaKPC-2"


# ------------------------------------------------------------------ report tiers + html
def test_report_tiers_and_html(tmp_path):
    df = pd.DataFrame([
        dict(species="Good bug", pathogen_class="Pathogen/opportunist", new_est_reads=500, rpm=900.0,
             presence_tier="Confirmed", contamination_flag="Passes NTC"),
        dict(species="Unconfirmed bug", pathogen_class="Pathogen/opportunist", new_est_reads=500, rpm=900.0,
             presence_tier="Doubtful", contamination_flag="Passes NTC"),
        dict(species="Reagent bug", pathogen_class="Pathogen/opportunist", new_est_reads=500, rpm=900.0,
             presence_tier="Confirmed", contamination_flag="Likely contaminant (NTC)", rpm_ratio_vs_ntc=1.2),
        dict(species="Gut bug", pathogen_class="Commensal", new_est_reads=500, rpm=900.0,
             presence_tier="Confirmed", contamination_flag="No NTC provided"),
    ])
    t = report.assign_report_tier(df)
    assert list(t["report_tier"]) == ["Reportable", "Review", "Excluded: likely contaminant", "Not reported"]
    page = report.build_html("S1", t, None, {"Sample": "S1"}, {"Min RPM": 10})
    assert "Good bug" in page and "Excluded as likely contaminants" in page and "sign-off" in page


def test_full_run_with_ntc_amr_and_report(tmp_path, monkeypatch):
    br = tmp_path / "s1.bracken"
    write_bracken(br)
    kr = tmp_path / "s1.kraken"
    kr.write_text("".join(f"C\tk{i}\t573\t100\t573:5\n" for i in range(150)))
    krep = tmp_path / "s1.kreport"
    krep.write_text(KREPORT)
    fq = tmp_path / "R1.fastq"
    fq.write_text("".join(f"@k{i}\nACGT\n+\nIIII\n" for i in range(150)))
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">k141_1 flag=1\nACGT\n")
    ck = tmp_path / "contigs.kraken"
    ck.write_text("C\tk141_1\t573\t4\t573:1\n")
    amr_tsv = tmp_path / "amr.tsv"
    amr_tsv.write_text(AMR_TSV)

    monkeypatch.setattr(score, "run_deepac", lambda fa, mode="rapid", binary="deepac": [0.9, 0.7, 0.2, 0.6])
    monkeypatch.setattr(score, "get_reference", lambda taxid, group, refs_dir, api_key=None: "ref.fna")
    monkeypatch.setattr(score, "map_and_cover", lambda ref, fq, wd, threads=4: score.parse_coverage(COV))
    monkeypatch.setattr(amr, "run_amrfinder", lambda contigs, outdir, threads=4, binary="amrfinder", plus=False: str(amr_tsv))

    out = tmp_path / "out"
    rc = cli.main(["run", "-o", str(out), "--sample", "s1", "--bracken", str(br), "--kraken-output", str(kr),
                   "--kraken-report", str(krep), "--r1", str(fq), "--contigs", str(contigs),
                   "--contig-kraken", str(ck), "--cache", str(tmp_path / "c.json"),
                   "--score-min-reads", "100", "--report-min-reads", "50", "--report-min-rpm", "1"])
    assert rc == 0
    df = pd.read_csv(out / "s1_final.csv")
    k = df[df["TaxID"] == 573].iloc[0]
    assert k["pathogen_class"] == "Pathogen/opportunist" and k["amr_genes"] == "blaKPC-2"
    assert k["report_tier"] in ("Reportable", "Review") and "rpm" in df.columns
    page = (out / "s1_report.html").read_text()
    assert "Klebsiella pneumoniae" in page and "blaKPC-2" in page
    assert (out / "s1_amr.csv").exists()
