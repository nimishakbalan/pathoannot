# pathoannot

**From Kraken2/Bracken results to a pathogen report.**
For every species in a metagenomic sample, `pathoannot` answers four questions and puts the answers in one table
and one HTML report:

| Question | How it is answered |
|---|---|
| Is it a pathogen, a commensal or an environmental organism? | WHO priority lists (built in) and your own lists such as ICMR, Europe PMC literature, NCBI host / isolation-source metadata |
| Is it a contaminant? | Reads per million (RPM) compared with a negative control (NTC), plus a list of common kit-contaminant genera |
| Is it really in the sample? | Its reads are aligned back to a reference genome; depth, breadth and observed-vs-expected breadth give a presence tier |
| Do its reads look pathogenic, and is it resistant? | DeePaC read-level scores (bacteria); MEGAHIT assembly and AMRFinderPlus resistance genes, attributed to species |

> **Evidence, not diagnosis.** Every call is meant to help an expert decide. Thresholds are defaults and have not
> been validated. Clinical use needs validation (for example CLIA or ISO 15189) and sign-off by a qualified person.

---

## The workflow

```
 raw reads --(fastp QC, host removal: your own steps)--> host-depleted FASTQ
        |
        +--> Kraken2 (--report --output [--report-minimizer-data]) --> Bracken
        |
        +--> (optional) MEGAHIT assembly --> contigs

                      pathoannot run
 ┌───────────────────────────────────────────────────────────────────────┐
 │ 1 Class        priority lists + literature + NCBI metadata            │
 │ 2 Contamination  RPM ratio vs NTC + kit-contaminant genera            │
 │ 3 Abundance    reads, RPM, minimizer counts (from the Kraken2 report) │
 │ 4 Confirmation extract species reads -> align to genome (minimap2)    │
 │                -> depth, breadth, expected breadth -> presence tier   │
 │ 5 ML score     DeePaC on the species' reads                           │
 │ 6 Resistance   AMRFinderPlus on contigs, attributed to species        │
 │ 7 Report       tiers + HTML report with methods and sign-off block    │
 └───────────────────────────────────────────────────────────────────────┘
```

Steps 1 and 2 need only the Bracken file. Steps 3 to 7 use the extra inputs when you give them, and are skipped
(with a message) when they are missing, so you can start simple and add inputs later.

### 1. Pathogen class (commensal / pathogen / environmental)

Three independent lines of evidence are collected for each species.

- **Priority lists.** The built-in WHO Bacterial Priority Pathogens List 2024 and WHO Fungal Priority Pathogens List 2022
  are matched by species name, alias (for example *Candida glabrata* = *Nakaseomyces glabratus*), genus entry
  (*Fusarium* spp.) and clade entry (Enterobacterales, Mucorales). Add ICMR or any other list with `--priority-list`.
- **Journals.** Europe PMC titles and abstracts are searched for the species alone and combined with disease terms,
  commensal / microbiome terms and environmental terms. The three shares give a literature label:
  `Disease-dominated` (disease share >= 50%), `Environmental-dominated` (environment >= 40% and disease < 25%),
  `Commensal/microbiome-dominated` (commensal >= 15% and disease < 35%), `Mixed`, or `Insufficient literature` (< 20 papers).
- **NCBI metadata.** Genome records (NCBI Virus for viruses) are summarised for the share from humans, from clinical
  sites and from environmental sources: `Strong human/clinical association` (>= 50% human or >= 40% clinical),
  `Moderate` (>= 15%), `Mostly environmental` (>= 50% environmental and < 5% human), `Weak/unclear`,
  `Insufficient metadata`, `No genome records`.

Each line that points to disease adds a point to `pathogen_support` (priority list 1, NCBI strong 1 / moderate 0.5,
disease-dominated literature 1). The `final_call` and `pathogen_class` follow:

| `pathogen_support` | `final_call` | `pathogen_class` |
|---|---|---|
| 2 or more | Likely pathogen/opportunist | Pathogen/opportunist |
| 1 to 1.5 | Possible pathogen/opportunist | Pathogen/opportunist (flagged `REVIEW`) |
| below 1, moderate human association or commensal-dominated literature | Likely commensal / human-associated | Commensal |
| below 1, environmental NCBI or literature evidence | Likely environmental | Environmental |
| below 1, nothing notable or no usable evidence | Low pathogen evidence / Unresolved | Unresolved |

### 2. Contamination

- **Negative-control ratio.** Each species' RPM in the sample is divided by its RPM in the NTC (the highest RPM if you give several).
  A ratio below `--ntc-ratio` (default 10) is a `Likely contaminant (NTC)`; below three times that is a `Possible contaminant (NTC)`.
  A species absent from the NTC is `Passes NTC`.
- **Kit-contaminant genera.** A built-in starter list of genera often reported as reagent contaminants (for example
  *Ralstonia*, *Bradyrhizobium*, *Cutibacterium*, *Malassezia*) is checked. Without an NTC this gives
  `Possible contaminant (kit-contaminant genus)`. Add genera with `--contaminant-list`. A match is a warning only: some of these
  genera also contain real pathogens.

Run a negative control with every batch. The ratio is only as good as the control.

### 3 and 4. Abundance, depth, breadth and alignment confirmation

- `new_est_reads` and `fraction_total_reads` come from Bracken. `rpm` is reads per million of **all** reads in the sample,
  taken from the Kraken2 report (unclassified plus classified). `reads_per_million` is per million Bracken-classified reads and is used for the NTC comparison.
- If Kraken2 was run with `--report-minimizer-data`, `minimizers_distinct` and `minimizer_dup_ratio` are added. Many distinct minimizers
  suggest reads spread across the genome; a high duplication ratio suggests reads piled on few regions. This works without a reference genome.
- The species' reads (including any strain below it) are pulled from your FASTQ files using the Kraken2 per-read output and
  aligned with **minimap2** to a reference genome. **samtools coverage** gives `mean_depth`, `breadth_pct`, `expected_breadth_pct`
  (`100 x (1 - e^-depth)`, the breadth expected if reads were spread randomly), `breadth_ratio` (observed / expected) and `rpkm`.

| `presence_tier` | Rule |
|---|---|
| Confirmed | breadth_ratio >= 0.7 and >= 50 mapped reads |
| Probable | breadth_ratio >= 0.3 and >= 20 mapped reads |
| Doubtful | otherwise, or fewer than 10 mapped reads (reads piled on a few regions, as with conserved genes or cross-mapping) |
| Not assessed | no reference genome or alignment not run |

A divergent strain can map poorly to a single reference, so Doubtful means "not confirmed", not "absent".

### 5. ML score for pathogenic reads

**DeePaC** (`deepac predict`) scores the species' reads. The table gets `deepac_n_reads`, `deepac_mean_score` and
`deepac_pct_above_threshold` (share of reads scoring at least 0.5). The built-in DeePaC models cover **bacteria**, so other groups are not scored.
Read-level scores are noisy: compare a species with commensals in the same run before reading much into a number.

### 6. Assembly and resistance genes

Give an existing assembly with `--contigs` (for example MEGAHIT's `final.contigs.fa`), or add `--assemble` to run MEGAHIT on `--r1/--r2`.
**AMRFinderPlus** runs on the contigs; hits below 90% identity or 80% coverage, and non-AMR hits, are dropped (`--amr-min-identity`, `--amr-min-coverage`).
To say which species carries a gene, the contigs are classified with Kraken2: pass `--contig-kraken` (an existing Kraken2 output on the contigs) or `--kraken-db`.
Without either, genes are listed as `unassigned`. Plasmid-borne genes can sit on contigs assigned to another species, so attribution is indicative.

### 7. Reporting tiers and the HTML report

| `report_tier` | Rule |
|---|---|
| **Reportable** | class Pathogen/opportunist, reads >= `--report-min-reads` (50), RPM >= `--report-min-rpm` (10), presence Confirmed or Probable, no contamination warning |
| **Review** | pathogen evidence but low abundance, alignment not confirmed or not run, or a possible-contaminant warning; also unresolved species that are Confirmed present |
| **Excluded: likely contaminant** | NTC ratio says likely contaminant |
| **Not reported** | commensal, environmental or low-evidence species |

`<sample>_report.html` lists the reportable and review organisms with all the evidence side by side, resistance genes, excluded contaminants, methods,
the thresholds used, limitations and a sign-off block. Open it in a browser and print to PDF.

---

## Install

```bash
git clone https://github.com/<your-user>/pathoannot.git
cd pathoannot
pip install .
```

Python 3.9+, pandas, numpy and requests are installed automatically. The read-level steps need external programs on your `PATH`:

```bash
conda create -n pathoannot_env -c conda-forge -c bioconda python=3.9 \
      minimap2 samtools kraken2 bracken megahit ncbi-amrfinderplus fastp bowtie2
conda activate pathoannot_env
amrfinder -u                                            # download the AMRFinderPlus database once
pip install .
pip install "protobuf==3.20.3" "h5py<3.0" deepac       # DeePaC; see its own documentation if this fails
```

Check the setup before a real run:

```bash
pathoannot check --taxid 573 10298 5476 --name "Klebsiella pneumoniae"
```

It lists which tools were found and tests the NCBI and Europe PMC lookups for a bacterium (573), a virus (10298) and a fungus (5476).
`pathoannot annotate` with `--sources priority` works fully offline.

---

## Prepare the inputs

```bash
# 1. QC and host removal (your own settings)
fastp -i raw_R1.fq.gz -I raw_R2.fq.gz -o qc_R1.fq.gz -O qc_R2.fq.gz --dedup
bowtie2 -x human_index -1 qc_R1.fq.gz -2 qc_R2.fq.gz --very-sensitive -p 16 \
        --un-conc-gz sample_nohost_R%.fq.gz -S /dev/null

# 2. Classify the host-depleted reads; keep BOTH the report and the per-read output
DB=/path/to/k2_pluspf        # PlusPF includes fungi and protozoa
kraken2 --db $DB --threads 16 --confidence 0.1 --paired --report-minimizer-data \
        --report sample.kreport --output sample.kraken sample_nohost_R1.fq.gz sample_nohost_R2.fq.gz
bracken -d $DB -i sample.kreport -o sample.bracken -r 150 -l S -t 10

# 3. Assemble (or let pathoannot do it with --assemble)
megahit -1 sample_nohost_R1.fq.gz -2 sample_nohost_R2.fq.gz -o assembly/sample -t 16
```

Run the **negative control** (extraction blank / no-template control) through exactly the same steps and keep its `.bracken` file.
`--report-minimizer-data` is optional; it only adds the minimizer columns.

---

## Commands

### `pathoannot annotate`: class and contamination from Bracken only

```bash
pathoannot annotate -b sample1.bracken sample2.bracken -o results --ntc-bracken ntc.bracken
```

### `pathoannot run`: the full pipeline

One sample:

```bash
pathoannot run -o results --sample PT001 \
  --bracken PT001.bracken --kraken-report PT001.kreport --kraken-output PT001.kraken \
  --r1 PT001_nohost_R1.fq.gz --r2 PT001_nohost_R2.fq.gz \
  --contigs assembly/PT001/final.contigs.fa --kraken-db $DB \
  --ntc-bracken NTC01.bracken --priority-list icmr_priority.csv \
  --refs-dir refs --host-pct 92.5 --threads 16
```

Many samples: write a sample sheet (tab or comma separated; only `sample` and `bracken` are required; see `examples/samples.tsv`):

```
sample  bracken  kraken_report  kraken_output  r1  r2  contigs  contig_kraken  host_pct
```

```bash
pathoannot run -o results --sample-sheet samples.tsv --ntc-bracken NTC01.bracken --kraken-db $DB --refs-dir refs --threads 16
```

| Option | Meaning |
|---|---|
| `-o, --outdir` | output folder (default `pathoannot_out`) |
| `-m, --min-reads` | ignore species below this Bracken read count (default 10) |
| `--sources` | class evidence: any of `priority,literature,ncbi` (default all) |
| `--priority-list CSV` | add a list such as ICMR (repeatable) |
| `--ntc-bracken FILE...` | negative-control Bracken file(s) |
| `--ntc-ratio` | sample/NTC RPM ratio below which a species is a likely contaminant (default 10) |
| `--contaminant-list CSV` | extra kit-contaminant genera (column `genus`) |
| `--steps` | steps after annotation: any of `score,amr,report` (default all; skipped when inputs are missing) |
| `--refs-dir` | reference genomes named `<taxid>.fna` (default `<outdir>/refs`; missing ones are downloaded from NCBI, best effort) |
| `--score-min-reads`, `--score-top` | only align and score species with at least this many reads (default 100), at most this many per sample (default 50) |
| `--deepac-groups`, `--deepac-mode`, `--deepac-max-reads`, `--no-deepac` | DeePaC settings (default Bacteria, `rapid`, 50000 reads per species) |
| `--no-coverage` | skip alignment, depth and breadth |
| `--contigs`, `--assemble`, `--contig-kraken`, `--kraken-db` | assembly input and species attribution for resistance genes |
| `--amr-min-identity`, `--amr-min-coverage` | AMRFinderPlus filters (default 90 and 80) |
| `--report-min-reads`, `--report-min-rpm` | reporting thresholds (default 50 reads, 10 RPM) |
| `--host-pct` | host percentage shown in the report |
| `--threads`, `--keep-work`, `--cache`, `--api-key` | threads; keep extracted reads and BAMs; lookup cache; NCBI API key (or `NCBI_API_KEY`) |

### `pathoannot report`: rebuild the HTML after changing thresholds

```bash
pathoannot report --final results/PT001_final.csv --amr results/PT001_amr.csv --report-min-rpm 5 -o results
```

### `pathoannot check`: self-test

```bash
pathoannot check --taxid 573 10298 --name "Candida albicans" --priority-list icmr_priority.csv
```

---

## Reference genomes

For every species aligned, the tool looks for `<refs-dir>/<taxid>.fna` (also `.fa`, `.fasta`, and `.gz`). If missing it tries to download a RefSeq
reference from NCBI (viruses via Entrez). This is best effort; for reliable runs, fetch references yourself:

```bash
datasets download genome taxon 573 --reference --include genome --filename kp.zip
unzip kp.zip && cp ncbi_dataset/data/*/*.fna refs/573.fna
```

Species without a reference are still annotated and scored; their `presence_tier` is `Not assessed`, so they land in **Review** rather than **Reportable**.

## Adding the ICMR (or any other) priority list

Make a CSV with a `taxon` column and, optionally, `aliases` (separated by `;`), `taxid` (matches that species or anything below it),
`list` (the name shown in the output), `priority` and `note`. A template is in `examples/icmr_priority_template.csv`.

```
list,priority,taxon,aliases,taxid,note
ICMR syndromic list v1.0,priority,Streptococcus pneumoniae,,,ARI/AES
ICMR syndromic list v1.0,priority,Cryptococcus neoformans,,,AES
```

The tool ships the WHO lists only. ICMR content is not included, so supply it from the source document.

## Output files

| File | Content |
|---|---|
| `<sample>_annotated.csv` | evidence and contamination columns (annotate) |
| `<sample>_final.csv` | everything: evidence, contamination, abundance, alignment, DeePaC, resistance, `report_tier` |
| `<sample>_amr.csv` | each resistance gene with class, identity, coverage, contig and attributed species |
| `<sample>_report.html` | the report |
| `ALL_samples_*.csv` | the same tables for all samples together |

Key columns of `_final.csv`:

| Column | Meaning |
|---|---|
| `species`, `TaxID`, `group` | species from Bracken; group is Virus, Bacteria, Archaea, Fungi or Eukaryota (non-fungal) |
| `new_est_reads`, `fraction_total_reads`, `rpm`, `reads_per_million` | abundance |
| `priority_level`, `priority_hits`, `ncbi_label`, `lit_label` (and their percentages) | the three class evidence lines |
| `pathogen_support`, `evidence_lines`, `final_call`, `pathogen_class`, `review_flag` | combined class |
| `ntc_reads`, `ntc_rpm`, `rpm_ratio_vs_ntc`, `known_contaminant_genus`, `contamination_flag` | contamination |
| `minimizers_distinct`, `minimizer_dup_ratio` | breadth proxy from Kraken2 (if available) |
| `reads_extracted`, `mapped_reads`, `mean_depth`, `breadth_pct`, `expected_breadth_pct`, `breadth_ratio`, `rpkm`, `presence_tier` | alignment confirmation |
| `deepac_n_reads`, `deepac_mean_score`, `deepac_pct_above_threshold` | ML read scores |
| `n_amr_genes`, `amr_genes`, `amr_classes` | resistance genes attributed to the species |
| `report_tier`, `report_reason` | Reportable / Review / Excluded / Not reported, and why |

## Limitations

- **WHO lists concern drug-resistant pathogens.** A match shows the species is of concern, not that this strain is resistant. The lists are species-level and cannot separate strains such as *Salmonella* Typhi; *eumycetoma* agents cannot be matched automatically.
- **Literature counts are crude.** Heavily studied species can look disease-dominated, viruses are often published under abbreviations, and recently renamed species may need an alias.
- **NCBI metadata is biased** toward sequenced disease isolates (about 1,000 records per species are sampled). For viruses, "human host" means the virus infects people, not how severe it is.
- **DeePaC** covers bacteria only and is not applied to contigs. Fungi and viruses are annotated but not ML-scored.
- **Alignment** uses one reference genome; divergent strains and accessory genes lower breadth.
- **Contamination checks** depend on a good negative control; the genus list is a starter list.
- **Presence is not infection.** Interpret with clinical findings, and remove host reads first.
- The NCBI and Europe PMC steps need internet access.

## Troubleshooting

| Problem | What to do |
|---|---|
| `group=Unknown` or zero records in `pathoannot check` | no access to `api.ncbi.nlm.nih.gov`, or NCBI changed its response format; report the printed output |
| `NOT FOUND` for a tool | install it (see Install) and activate the environment |
| `no reference genome for <taxid>` | put `<taxid>.fna` in the refs folder |
| Resistance genes all `unassigned` | give `--contig-kraken` or `--kraken-db` |
| `amrfinder` database error | run `amrfinder -u` once |
| DeePaC import errors about protobuf or h5py | `pip install "protobuf==3.20.3" "h5py<3.0"` |
| Slow literature step | first run only; results are cached. Use `--sources priority,ncbi` to skip it |

## Development

```bash
pip install -e ".[test]"
pytest
```

The tests use mocked NCBI and Europe PMC responses and mocked DeePaC, minimap2, Kraken2 and AMRFinderPlus calls, so they run offline.
