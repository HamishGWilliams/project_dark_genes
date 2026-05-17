# Project Dark Genes annotation and candidate-tracking status

## Current phase

The `revise-homology-filtering` branch has been pushed successfully with large generated TSVs excluded from Git. The core functional annotation, master-table construction, dark-candidate extraction, genome linking, genome-linking multiplicity QC, candidate-level duplication/prioritisation, BUSCO-backed duplication validation, BUSCO-backed figure refresh, filtered repeat/TE-overlap integration, two-experiment methylation/DMR overlap integration, RNA-seq differential-expression integration, and stress-responsive dark-gene plotting are complete.

The active phase is now final candidate synthesis: combine annotation darkness, genome context, BUSCO/duplication validation, repeat/TE overlap, DMR results, RNA-seq DE evidence, and stress-responsive plotting outputs into a final shortlist and interpretation table.

## Repository status

- Branch: `revise-homology-filtering`
- Remote branch visible on GitHub: yes
- Large generated full TSVs excluded from pushed diff: yes
- Full generated genome-linked table remains local/generated and ignored: `03_dark_candidates/genomic_context/equina_dark_candidates.genome_linked.tsv`
- First-pass genome-linking summary is tracked.
- Multiplicity-QC summary is tracked.
- BUSCO-backed duplication/prioritisation summary is tracked.
- BUSCO-backed validation wrapper is tracked: `scripts/run_busco_duplication_validation.sh`.
- BUSCO-backed compact figure tables and figures are tracked.
- Repeat-overlap workflow script is tracked: `scripts/add_repeat_overlap_context.py`.
- Filtered repeat/TE-overlap table and summary are tracked.
- DMR overlap script is tracked: `scripts/add_dmr_overlap_context.py`.
- Two-experiment DMR wrapper is tracked: `scripts/run_dmr_overlap_experiments.sh`.
- Two-experiment DMR overlap summaries, tables, and manifest are tracked.
- RNA-seq DE integration script and outputs are tracked: `scripts/add_de_expression_context.py` and `11_expression_context/`.
- Stress-responsive dark-gene plots, plot tables, manifest, and plotting script are tracked under `12_final_candidates/stress_responsive_plots/` and `scripts/make_stress_responsive_dark_gene_plots.R`.
- Final shortlist synthesis script is tracked: `scripts/build_final_dark_candidate_shortlist.py`.

The branch is currently diverged from `main` because `main` contains a small `.gitignore` update made while the branch was being cleaned. This can be resolved by merging or rebasing after the current branch state is stable.

## Key confirmed results

### Dark-candidate extraction

```text
Total input rows scanned: 47671
Total FASTA records indexed: 47671
Total dark candidates: 4531
function_dark_but_signalp_secretory_candidate: 0
function_dark_no_current_annotation: 4531
Candidates with retained review evidence: 4531
Evidence columns retained for manual review: 78
FASTA records written: 4531
Candidates missing sequences: 0
```

### Genome linking and multiplicity QC

```text
First-pass linked rows: 574504
Unique candidate IDs after primary collapse: 4531
Primary rows written: 4531
matched_primary_from_multiple_strict: 4531
rank_1: 4531
```

### BUSCO-backed duplication validation

```text
BUSCO full table: /uoa/scratch/users/r02hw22/project_dark_genes/05_busco/equina_representative_longest_per_gene_metazoa_odb10/full_table.tsv
Total BUSCO rows parsed: 1497
Duplicated BUSCO loci parsed: 950
BUSCO-duplication-rich scaffolds: 74
single_copy_or_no_near_identical_dark_duplicate: 2926
possible_biological_gene_family_expansion: 1385
possible_assembly_redundancy_or_haplotig_duplication: 220
high_priority: 3903
medium_priority: 420
low_priority_or_manual_review: 208
```

### Repeat/TE overlap

```text
Candidate rows read: 4531
GFF3 feature filtering enabled: True
parsed_gff3: 978221
gff3_feature_type_filtered: 1069899
repeat/TE contigs indexed: 1486
repeat_overlap: 3329
no_repeat_overlap: 1202
```

### Methylation/DMR overlap

```text
exp1_acute_naive_diesel:
parsed_table: 16590
DMR contigs indexed: 16590
no_dmr_overlap: 4531

exp2_primed_acclimated_diesel:
parsed_table: 22813
DMR contigs indexed: 22813
no_dmr_overlap: 4531
```

Interpretation: after parser correction, this is a valid negative result. None of the 4,531 dark candidates overlap the acute/naive or primed/acclimated diesel-response DMR intervals.

### RNA-seq differential-expression integration

```text
Candidate rows read: 4531
DE rows read: 42101
DE IDs indexed: 42091
padj threshold: 0.05
abs(logFC) threshold: 1.0
de_id_column: gene
de_format_used: wide
wide_pairs_detected: 10
parsed_de_record: 394453
missing_lfc_or_padj: 26557
de_record_matched_not_significant: 2263
no_de_record_matched: 1542
de_significant: 726
```

Significant dark-candidate counts by contrast:

```text
combined_wald: 406
salinity_added_wald: 398
full_model_LRT: 301
salinity_only_wald: 270
salinity_only_LRT: 269
diesel_added_wald: 69
interactive_only_wald: 7
interactive_only_LRT: 3
```

### Stress-responsive dark-gene plotting

Plots are tracked under:

```text
12_final_candidates/stress_responsive_plots/
```

Confirmed plot summary:

```text
Candidates read: 4531
Significant DE rows read: 1723
Unique significant dark candidates: 726
Unique focal stress-responsive dark candidates: 554
Diesel-added significant dark candidates: 69
```

Generated plot set:

```text
stress_figure_16_de_candidates_by_contrast
stress_figure_17_de_direction_by_contrast
stress_figure_18_priority_vs_focal_stress
stress_figure_19_repeat_context_vs_focal_stress
stress_figure_20_focal_stress_signature_combinations
stress_figure_21_diesel_added_dark_gene_response
stress_figure_22_dark_gene_evidence_stack
```

The plotting script was adjusted after the first pass to reduce text size, wrap long labels, tighten margins, and widen export canvases so the figures fit the plotting space more cleanly.

## Active issue: final candidate synthesis and shortlist

A final synthesis script has been added:

```text
scripts/build_final_dark_candidate_shortlist.py
```

Recommended input:

```text
11_expression_context/equina_dark_candidates.de_context.tsv
11_expression_context/equina_dark_candidates.de_significant_long.tsv
```

Recommended outputs:

```text
12_final_candidates/equina_dark_candidates.final_integrated.tsv
12_final_candidates/equina_dark_candidates.final_shortlist.tsv
12_final_candidates/equina_dark_candidates.final_summary.txt
```

## Immediate next action

Run the final shortlist synthesis:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

git pull

PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes

python3 scripts/build_final_dark_candidate_shortlist.py \
  --de-context "${PROJECT_DIR}/11_expression_context/equina_dark_candidates.de_context.tsv" \
  --sig-long "${PROJECT_DIR}/11_expression_context/equina_dark_candidates.de_significant_long.tsv" \
  --outdir "${PROJECT_DIR}/12_final_candidates" \
  --prefix equina_dark_candidates \
  --top-n 250
```

Then inspect:

```bash
cat "${PROJECT_DIR}/12_final_candidates/equina_dark_candidates.final_summary.txt"
head -n 20 "${PROJECT_DIR}/12_final_candidates/equina_dark_candidates.final_shortlist.tsv"
```

## Remaining steps

- [ ] Run final integrated candidate synthesis.
- [ ] Review top 250 shortlist.
- [ ] Decide whether the final shortlist should be narrowed to diesel-only candidates, focal-stress candidates, or all multi-stressor dark candidates.
- [ ] Produce a final written interpretation of the stress-responsive dark-gene set.
