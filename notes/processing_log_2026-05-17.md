# Project Dark Genes processing log — 2026-05-17

This log records the commands, scripts, and processing decisions used during the 2026-05-17 Project Dark Genes session. It is intended as a reproducibility aid rather than a full methods section.

## Branch

```bash
git checkout revise-homology-filtering
```

Primary working branch:

```text
revise-homology-filtering
```

Primary scratch project directory used by pipeline commands:

```bash
PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes
```

Repository working copy on Maxwell:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes
```

## 1. Large generated files and `.gitignore` clean-up

Large generated TSVs were removed from Git tracking and kept local. The important pattern was:

```bash
git rm --cached --ignore-unmatch <large_generated_file.tsv>
git rm --cached --ignore-unmatch <large_generated_directory>/*.tsv
```

Large/generated files were excluded from future commits with `.gitignore` entries for:

```text
03_dark_candidates/*.tsv
03_dark_candidates/*.fa
03_dark_candidates/*.faa
03_dark_candidates/*.fasta
03_dark_candidates/genomic_context/*.tsv
08_figures/figure_tables/sanitised_inputs/*.tsv
05_busco/**
```

Concise summaries and compact figure tables were allowed to remain tracked where appropriate.

## 2. Genome-linking multiplicity QC

Problem identified:

```text
Candidate rows read: 4531
Linked output rows: 574504
Unmatched candidates: 0
matched_multiple: 4531
```

Script added:

```text
scripts/qc_genome_linking_multiplicity.py
```

Command used:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

python3 scripts/qc_genome_linking_multiplicity.py \
  --linked /uoa/scratch/users/r02hw22/project_dark_genes/06_genome_lookup/equina_dark_candidates.genome_linked.tsv \
  --outdir /uoa/scratch/users/r02hw22/project_dark_genes/06_genome_lookup \
  --candidate-id-column protein_id \
  --prefix equina_dark_candidates
```

Expected/confirmed output:

```text
06_genome_lookup/equina_dark_candidates.genome_linked.primary.tsv
06_genome_lookup/equina_dark_candidates.genome_linking_multiplicity.tsv
06_genome_lookup/equina_dark_candidates.genome_linking_multiplicity.summary.txt
```

Confirmed summary:

```text
Total linked rows read: 574504
Unique candidate IDs: 4531
Primary rows written: 4531
matched_primary_from_multiple_strict: 4531
rank_1: 4531
```

## 3. Candidate-level duplication/prioritisation rerun

The old prioritised table had been generated from the 574,504-row genome-linked table. It was regenerated from the primary one-row-per-candidate table.

Representative command pattern:

```bash
PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes

PRIMARY_LINKED="${PROJECT_DIR}/06_genome_lookup/equina_dark_candidates.genome_linked.primary.tsv"
LOOKUP="${PROJECT_DIR}/06_genome_lookup/equina_gff3_gene_transcript_protein_contig_lookup.tsv"
DARK_FASTA="${PROJECT_DIR}/03_dark_candidates/equina_dark_candidates.all.fa"
DUP_OUTDIR="${PROJECT_DIR}/07_duplication_context"

python3 scripts/add_busco_duplication_context.py \
  --linked "${PRIMARY_LINKED}" \
  --lookup "${LOOKUP}" \
  --dark-fasta "${DARK_FASTA}" \
  --busco-full-table "NA" \
  --cluster-file "${DUP_OUTDIR}/equina_dark_candidate_clusters.exact.tsv" \
  --cluster-method exact_python_fallback \
  --outdir "${DUP_OUTDIR}"
```

Initial no-BUSCO candidate-level result:

```text
high_priority: 4209
medium_priority: 322
low_priority_or_manual_review: 0
```

## 4. BUSCO-backed duplication validation

Reason for this step: the previous duplication/prioritisation run used `BUSCO full table: NA`, while prior BUSCO interpretation showed a high duplicated BUSCO fraction.

Script added:

```text
scripts/run_busco_duplication_validation.sh
```

Command used:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes
sbatch scripts/run_busco_duplication_validation.sh
```

If an existing BUSCO full table needed to be forced:

```bash
BUSCO_FULL_TABLE=/path/to/full_table.tsv \
sbatch scripts/run_busco_duplication_validation.sh
```

Stable BUSCO output path used:

```text
/uoa/scratch/users/r02hw22/project_dark_genes/05_busco/equina_representative_longest_per_gene_metazoa_odb10/full_table.tsv
```

Confirmed BUSCO-backed summary values:

```text
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

## 5. Figure regeneration from BUSCO-backed prioritisation

Figures were regenerated after the BUSCO-backed prioritisation step.

Command used:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes

GENOME_LINKED_TSV="${PROJECT_DIR}/06_genome_lookup/equina_dark_candidates.genome_linked.primary.tsv" \
PRIORITISED_TSV="${PROJECT_DIR}/07_duplication_context/equina_dark_candidates.prioritised.tsv" \
BUSCO_FULL_TABLE="${PROJECT_DIR}/05_busco/equina_representative_longest_per_gene_metazoa_odb10/full_table.tsv" \
bash scripts/make_figures.sh
```

Key checks:

```bash
cat "${PROJECT_DIR}/08_figures/figure_tables/figure_13_duplication_interpretation_counts.tsv"
cat "${PROJECT_DIR}/08_figures/figure_tables/figure_14_priority_tier_counts.tsv"
```

Confirmed Figure 13:

```text
single_copy_or_no_near_identical_dark_duplicate       2926
possible_biological_gene_family_expansion             1385
possible_assembly_redundancy_or_haplotig_duplication   220
```

Confirmed Figure 14:

```text
high_priority                    3903
medium_priority                   420
low_priority_or_manual_review     208
```

## 6. Repeat/TE-overlap integration

Script added:

```text
scripts/add_repeat_overlap_context.py
```

Initial diagnostic run used all GFF3 features and was rejected as too permissive because labels included ordinary gene features such as `downstream_region`, `cds`, and GMAP/gene-model features.

The script was patched to restrict GFF3 input to repeat-like feature types by default.

Final command used:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes

python3 scripts/add_repeat_overlap_context.py \
  --candidates "${PROJECT_DIR}/07_duplication_context/equina_dark_candidates.prioritised.tsv" \
  --repeats "${PROJECT_DIR}/00_raw/combined_annotations.gff3" \
  --outdir "${PROJECT_DIR}/09_repeat_overlap" \
  --prefix equina_dark_candidates \
  --candidate-id-column protein_id \
  --contig-column contig \
  --start-column transcript_start \
  --end-column transcript_end
```

Confirmed filtered repeat/TE result:

```text
Candidate rows read: 4531
GFF3 feature filtering enabled: True
parsed_gff3: 978221
gff3_feature_type_filtered: 1069899
comment_or_blank: 1490
repeat/TE contigs indexed: 1486
repeat_overlap: 3329
no_repeat_overlap: 1202
```

Output files:

```text
09_repeat_overlap/equina_dark_candidates.repeat_overlap.tsv
09_repeat_overlap/equina_dark_candidates.repeat_overlap.summary.txt
```

## 7. DMR overlap integration

Two methylation/DMR experiments were added:

```text
exp1_acute_naive_diesel: acute/naive response to diesel exposure
exp2_primed_acclimated_diesel: primed/acclimated response to diesel exposure
```

Input directory:

```text
/uoa/scratch/users/r02hw22/project_dark_genes/10_methylation_overlap
```

Scripts added:

```text
scripts/add_dmr_overlap_context.py
scripts/run_dmr_overlap_experiments.sh
```

The DMR parser was patched to handle:

```text
seqnames,start,end
ID
meth.diff
comma/tab delimiter auto-detection
quoted fields
```

Wrapper command used:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes
sbatch scripts/run_dmr_overlap_experiments.sh
```

Explicit-path pattern if needed:

```bash
EXP1_DMR=/uoa/scratch/users/r02hw22/project_dark_genes/10_methylation_overlap/exp1_DMRs.txt \
EXP2_DMR=/uoa/scratch/users/r02hw22/project_dark_genes/10_methylation_overlap/exp2_DMRs.txt \
sbatch scripts/run_dmr_overlap_experiments.sh
```

Confirmed exp1 result:

```text
Candidate rows read: 4531
dmr_format_used: delimited_table
dmr_contig_column: seqnames
dmr_start_column: start
dmr_end_column: end
dmr_id_column: ID
dmr_effect_column: meth.diff
parsed_table: 16590
DMR contigs indexed: 16590
no_dmr_overlap: 4531
```

Confirmed exp2 result:

```text
Candidate rows read: 4531
dmr_format_used: delimited_table
dmr_contig_column: seqnames
dmr_start_column: start
dmr_end_column: end
dmr_id_column: ID
dmr_effect_column: meth.diff
parsed_table: 22813
DMR contigs indexed: 22813
no_dmr_overlap: 4531
```

Interpretation:

```text
No dark candidates overlap diesel-response DMR intervals in either experiment.
```

Output files:

```text
10_methylation_overlap/equina_dark_candidates.dmr_overlap.experiment_manifest.tsv
10_methylation_overlap/exp1_acute_naive_diesel/equina_dark_candidates.exp1_acute_naive_diesel.dmr_overlap.tsv
10_methylation_overlap/exp1_acute_naive_diesel/equina_dark_candidates.exp1_acute_naive_diesel.dmr_overlap.summary.txt
10_methylation_overlap/exp2_primed_acclimated_diesel/equina_dark_candidates.exp2_primed_acclimated_diesel.dmr_overlap.tsv
10_methylation_overlap/exp2_primed_acclimated_diesel/equina_dark_candidates.exp2_primed_acclimated_diesel.dmr_overlap.summary.txt
```

## 8. RNA-seq differential-expression integration

Input DE file:

```text
/uoa/scratch/users/r02hw22/project_dark_genes/05_rnaseq/05_de/Multi_Stressor_all_Differential_Expression_Analysis_results.csv
```

Script added:

```text
scripts/add_de_expression_context.py
```

Command used:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes

python3 scripts/add_de_expression_context.py \
  --candidates "${PROJECT_DIR}/10_methylation_overlap/exp1_acute_naive_diesel/equina_dark_candidates.exp1_acute_naive_diesel.dmr_overlap.tsv" \
  --de "${PROJECT_DIR}/05_rnaseq/05_de/Multi_Stressor_all_Differential_Expression_Analysis_results.csv" \
  --outdir "${PROJECT_DIR}/11_expression_context" \
  --prefix equina_dark_candidates \
  --padj-threshold 0.05 \
  --lfc-threshold 1
```

Confirmed DE integration summary:

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

Output files:

```text
11_expression_context/equina_dark_candidates.de_context.tsv
11_expression_context/equina_dark_candidates.de_significant_long.tsv
11_expression_context/equina_dark_candidates.de_context.summary.txt
```

## 9. Stress-responsive dark-gene plotting

Script added and then adjusted for smaller text, wrapped labels, tighter margins, and wider export canvases:

```text
scripts/make_stress_responsive_dark_gene_plots.R
```

Command used:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes \
Rscript scripts/make_stress_responsive_dark_gene_plots.R
```

The plotting script applied `bbplot::bbc_style()` when available, with a fallback clean minimal theme.

Confirmed plotting summary:

```text
Candidates read: 4531
Significant DE rows read: 1723
Unique significant dark candidates: 726
Unique focal stress-responsive dark candidates: 554
Diesel-added significant dark candidates: 69
```

Generated plot set:

```text
12_final_candidates/stress_responsive_plots/stress_figure_16_de_candidates_by_contrast.png
12_final_candidates/stress_responsive_plots/stress_figure_17_de_direction_by_contrast.png
12_final_candidates/stress_responsive_plots/stress_figure_18_priority_vs_focal_stress.png
12_final_candidates/stress_responsive_plots/stress_figure_19_repeat_context_vs_focal_stress.png
12_final_candidates/stress_responsive_plots/stress_figure_20_focal_stress_signature_combinations.png
12_final_candidates/stress_responsive_plots/stress_figure_21_diesel_added_dark_gene_response.png
12_final_candidates/stress_responsive_plots/stress_figure_22_dark_gene_evidence_stack.png
```

Note: these plots used the earlier focal stress set for signature visualisation. The final synthesis below uses all stress contrasts.

## 10. Final all-stress candidate synthesis

Script added:

```text
scripts/build_final_dark_candidate_shortlist.py
```

The script was patched so the default is:

```text
--focal-contrasts all
```

This means every contrast present in the significant-long DE table is treated as stress-response evidence.

Command used for the final all-stress run:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes

python3 scripts/build_final_dark_candidate_shortlist.py \
  --de-context "${PROJECT_DIR}/11_expression_context/equina_dark_candidates.de_context.tsv" \
  --sig-long "${PROJECT_DIR}/11_expression_context/equina_dark_candidates.de_significant_long.tsv" \
  --outdir "${PROJECT_DIR}/12_final_candidates" \
  --prefix equina_dark_candidates \
  --top-n 250
```

Confirmed final all-stress synthesis summary:

```text
Candidate rows read: 4531
Significant-DE candidate IDs indexed: 726
Shortlist rows written: 250
Stress contrast scope: all
Stress contrasts used: combined_wald,diesel_added_wald,full_model_LRT,interactive_only_LRT,interactive_only_wald,salinity_added_wald,salinity_only_LRT,salinity_only_wald
```

Final shortlist category counts:

```text
retain_context_only: 3767
top_stress_responsive_dark_candidate: 684
stress_responsive_dark_candidate: 67
strong_stress_responsive_dark_candidate: 13
```

Most common shortlist evidence reasons:

```text
no_DMR_overlap: 4531
high_priority: 3903
repeat_TE_overlap_context: 3329
single_copy_or_no_near_identical_dark_duplicate: 2926
possible_biological_gene_family_expansion: 1385
no_repeat_TE_overlap: 1202
DE_significant_any_contrast: 764
stress_contrast_DE: 764
medium_priority: 420
abs_logFC_ge_2: 283
padj_le_0.001: 235
padj_le_0.01: 227
assembly_redundancy_or_haplotig_caution: 220
low_priority_or_manual_review: 208
abs_logFC_ge_1.5: 193
diesel_added_DE: 69
```

Final outputs:

```text
12_final_candidates/equina_dark_candidates.final_integrated.tsv
12_final_candidates/equina_dark_candidates.final_shortlist.tsv
12_final_candidates/equina_dark_candidates.final_summary.txt
```

## 11. Tracker updates

Main tracker updated repeatedly during the session:

```text
notes/annotation_threshold_revision_plan.md
```

Current status at end of this log:

```text
Functional annotation complete.
Dark-candidate extraction complete.
Genome-linking and multiplicity QC complete.
BUSCO-backed duplication validation complete.
Repeat/TE-overlap integration complete.
DMR integration complete with valid negative overlap result.
RNA-seq DE integration complete.
Stress-responsive dark-gene plots complete.
Final all-stress candidate synthesis complete.
Next step: review and interpretation of the final all-stress-responsive dark-gene shortlist.
```

## 12. Next planned step

Review the final top 250 shortlist and decide the reporting subset size:

```bash
PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes

cat "${PROJECT_DIR}/12_final_candidates/equina_dark_candidates.final_summary.txt"
cut -f1-12 "${PROJECT_DIR}/12_final_candidates/equina_dark_candidates.final_shortlist.tsv" | head -n 30
```

Possible reporting subsets:

```text
top 20 all-stress-responsive dark candidates
top 50 all-stress-responsive dark candidates
top 100 all-stress-responsive dark candidates
full top 250 supplement
```

Optional next analysis: regenerate the stress-responsive plots with all stress contrasts rather than the earlier focal subset.
