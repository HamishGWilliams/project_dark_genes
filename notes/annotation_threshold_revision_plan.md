# Project Dark Genes annotation and candidate-tracking status

## Current phase

The `revise-homology-filtering` branch has been pushed successfully with large generated TSVs excluded from Git. The core functional annotation, master-table construction, dark-candidate extraction, genome linking, genome-linking multiplicity QC, candidate-level duplication/prioritisation, BUSCO-backed duplication validation, BUSCO-backed figure refresh, filtered repeat/TE-overlap integration, two-experiment methylation/DMR overlap integration, RNA-seq differential-expression integration, stress-responsive dark-gene plotting, and final all-stress candidate synthesis are complete.

The active phase is now review and interpretation of the final all-stress-responsive dark-gene shortlist.

## Repository status

- Branch: `revise-homology-filtering`
- Remote branch visible on GitHub: yes
- Large generated full TSVs excluded from pushed diff: yes
- Stress-responsive plot set is tracked under `12_final_candidates/stress_responsive_plots/`.
- Final synthesis script is tracked: `scripts/build_final_dark_candidate_shortlist.py`.
- Final synthesis script defaults to `--focal-contrasts all`, meaning every contrast present in `11_expression_context/equina_dark_candidates.de_significant_long.tsv` is treated as stress-response evidence.
- Final all-stress synthesis outputs are tracked under `12_final_candidates/`.

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
wide_pairs_detected: 10
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

Note: the plotting summary still labels 554 as `focal stress-responsive` because the plotting workflow used the earlier focal subset for signature visualisation. The final synthesis below is broader and uses all stress contrasts.

### Final all-stress candidate synthesis

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

Important interpretation: under the all-stress-contrast scope, 764 dark candidates are stress-responsive by the final synthesis scoring categories, including 684 top stress-responsive dark candidates. The top 250 are written to the final shortlist table.

Final outputs:

```text
12_final_candidates/equina_dark_candidates.final_integrated.tsv
12_final_candidates/equina_dark_candidates.final_shortlist.tsv
12_final_candidates/equina_dark_candidates.final_summary.txt
```

## Active issue: review and interpretation

The computational integration is now complete. The next work is interpretive and reporting-focused:

1. Review the top 250 final shortlist.
2. Decide whether to present a smaller top 20, top 50, or top 100 subset in the written results.
3. Produce a final written interpretation of the all-stress-responsive dark-gene set.
4. Optionally regenerate stress-responsive plots so the plotting workflow also uses all stress contrasts rather than the earlier focal subset.

## Suggested next commands

Inspect the final all-stress summary and shortlist:

```bash
PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes

cat "${PROJECT_DIR}/12_final_candidates/equina_dark_candidates.final_summary.txt"

cut -f1-12 "${PROJECT_DIR}/12_final_candidates/equina_dark_candidates.final_shortlist.tsv" | head -n 30
```

If desired, regenerate plots with all stress contrasts by setting:

```bash
FOCAL_CONTRASTS="combined_wald,diesel_added_wald,full_model_LRT,interactive_only_LRT,interactive_only_wald,salinity_added_wald,salinity_only_LRT,salinity_only_wald" \
PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes \
Rscript scripts/make_stress_responsive_dark_gene_plots.R
```

## Remaining steps

- [ ] Review top 250 shortlist.
- [ ] Decide the final reporting subset size.
- [ ] Optionally regenerate stress-responsive plots using all stress contrasts.
- [ ] Produce a final written interpretation of the all-stress-responsive dark-gene set.
