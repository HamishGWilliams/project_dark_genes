# Project Dark Genes annotation and candidate-tracking status

## Current phase

The `revise-homology-filtering` branch has been pushed successfully with large generated TSVs excluded from Git. The core functional annotation, master-table construction, dark-candidate extraction, genome linking, genome-linking multiplicity QC, candidate-level duplication/prioritisation, BUSCO-backed duplication validation, BUSCO-backed figure refresh, filtered repeat/TE-overlap integration, two-experiment methylation/DMR overlap integration, and RNA-seq differential-expression integration are complete.

The active phase is now final candidate synthesis: combine annotation darkness, genome context, BUSCO/duplication validation, repeat/TE overlap, DMR results, and RNA-seq DE evidence into a final shortlist and interpretation table.

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

The branch is currently diverged from `main` because `main` contains a small `.gitignore` update made while the branch was being cleaned. This can be resolved by merging or rebasing after the current branch state is stable.

## Annotation threshold decision

DIAMOND and BLASTp evidence is filtered before assignment to the sequence-supported annotation classes. The project uses:

```text
e-value <= 1e-5
```

Identity, alignment length, bitscore, query coverage, and subject coverage are retained for QC/reporting rather than being used as hard filters for the main six annotation classes.

## Annotation hierarchy

The master table assigns one primary annotation class per representative protein using this fixed hierarchy:

1. `annotated_swissprot_supported`
2. `sequence_supported_trembl_cnidaria`
3. `domain_supported_interpro`
4. `orthology_supported_eggnog`
5. `function_dark_but_signalp_secretory_candidate`
6. `function_dark_no_current_annotation`

Only threshold-passing sequence hits are allowed to assign categories 1 or 2.

## Completed checkpoints

### Homology filtering and annotation evidence

- [X] Filter DIAMOND Swiss-Prot evidence at `e-value <= 1e-5`.
- [X] Filter DIAMOND Cnidaria-TrEMBL evidence at `e-value <= 1e-5`.
- [X] Filter BLASTp Swiss-Prot evidence at `e-value <= 1e-5`.
- [X] Filter BLASTp Cnidaria-TrEMBL evidence at `e-value <= 1e-5`.
- [X] Generate filtered top-hit/best-hit files after e-value filtering.
- [X] Confirm full InterProScan output exists.
- [X] Confirm full eggNOG output exists.
- [X] Confirm full SignalP output exists.

### Master annotation and dark-candidate extraction

- [X] Build master annotation script.
- [X] Rebuild full representative master annotation from filtered evidence.
- [X] Run full master-table QC.
- [X] Extract dark candidates from the full master table.
- [X] Generate dark-candidate TSV outputs.
- [X] Generate dark-candidate FASTA outputs.
- [X] Generate dark-candidate extraction summary.

Confirmed extraction summary:

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

- [X] Run first-pass GFF3/genome lookup.
- [X] Generate first-pass genome-linked output locally.
- [X] Track first-pass genome-linking summary.
- [X] Add `scripts/qc_genome_linking_multiplicity.py`.
- [X] Run genome-linking multiplicity QC.
- [X] Track the multiplicity-QC summary.
- [X] Produce local one-row-per-candidate primary mapping table.
- [X] Confirm primary table contains 4,531 dark-candidate rows.

First-pass genome-linking summary:

```text
Candidate rows read: 4531
Linked output rows: 574504
Unmatched candidates: 0
matched_multiple: 4531
```

Confirmed multiplicity-QC summary:

```text
Total linked rows read: 574504
Unique candidate IDs: 4531
Primary rows written: 4531
matched_primary_from_multiple_strict: 4531
rank_1: 4531
```

### Corrected candidate-level duplication/prioritisation and figures

- [X] Rerun duplication/prioritisation using the primary genome-linked table.
- [X] Regenerate Figures 05–15 using candidate-level inputs.
- [X] Confirm Figure 05 sums to 4,531 genome-linked candidates.
- [X] Confirm Figure 14 priority tiers sum to 4,531 candidates in the no-BUSCO candidate-level run.
- [X] Track corrected duplication/prioritisation summary.
- [X] Track corrected compact figure tables and figures.

Previous no-BUSCO candidate-level counts:

```text
Figure 05 genome-linking status:
matched_multiple: 4531

Figure 14 priority tiers:
high_priority: 4209
medium_priority: 322
low_priority_or_manual_review: 0
```

### BUSCO-backed duplication validation and figure refresh

- [X] Add `scripts/run_busco_duplication_validation.sh`.
- [X] Run BUSCO-backed duplication validation.
- [X] Confirm duplication/prioritisation summary now uses a real BUSCO full table.
- [X] Confirm BUSCO duplicated loci and duplication-rich scaffolds are incorporated in the duplication summary.
- [X] Regenerate and push compact Figures 13–15 from the BUSCO-backed prioritised table.

Confirmed BUSCO-backed duplication/prioritisation summary:

```text
Genome-linked candidate TSV: /uoa/scratch/users/r02hw22/project_dark_genes/06_genome_lookup/equina_dark_candidates.genome_linked.primary.tsv
BUSCO full table: /uoa/scratch/users/r02hw22/project_dark_genes/05_busco/equina_representative_longest_per_gene_metazoa_odb10/full_table.tsv
Total BUSCO rows parsed: 1497
Duplicated BUSCO loci parsed: 950
BUSCO-duplication-rich scaffolds: 74
Cluster rows written: 4531
Unique clusters: 3371
large_near_identical_clusters_5_plus: 104
singleton_clusters: 2926
small_near_identical_clusters_2_to_4: 341
single_copy_or_no_near_identical_dark_duplicate: 2926
possible_biological_gene_family_expansion: 1385
possible_assembly_redundancy_or_haplotig_duplication: 220
high_priority: 3903
medium_priority: 420
low_priority_or_manual_review: 208
```

Confirmed BUSCO-backed figure tables:

```text
Figure 13 duplication interpretation:
single_copy_or_no_near_identical_dark_duplicate: 2926
possible_biological_gene_family_expansion: 1385
possible_assembly_redundancy_or_haplotig_duplication: 220

Figure 14 priority tiers:
high_priority: 3903
medium_priority: 420
low_priority_or_manual_review: 208
```

Interpretation: BUSCO validation moved 220 candidates into an assembly-redundancy/haplotig-duplication interpretation and reduced the high-priority candidate set from 4,209 to 3,903.

### Repeat/TE-overlap integration

- [X] Add `scripts/add_repeat_overlap_context.py`.
- [X] Run first diagnostic repeat-overlap pass using `00_raw/combined_annotations.gff3`.
- [X] Confirm diagnostic pass produced 4,531 candidate rows.
- [X] Identify that the diagnostic pass was not TE-specific because the combined GFF3 contributed non-repeat features.
- [X] Patch `scripts/add_repeat_overlap_context.py` so GFF3 input is restricted to repeat-like feature types by default.
- [X] Rerun repeat/TE overlap with the patched script.
- [X] Confirm filtered repeat/TE-overlap summary and candidate table are tracked.

Confirmed filtered repeat/TE-overlap summary:

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

Interpretation: after restricting to repeat-like GFF3 feature types, 3,329 dark candidates overlap repeat/TE-like intervals and 1,202 do not. The previous all-features result should be treated as diagnostic only.

### Methylation/DMR overlap integration

- [X] Add `scripts/add_dmr_overlap_context.py`.
- [X] Add `scripts/run_dmr_overlap_experiments.sh`.
- [X] Run DMR overlap for `exp1_acute_naive_diesel`.
- [X] Run DMR overlap for `exp2_primed_acclimated_diesel`.
- [X] Confirm DMR parser successfully indexed intervals for both experiments.
- [X] Confirm no dark candidates overlap DMR intervals in either diesel-response experiment.

Confirmed exp1 acute/naive diesel DMR-overlap summary:

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

Confirmed exp2 primed/acclimated diesel DMR-overlap summary:

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

Interpretation: after parser correction, this is a valid negative result. None of the 4,531 dark candidates overlap the acute/naive or primed/acclimated diesel-response DMR intervals.

### RNA-seq differential-expression integration

- [X] Add `scripts/add_de_expression_context.py`.
- [X] Run DE integration against `05_rnaseq/05_de/Multi_Stressor_all_Differential_Expression_Analysis_results.csv`.
- [X] Confirm DE table parsing and contrast-pair detection.
- [X] Confirm DE context outputs are tracked under `11_expression_context/`.

Confirmed DE-context summary:

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

Interpretation: 2,989 of 4,531 dark candidates matched RNA-seq DE records, and 726 were differentially expressed in at least one multi-stressor contrast. The strongest dark-candidate DE signal is not the pure diesel contrast but the combined/salinity-associated contrasts; `diesel_added_wald` identifies 69 significant dark candidates.

## Active issue: final candidate synthesis and shortlist

Next, generate a final candidate synthesis table that combines:

```text
annotation darkness
BUSCO-backed duplication interpretation
priority tier
repeat/TE overlap
DMR overlap from exp1 and exp2
RNA-seq DE status and significant contrasts
```

Recommended input:

```text
11_expression_context/equina_dark_candidates.de_context.tsv
```

Recommended outputs:

```text
12_final_candidates/equina_dark_candidates.final_integrated.tsv
12_final_candidates/equina_dark_candidates.final_shortlist.tsv
12_final_candidates/equina_dark_candidates.final_summary.txt
```

## Remaining biological validation steps

- [ ] Generate final integrated candidate table.
- [ ] Generate final shortlist of strongest dark candidates.
- [ ] Prioritise candidates with high priority, non-assembly-redundant duplication status, no DMR overlap caveat, and significant diesel/combined stress DE evidence.
- [ ] Produce final candidate shortlist with genome, duplication, BUSCO, repeat/TE, DMR, expression, diesel-response, methylation, and annotation-evidence fields.
