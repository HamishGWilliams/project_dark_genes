# Project Dark Genes annotation and candidate-tracking status

## Current phase

The `revise-homology-filtering` branch has been pushed successfully with large generated TSVs excluded from Git. The core functional annotation, master-table construction, dark-candidate extraction, genome linking, genome-linking multiplicity QC, candidate-level duplication/prioritisation, BUSCO-backed duplication validation, and BUSCO-backed figure refresh are complete.

The active phase is repeat/TE-overlap integration for the 4,531 dark candidates. The first repeat-overlap run completed, but it was diagnostic rather than final because the original script treated all GFF3 features as repeat intervals. The script has now been patched to restrict GFF3 input to repeat-like feature types by default, so the repeat-overlap step should be rerun before moving to methylation/DMR integration.

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
- First repeat-overlap summary/table are tracked, but should be considered diagnostic because they include non-repeat GFF3 feature types such as `downstream_region`, `cds`, and GMAP/gene features.

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
- [X] Run first repeat-overlap pass using `00_raw/combined_annotations.gff3`.
- [X] Confirm first pass produced 4,531 candidate rows.
- [X] Identify that the first pass was not TE-specific because the combined GFF3 contributed non-repeat features.
- [X] Patch `scripts/add_repeat_overlap_context.py` so GFF3 input is restricted to repeat-like feature types by default.
- [ ] Rerun repeat/TE overlap with the patched script.

Diagnostic first-pass repeat-overlap summary:

```text
Candidate rows read: 4531
Repeat records parsed: 2048120
Repeat records skipped: 1490
repeat_overlap: 4531
```

Reason this pass is diagnostic only: top labels included `downstream_region`, `cds`, and GMAP/gene-model features, indicating that all GFF3 features were treated as repeats rather than only repeat/TE features.

## Active issue: rerun repeat/TE-overlap using repeat-like feature filter

The patched script now defaults to these GFF3 feature types:

```text
dispersed_repeat,repeat_region,repeat,transposable_element,transposable_element_gene,mobile_genetic_element,ltr_retrotransposon,non_ltr_retrotransposon,retrotransposon,dna_transposon,helitron,sine,line,ltr,long_terminal_repeat,tandem_repeat,satellite_dna,simple_repeat,low_complexity_region
```

Use the same command as before; do not pass `--include-all-gff3-features`.

## Immediate next action

Rerun repeat-overlap context using the patched script:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

git pull

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

After running, inspect:

```bash
cat "${PROJECT_DIR}/09_repeat_overlap/equina_dark_candidates.repeat_overlap.summary.txt"
awk -F '\t' 'NR>1 {sum++} END {print sum}' "${PROJECT_DIR}/09_repeat_overlap/equina_dark_candidates.repeat_overlap.tsv"
```

The repeat-overlap table should still contain 4,531 candidate rows, but overlap counts should now reflect only repeat/TE-like features.

## Remaining biological validation steps

- [ ] Rerun repeat/TE-overlap context with patched feature filtering.
- [ ] Add methylation/DMR overlap context.
- [ ] Add RNA-seq expression evidence.
- [ ] Flag expressed versus unsupported dark candidates.
- [ ] Prioritise stress-responsive dark candidates.
- [ ] Produce final candidate shortlist with genome, duplication, BUSCO, repeat/TE, expression, methylation, and annotation-evidence fields.
