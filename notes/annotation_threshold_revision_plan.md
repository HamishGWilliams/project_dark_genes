# Project Dark Genes annotation and candidate-tracking status

## Current phase

The `revise-homology-filtering` branch has been pushed successfully with large generated TSVs excluded from Git. The core functional annotation, master-table construction, dark-candidate extraction, genome linking, genome-linking multiplicity QC, duplication/prioritisation rerun, and candidate-level figure regeneration are complete.

Before moving to repeat/TE-overlap integration, the active phase is now BUSCO-backed duplication validation. This is needed because the current duplication/prioritisation run used `BUSCO_FULL_TABLE=NA`, so priority tiers currently reflect genome linking and exact-sequence clustering, but not BUSCO-validated scaffold duplication.

## Repository status

- Branch: `revise-homology-filtering`
- Remote branch visible on GitHub: yes
- Large generated full TSVs excluded from pushed diff: yes
- Full generated genome-linked table remains local/generated and ignored: `03_dark_candidates/genomic_context/equina_dark_candidates.genome_linked.tsv`
- First-pass genome-linking summary is tracked.
- Multiplicity-QC summary is tracked.
- Duplication/prioritisation summary is tracked.
- Corrected figures, figure manifest, compact figure tables, scripts, and notes are tracked.
- BUSCO-backed validation wrapper has been added: `scripts/run_busco_duplication_validation.sh`.

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

### Corrected duplication/prioritisation and figures

- [X] Rerun duplication/prioritisation using the primary genome-linked table.
- [X] Regenerate Figures 05–15 using candidate-level inputs.
- [X] Confirm Figure 05 now sums to 4,531 genome-linked candidates.
- [X] Confirm Figure 14 priority tiers now sum to 4,531 candidates.
- [X] Track corrected duplication/prioritisation summary.
- [X] Track corrected compact figure tables and figures.

Confirmed corrected candidate-level figure counts:

```text
Figure 05 genome-linking status:
matched_multiple: 4531

Figure 14 priority tiers:
high_priority: 4209
medium_priority: 322
low_priority_or_manual_review: 0
```

Confirmed duplication/prioritisation summary:

```text
Genome-linked candidate TSV: /uoa/scratch/users/r02hw22/project_dark_genes/06_genome_lookup/equina_dark_candidates.genome_linked.primary.tsv
BUSCO full table: NA
Cluster method: exact_python_fallback
Cluster rows written: 4531
Unique clusters: 3371
large_near_identical_clusters_5_plus: 104
singleton_clusters: 2926
small_near_identical_clusters_2_to_4: 341
single_copy_or_no_near_identical_dark_duplicate: 2926
possible_biological_gene_family_expansion: 1605
high_priority: 4209
medium_priority: 322
```

Important caveat: BUSCO context is currently unavailable in this run because `BUSCO_FULL_TABLE` was `NA`. Duplication/prioritisation currently reflects exact-sequence clustering and genome context, not BUSCO-validated scaffold duplication.

## Active issue: BUSCO-backed duplication validation

A new wrapper has been added:

```text
scripts/run_busco_duplication_validation.sh
```

This script will:

1. Locate an existing BUSCO `full_table.tsv`, if present.
2. Otherwise run BUSCO in protein mode on the representative no-stop proteome.
3. Copy the BUSCO `full_table.tsv` and short summary into a stable `05_busco/` location.
4. Rerun `scripts/add_busco_duplication_context.py` using the primary genome-linked dark-candidate table and the BUSCO full table.
5. Rerun Figures 05–15 with BUSCO-backed prioritisation.

Previous BUSCO interpretation note:

```text
BUSCO v5.3.2
mode: proteins
lineage: metazoa_odb10
C:95.6% [S:52.9%, D:42.7%], F:2.6%, M:1.8%, n:954
```

The high duplicated BUSCO fraction is biologically/technically important and should be used as a conservative validation layer for dark-candidate prioritisation.

## Immediate next action

Run BUSCO-backed validation on Maxwell:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

git pull

sbatch scripts/run_busco_duplication_validation.sh
```

If BUSCO is not available as `busco`, load the relevant module or activate the correct conda environment, then rerun. If an existing BUSCO full table already exists, the script should reuse it. To force a known BUSCO full table, run:

```bash
BUSCO_FULL_TABLE=/path/to/full_table.tsv \
sbatch scripts/run_busco_duplication_validation.sh
```

After completion, inspect:

```bash
PROJECT_DIR=/uoa/scratch/users/r02hw22/project_dark_genes

cat "${PROJECT_DIR}/05_busco/equina_representative_longest_per_gene_metazoa_odb10/short_summary.txt"
cat "${PROJECT_DIR}/07_duplication_context/equina_duplication_context.summary.txt"
cat "${PROJECT_DIR}/08_figures/figure_tables/figure_14_priority_tier_counts.tsv"
```

The updated duplication summary should no longer say:

```text
BUSCO full table: NA
```

## Next phase after BUSCO-backed validation

After BUSCO-backed prioritisation has been rerun and committed, return to repeat/TE-overlap integration using:

```text
scripts/add_repeat_overlap_context.py
```

## Remaining biological validation steps

- [ ] Run BUSCO-backed duplication validation.
- [ ] Run repeat/TE-overlap context.
- [ ] Add methylation/DMR overlap context.
- [ ] Add RNA-seq expression evidence.
- [ ] Flag expressed versus unsupported dark candidates.
- [ ] Prioritise stress-responsive dark candidates.
- [ ] Produce final candidate shortlist with genome, duplication, BUSCO, repeat/TE, expression, methylation, and annotation-evidence fields.
