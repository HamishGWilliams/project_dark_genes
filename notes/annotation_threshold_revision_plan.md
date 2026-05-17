# Project Dark Genes annotation and candidate-tracking status

## Current phase

The `revise-homology-filtering` branch has been pushed successfully with large generated TSVs excluded from Git. The core functional annotation, master-table construction, dark-candidate extraction, first-pass genome linking, first-pass figure generation, and genome-linking multiplicity QC are now complete.

The active issue is now figure/prioritisation regeneration: the multiplicity QC produced a clean one-row-per-candidate primary genome-linking table, but the committed first-pass figure tables still reflect the pre-collapse many-to-many linked table. Figures 05–15 and priority summaries should therefore be regenerated from the primary mapping table before biological interpretation.

## Repository status

Checked after the successful push and multiplicity-QC update:

- Branch: `revise-homology-filtering`
- Remote branch visible on GitHub: yes
- Large rejected full genome-linked table removed from pushed diff: yes
- Full generated table remains local/generated and ignored: `03_dark_candidates/genomic_context/equina_dark_candidates.genome_linked.tsv`
- First-pass genome-linking summary is tracked.
- Multiplicity-QC summary is tracked.
- Figure outputs, figure manifest, compact figure tables, scripts, and notes are tracked.

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

### First-pass genome linking and figures

- [X] Run first-pass GFF3/genome lookup.
- [X] Generate first-pass genome-linked output locally.
- [X] Track first-pass genome-linking summary.
- [X] Generate first-pass figure set, Figures 02–15.
- [X] Generate figure manifest and compact figure tables.
- [X] Push cleaned branch with large generated TSVs excluded.

First-pass genome-linking summary:

```text
Candidate rows read: 4531
Linked output rows: 574504
Unmatched candidates: 0
matched_multiple: 4531
```

Interpretation: every dark candidate found at least one genome match, but the first-pass table is many-to-many and inflates downstream linked-row counts.

### Genome-linking multiplicity QC

- [X] Add `scripts/qc_genome_linking_multiplicity.py`.
- [X] Run genome-linking multiplicity QC.
- [X] Track the multiplicity-QC summary.
- [X] Produce local one-row-per-candidate primary mapping table.
- [X] Confirm primary table contains 4,531 dark-candidate rows.

Confirmed multiplicity-QC summary:

```text
Total linked rows read: 574504
Unique candidate IDs: 4531
Primary rows written: 4531
matched_primary_from_multiple_strict: 4531
rank_1: 4531
```

Expected local outputs:

```text
06_genome_lookup/equina_dark_candidates.genome_linked.primary.tsv
06_genome_lookup/equina_dark_candidates.genome_linking_multiplicity.tsv
06_genome_lookup/equina_dark_candidates.genome_linking_multiplicity.summary.txt
```

Only the summary should normally be committed; the full primary and multiplicity TSVs should remain ignored/local unless deliberately force-added.

## Active issue: regenerate downstream figures and prioritisation

The committed compact figure tables still reflect the first-pass many-to-many mapping. In particular:

```text
Figure 05: matched_multiple = 574504
Figure 14: high_priority + medium_priority + low_priority_or_manual_review = 574504
```

These should be treated as first-pass diagnostic outputs, not final candidate-level biological summaries.

## Immediate next action

Regenerate duplication/prioritisation and Figures 05–15 using the primary genome-linking table:

```text
/uoa/scratch/users/r02hw22/project_dark_genes/06_genome_lookup/equina_dark_candidates.genome_linked.primary.tsv
```

The regenerated candidate-level counts should sum to:

```text
4531
```

not:

```text
574504
```

Practical run plan:

1. Retarget any duplication/prioritisation command that currently uses `equina_dark_candidates.genome_linked.tsv` so that it uses `equina_dark_candidates.genome_linked.primary.tsv` instead.
2. Rebuild prioritisation outputs locally.
3. Rerun figure generation with `GENOME_LINKED_TSV` pointing to the primary table.
4. Confirm `08_figures/figure_tables/figure_05_genome_linking_status_counts.tsv` sums to 4,531.
5. Confirm `08_figures/figure_tables/figure_14_priority_tier_counts.tsv` sums to 4,531.
6. Commit only scripts, notes, summaries, figures, and compact figure tables.
7. Keep large full TSVs ignored/local.

## Remaining biological validation steps

- [ ] Add TE-overlap context for dark candidates.
- [ ] Add methylation/DMR overlap context.
- [ ] Add RNA-seq expression evidence.
- [ ] Flag expressed versus unsupported dark candidates.
- [ ] Prioritise stress-responsive dark candidates.
- [ ] Produce final candidate shortlist with genome, duplication, expression, methylation, and annotation-evidence fields.
