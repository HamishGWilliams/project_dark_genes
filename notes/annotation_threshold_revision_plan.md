# Project Dark Genes annotation and candidate-tracking status

## Current phase

The core functional annotation, master-table construction, dark-candidate extraction, and first-pass figure generation have been completed on the `revise-homology-filtering` branch. The active issue is now genome-linking multiplicity: all 4,531 dark candidates were linked to the genome, but the current genome-linked table contains many-to-many matches and should be collapsed/QC-checked before interpreting scaffold, exon/span, BUSCO-duplication, cluster, and priority-tier figures.

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
- [X] Generate first-pass genome-linked output.
- [X] Generate first-pass figure set, Figures 02–15.
- [X] Generate figure manifest and compact figure tables.

First-pass genome-linking summary:

```text
Candidate rows read: 4531
Linked output rows: 574504
Unmatched candidates: 0
matched_multiple: 4531
```

Interpretation: every dark candidate found at least one genome match, but the mapping is currently many-to-many and inflates downstream linked-row counts.

## Active issue: genome-linking multiplicity

The current genome-linked table should not yet be used directly for final scaffold, exon/span, BUSCO-duplication, cluster, or priority-tier interpretation because one candidate can contribute many linked rows.

A new QC/collapse script has been added:

```text
scripts/qc_genome_linking_multiplicity.py
```

This script keeps the full multi-match genome-linked table but writes a deterministic one-row-per-candidate primary mapping table plus multiplicity diagnostics.

Expected outputs:

```text
06_genome_lookup/equina_dark_candidates.genome_linked.primary.tsv
06_genome_lookup/equina_dark_candidates.genome_linking_multiplicity.tsv
06_genome_lookup/equina_dark_candidates.genome_linking_multiplicity.summary.txt
```

## Immediate next action

Run genome-linking multiplicity QC on Maxwell:

```bash
cd /uoa/home/r02hw22/sharedscratch/project_dark_genes

git pull

python3 scripts/qc_genome_linking_multiplicity.py \
  --linked /uoa/scratch/users/r02hw22/project_dark_genes/06_genome_lookup/equina_dark_candidates.genome_linked.tsv \
  --outdir /uoa/scratch/users/r02hw22/project_dark_genes/06_genome_lookup \
  --candidate-id-column protein_id \
  --prefix equina_dark_candidates
```

Then inspect:

```bash
cat /uoa/scratch/users/r02hw22/project_dark_genes/06_genome_lookup/equina_dark_candidates.genome_linking_multiplicity.summary.txt
wc -l /uoa/scratch/users/r02hw22/project_dark_genes/06_genome_lookup/equina_dark_candidates.genome_linked.primary.tsv
```

The primary mapping table should contain 4,532 lines: one header plus 4,531 dark candidates.

## Next project step after QC

If the primary mapping table is clean:

1. Update figure-generation scripts to use `equina_dark_candidates.genome_linked.primary.tsv` for candidate-level summaries.
2. Regenerate Figures 05–15.
3. Re-check that priority-tier counts sum to 4,531 candidates, not hundreds of thousands of linked rows.
4. Commit only scripts, notes, summaries, figures, and compact figure tables.
5. Keep large full TSVs ignored/local.

## Remaining biological validation steps

- [ ] Add TE-overlap context for dark candidates.
- [ ] Add methylation/DMR overlap context.
- [ ] Add RNA-seq expression evidence.
- [ ] Flag expressed versus unsupported dark candidates.
- [ ] Prioritise stress-responsive dark candidates.
- [ ] Produce final candidate shortlist with genome, duplication, expression, methylation, and annotation-evidence fields.
