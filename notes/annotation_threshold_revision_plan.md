# Annotation threshold revision plan

## Purpose

This note records the planned revision to the Project Dark Genes functional annotation workflow. The immediate goal is to move from unfiltered top-hit sequence similarity evidence to a defensible, paper-aligned filtering approach before rebuilding the master annotation tables.

The key methodological change is that DIAMOND and BLASTp evidence should be filtered before being used in the master annotation table.

## Basis for threshold choice

The coral dark-gene paper by Stephens et al. used DIAMOND BLASTP against NCBI nr with an e-value threshold of `1e-5` for protein functional annotation. Proteins were considered functionally annotated if at least one qualifying hit had an informative description. Proteins were considered unknown-function if all qualifying hits had ambiguous descriptions such as `uncharacterized protein`, `hypothetical protein`, `predicted protein`, `expressed protein`, or `unnamed protein product`.

For this project, the primary sequence-similarity threshold will therefore be:

```text
e-value <= 1e-5
```

Percent identity, alignment length, bitscore, query coverage, and subject coverage should be retained where possible as QC/reporting fields, but they should not be used as hard filters for the main annotation classes unless a separate high-confidence subcategory is created.

## Current status

Updated from the HPC file audit supplied on 2026-05-17.

The major annotation evidence files are now present for both the 100-protein test set and the full representative proteome:

- Filtered DIAMOND Swiss-Prot and Cnidaria-TrEMBL outputs are present in `02_annotation/diamond/filtered/`.
- Filtered BLASTp Swiss-Prot and Cnidaria-TrEMBL outputs are present in `02_annotation/blastp/filtered/`.
- Full and test InterProScan TSV/GFF3/XML outputs are present in `02_annotation/interproscan/raw/`.
- Full and test eggNOG-mapper outputs are present in `02_annotation/eggnog/raw/`.
- Full and test SignalP summary and positive-prediction outputs are present in `02_annotation/signalp/summary/`.
- The master annotation script has been built.

The active phase is now master-table rebuild, QC, and dark-candidate extraction.

## Filtering issue addressed

The previous DIAMOND/BLASTp workflow produced top-hit files by taking the first hit per query. This was not sufficient for final dark-gene classification because weak hits could be retained. The revised workflow filters hits by e-value before selecting the retained top hit.

Previous problematic pattern:

```bash
awk -F '\t' '!seen[$1]++' raw_hits.tsv > top_hits.tsv
```

Revised pattern:

```bash
awk -F '\t' 'BEGIN {OFS="\t"} $11 <= 1e-5 {print}' raw_hits.tsv > filtered_hits.tsv
awk -F '\t' '!seen[$1]++' filtered_hits.tsv > filtered_top_hits.tsv
```

This assumes the standard outfmt field order:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
```

For future reruns, include query and subject lengths:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
```

## Revised annotation hierarchy

The master table should assign one primary annotation class per representative protein using this hierarchy:

1. `annotated_swissprot_supported`
2. `sequence_supported_trembl_cnidaria`
3. `domain_supported_interpro`
4. `orthology_supported_eggnog`
5. `function_dark_but_signalp_secretory_candidate`
6. `function_dark_no_current_annotation`

Only threshold-passing sequence hits should be allowed to assign categories 1 or 2.

## Checklist of tasks

### 1. Documentation

- [X] Record the decision to use `e-value <= 1e-5` for DIAMOND/BLASTp annotation evidence.
- [X] Note that this threshold is aligned with the Stephens et al. coral dark-gene workflow.
- [X] Document that identity and coverage are retained for QC/reporting, not as hard filters in the main classification.
- [ ] Add a short methods note explaining ambiguous-description filtering.

### 2. DIAMOND script/output updates

- [X] Add explicit `--evalue 1e-5` to DIAMOND searches or apply equivalent post-search e-value filtering.
- [ ] Consider using `--ultra-sensitive` for final reruns, if reruns become necessary.
- [ ] Consider using `--max-target-seqs 0` or another sufficiently exhaustive setting if retaining multiple hits for description filtering.
- [ ] Confirm DIAMOND output format includes `qlen` and `slen` where needed for coverage QC:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
```

- [ ] Confirm raw DIAMOND outputs are retained in `02_annotation/diamond/raw/`.
- [X] Write e-value-filtered DIAMOND outputs to `02_annotation/diamond/filtered/`.
- [X] Write filtered top-hit DIAMOND files after applying the e-value threshold.
- [X] Confirm test100 filtered DIAMOND Swiss-Prot and Cnidaria-TrEMBL outputs exist.
- [X] Confirm full filtered DIAMOND Swiss-Prot and Cnidaria-TrEMBL outputs exist.

### 3. BLASTp script/output updates

- [X] Add explicit `-evalue 1e-5` to BLASTp searches or apply equivalent post-search e-value filtering.
- [ ] Confirm BLASTp output format includes `qlen` and `slen` where needed for coverage QC:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
```

- [ ] Confirm raw BLASTp outputs are retained in `02_annotation/blastp/raw/`.
- [X] Write e-value-filtered BLASTp outputs to `02_annotation/blastp/filtered/`.
- [X] Write filtered best-hit BLASTp files after applying the e-value threshold.
- [X] Confirm test100 filtered BLASTp Swiss-Prot and Cnidaria-TrEMBL outputs exist.
- [X] Confirm full filtered BLASTp Swiss-Prot and Cnidaria-TrEMBL outputs exist.

### 4. Ambiguous-description filtering

Where subject descriptions are available, especially for TrEMBL-style hits:

- [ ] Flag hits with ambiguous descriptions.
- [ ] Treat the following as ambiguous:

```text
uncharacterized protein
hypothetical protein
predicted protein
expressed protein
unnamed protein product
```

- [ ] Decide whether ambiguous Cnidaria-TrEMBL hits should be excluded from `sequence_supported_trembl_cnidaria` or retained in a separate field such as `sequence_match_ambiguous_description`.
- [ ] Prefer excluding ambiguous TrEMBL descriptions from the main annotation class where descriptions are available.

### 5. Rebuild 100-protein test homology files

- [X] Regenerate or filter DIAMOND Swiss-Prot representative hits.
- [X] Regenerate or filter DIAMOND Cnidaria-TrEMBL representative hits.
- [X] Regenerate or filter BLASTp Swiss-Prot representative hits.
- [X] Regenerate or filter BLASTp Cnidaria-TrEMBL representative hits.
- [ ] Confirm filtered files contain only hits with `e-value <= 1e-5` by direct column check.

### 6. Master annotation compiler

- [X] Build the master annotation script.
- [X] Point the compiler to filtered DIAMOND/BLASTp top-hit files.
- [ ] Confirm unfiltered top-hit files are no longer used for classification.
- [ ] Add or retain columns for:

```text
pident
alignment_length
evalue
bitscore
qlen
slen
query_coverage
subject_coverage
```

- [ ] Add pass/fail fields where useful:

```text
diamond_swissprot_evalue_pass
diamond_trembl_cnidaria_evalue_pass
blastp_swissprot_evalue_pass
blastp_trembl_cnidaria_evalue_pass
```

- [ ] Keep the six-class hierarchy in a fixed order in the summary output.
- [ ] Keep counts, proportions, and percentages for all six classes, including categories with zero proteins.

### 7. Rebuild and QC the 100-protein master table

- [ ] Rebuild `02_annotation/master/test100/equina_representative_test100.master_annotation.tsv`.
- [ ] Rebuild `02_annotation/master/test100/equina_representative_test100.master_annotation.summary.txt`.
- [ ] Rerun `scripts/qc_master_annotation_test100.sh`.
- [ ] Confirm FASTA IDs match the master table.
- [ ] Confirm protein lengths match the FASTA.
- [ ] Confirm representative lookup fields match.
- [ ] Confirm DIAMOND fields match the filtered DIAMOND source files.
- [ ] Confirm BLASTp fields match the filtered BLASTp source files.
- [ ] Confirm InterProScan fields match source.
- [ ] Confirm eggNOG fields match source.
- [ ] Confirm SignalP fields match source.
- [ ] Confirm annotation classes match the hierarchy.
- [ ] Confirm summary counts match the master table.

### 8. Rebuild full representative master table

- [ ] Confirm the full representative FASTA exists:

```text
02_annotation/input/equina_representative_longest_per_gene.no_stop.fa
```

- [X] Confirm full filtered DIAMOND files exist.
- [X] Confirm full filtered BLASTp files exist.
- [X] Confirm full InterProScan output exists.
- [X] Confirm full eggNOG output exists.
- [X] Confirm full SignalP output exists.
- [ ] Run `scripts/build_master_annotation_full.py` using filtered homology files.
- [ ] Create a full QC script equivalent to the test100 QC script.
- [ ] Confirm the full master table has one row per representative protein.

### 9. Generate dark-candidate outputs

After the final full master table is rebuilt and QC-passed:

- [ ] Extract `function_dark_no_current_annotation` proteins.
- [ ] Extract `function_dark_but_signalp_secretory_candidate` proteins.
- [ ] Create dark-candidate TSV files.
- [ ] Create dark-candidate FASTA files.
- [ ] Summarise candidate counts and proportions.
- [ ] Record best subthreshold sequence hits, if any, for later manual review.

### 10. Future validation

- [ ] Map candidates back to GFF3/genome coordinates.
- [ ] Summarise exon count, CDS span, transcript span, and scaffold location.
- [ ] Add RNA-seq expression evidence.
- [ ] Flag expressed versus unsupported dark candidates.
- [ ] Prioritise stress-responsive dark candidates.

## Expected deliverables after revision

```text
02_annotation/diamond/raw/
02_annotation/diamond/filtered/
02_annotation/blastp/raw/
02_annotation/blastp/filtered/
02_annotation/master/test100/equina_representative_test100.master_annotation.tsv
02_annotation/master/test100/equina_representative_test100.master_annotation.summary.txt
02_annotation/master/test100/qc/
02_annotation/master/full/equina_representative_full.master_annotation.tsv
02_annotation/master/full/equina_representative_full.master_annotation.summary.txt
02_annotation/master/full/qc/
03_dark_candidates/
```

## Immediate next action

Run the master annotation compiler on the 100-protein test set first, QC the output, then run the same workflow on the full representative proteome. Only after the full master table passes QC should dark-candidate TSV and FASTA files be generated.
