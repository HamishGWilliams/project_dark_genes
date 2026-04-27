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

## Current issue to fix

The current DIAMOND/BLASTp workflow produced top-hit files by taking the first hit per query. This is not sufficient for final dark-gene classification because weak hits may be retained. The workflow must first filter hits by e-value and then select the top retained hit.

Current problematic pattern:

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

- [ ] Record the decision to use `e-value <= 1e-5` for DIAMOND/BLASTp annotation evidence.
- [ ] Note that this threshold is aligned with the Stephens et al. coral dark-gene workflow.
- [ ] Document that identity and coverage are retained for QC/reporting, not as hard filters in the main classification.
- [ ] Add a short methods note explaining ambiguous-description filtering.

### 2. DIAMOND script updates

- [ ] Add explicit `--evalue 1e-5` to DIAMOND searches.
- [ ] Consider using `--ultra-sensitive` for final searches.
- [ ] Consider using `--max-target-seqs 0` or another sufficiently exhaustive setting if retaining multiple hits for description filtering.
- [ ] Update DIAMOND output format to include `qlen` and `slen`:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
```

- [ ] Write raw DIAMOND outputs to `02_annotation/diamond/raw/`.
- [ ] Write e-value-filtered DIAMOND outputs to `02_annotation/diamond/filtered/`.
- [ ] Write filtered top-hit DIAMOND files after applying the e-value threshold.

### 3. BLASTp script updates

- [ ] Add explicit `-evalue 1e-5` to BLASTp searches.
- [ ] Update BLASTp output format to include `qlen` and `slen`:

```text
qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
```

- [ ] Write raw BLASTp outputs to `02_annotation/blastp/raw/`.
- [ ] Write e-value-filtered BLASTp outputs to `02_annotation/blastp/filtered/`.
- [ ] Write filtered top-hit BLASTp files after applying the e-value threshold.

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

- [ ] Regenerate or filter DIAMOND Swiss-Prot representative hits.
- [ ] Regenerate or filter DIAMOND Cnidaria-TrEMBL representative hits.
- [ ] Regenerate or filter BLASTp Swiss-Prot representative hits once the jobs finish.
- [ ] Regenerate or filter BLASTp Cnidaria-TrEMBL representative hits once the jobs finish.
- [ ] Confirm filtered files contain only hits with `e-value <= 1e-5`.

### 6. Update master annotation compiler

- [ ] Point the compiler to filtered DIAMOND/BLASTp top-hit files.
- [ ] Ensure unfiltered top-hit files are no longer used for classification.
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
- [ ] Confirm BLASTp fields match the filtered BLASTp source files when BLASTp finishes.
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

- [ ] Confirm full filtered DIAMOND files exist.
- [ ] Confirm full filtered BLASTp files exist once BLASTp jobs finish.
- [ ] Confirm full InterProScan output exists.
- [ ] Confirm full eggNOG output exists.
- [ ] Confirm full SignalP output exists.
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

## Current status note

The 100-protein master table structure and non-BLASTp source integration have already validated successfully. The remaining work is to revise homology filtering, regenerate filtered sequence-similarity files, then rebuild and QC the test and full master tables.
