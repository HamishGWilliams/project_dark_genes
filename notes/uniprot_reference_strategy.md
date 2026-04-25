# UniProt reference database strategy

Use **Swiss-Prot first** and treat it as the highest-confidence source for naming and functional interpretation.

Do **not** rely on Swiss-Prot alone for homology detection in this cnidarian dataset, because Swiss-Prot does not capture the full sequence diversity needed for non-model metazoan proteins.

Use **TrEMBL second**, but curate it rather than searching all TrEMBL entries indiscriminately.

Build the annotation workflow hierarchically:

1. Search against **UniProtKB/Swiss-Prot** first.
2. Search unresolved proteins against a **curated TrEMBL set**, preferably from **UniProt Reference Proteomes** and filtered to biologically relevant taxa such as **Metazoa** and, where possible, **Cnidaria-relevant lineages**.
3. Treat **Swiss-Prot hits as stronger functional evidence** than TrEMBL hits.
4. Treat **TrEMBL hits mainly as homology support**, unless InterPro, eggNOG, domain architecture, or other evidence corroborates the inferred function.

Apply this interpretation rule in the master annotation table:

- **Swiss-Prot hit** = preferred annotation source
- **TrEMBL hit only** = weaker annotation; describe conservatively
- **No sequence hit, but domain/orthology support** = not fully dark
- **No sequence hit and no domain/function support** = stronger dark-gene candidate

This project should therefore use a **hierarchical reference strategy**, not a single undifferentiated UniProt database.

# How to interpret the hits

""Swiss-Prot hit""

- Treat as your preferred annotation.
- If the hit is strong and coverage is good, this is your best sequence-based functional assignment.

**Cnidarian TrEMBL hit only**
- Treat as homology support, and sometimes family-level support.

Use conservative wording:

- “X-like protein”
- “putative homolog of X”
- “member of Y family”

Do not treat this as equal to reviewed Swiss-Prot evidence.

**Broad metazoan TrEMBL fallback hit only**
- Treat as weak-to-moderate sequence support unless corroborated.
- This can rescue a protein from being “sequence-dark,” but it should not automatically rescue it from being “function-dark.”

For a strong function claim, require corroboration from:
- InterPro domains
- eggNOG orthology
- domain architecture consistency
- ideally expression and genomic plausibility

**What this means for dark-gene logic**

classify results like this:

- Swiss-Prot hit ? annotated
- TrEMBL cnidarian/metazoan hit only ? sequence-supported, but lower confidence
- No sequence hit, but InterPro/eggNOG support ? not fully dark
- No Swiss-Prot, no curated TrEMBL, no InterPro, no eggNOG ? strong dark-gene candidate