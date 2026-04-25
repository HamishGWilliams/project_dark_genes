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
