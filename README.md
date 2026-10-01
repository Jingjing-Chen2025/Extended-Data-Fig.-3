# Extended Data Figure 3 — specificity controls for stromal–epithelial inference

Standalone panel code. The source data are **bundled** in `data/`, so the figure runs on any
machine with R and four packages — no access to the analysis volume is needed.

```bash
Rscript "Extended Data Figure 3.R"               # writes output/Extended Data Figure 3.pdf
ED3_OUT=/some/path Rscript "Extended Data Figure 3.R"
```

Requirements: R >= 4.3 with `dplyr`, `ggplot2`, `patchwork`, `scales`.

## Panels

| Panel | Content | Input |
|---|---|---|
| a | Ligand breadth confounds interaction-strength scoring | `ligand_breadth.csv` |
| b | Ligand–receptor panel coupling against a matched random-panel null | `random_panel_null.csv` |

## Source data (`data/`)

| File | Rows | Contents |
|---|---|---|
| `ligand_breadth.csv` | 151 | % of 40 secreted/ECM ligands detected per cell, per dataset x population |
| `random_panel_null.csv` | 96 | Observed partial rho vs 200 matched random panel pairs, per panel combination |

## Why these two controls

Both are **negative** controls, and they are what licenses the main claim rather than supporting it
directly.

**a — repertoire breadth, not biology.** Interaction-strength methods such as CellChat sum
expression products over ligand–receptor pairs, so the population expressing the broadest ligand
repertoire ranks highest by construction. Fibroblast-lineage cells detect ~45% of a 40-gene panel
per cell against 7% for luminal epithelium and 2–5% for immune populations. A ranking of senders
that does not normalise for this cannot distinguish biology from breadth, which is why interaction
strength was not used to nominate a dominant sender.

**b — correlation is not an axis.** For every sender ligand panel by epithelial receptor panel, the
observed dataset-blocked partial Spearman correlation is compared with 200 random panel pairs
matched for gene count and expression bin and scored through the identical control-gene-corrected
path. **No coupling survives correction for the 144 tests in either species**; the strongest,
MDK/TNFSF12 to the midkine-receptor panel in mouse, reaches z = 2.18. So coordinated ligand and
receptor scores are not by themselves evidence of signalling, and the claims in the paper rest
instead on differential expression of matched ligand–receptor pairs and on the functional
experiments.

## Caveat carried on the figure

The breadth shown for the PI16+ progenitor-like stroma in panel a is currently that of the whole
fibroblast lineage. `ligand_breadth.csv` was built from cluster labels under which a relabelling
rule renamed every fibroblast cluster, so the table has no separate fibroblast row. Regenerating it
from the per-cell calls (pipeline v3.7, `cfg$aspc_cluster_rule = FALSE`) would separate the two.
This does not affect the panel's point, which is about breadth being a confounder.
