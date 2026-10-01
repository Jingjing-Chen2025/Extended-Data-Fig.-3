# Extended Data Figure 3 | Specificity controls for stromal-epithelial inference
#
# Standalone: reads the two bundled tables in data/ and writes output/Extended Data Figure 3.pdf.
# No access to the analysis volume is needed. Both tables are per-panel exports of upstream control
# experiments; this script performs no statistics, only layout.
#
#   a  ligand breadth confounds interaction-strength scoring
#   b  ligand-receptor panel coupling against a matched random-panel null
#
# These are the two negative controls that license the main claim. Panel a is why interaction
# strength was not used to nominate a dominant sender; panel b is why the pair-level evidence does
# not rest on a panel-level correlation that a random panel would reproduce.
#
# Requirements: R >= 4.3 with dplyr, ggplot2, patchwork, scales.
# Override the output directory with ED3_OUT=/some/path Rscript "Extended Data Figure 3.R"
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(patchwork); library(scales)})
DATA <- "data"; OUT <- Sys.getenv("ED3_OUT", unset = file.path(getwd(), "output"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.csv(file.path(DATA, f), check.names = FALSE)
set.seed(1)

BLUE <- "#3B6FA8"; RED <- "#C6553B"
COL_SPECIES <- c(human = BLUE, mouse = RED)
# "ASPC_adipose_progenitor" is the pipeline's internal name for cells passing the per-cell
# progenitor score. It is displayed by its markers, not as a lineage claim.
POP_LABEL <- c(ASPC_adipose_progenitor      = "PI16+ progenitor-like stroma",
               Fibroblast                   = "Fibroblast (other)",
               Smooth_muscle_myofibroblast  = "Smooth muscle / myofibroblast",
               Club_Hillock_epithelial      = "Club / hillock epithelial",
               Contaminant_adrenal_cortex   = "Adrenal cortex (contaminant)")
relabel_pop <- function(p) ifelse(p %in% names(POP_LABEL), unname(POP_LABEL[p]), gsub("_", " ", p))
th <- theme_bw(base_size = 6) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95", colour = NA),
        plot.title = element_text(face = "bold", size = 7), plot.subtitle = element_text(size = 5.2, colour = "grey30"),
        legend.key.size = unit(3, "mm"), legend.text = element_text(size = 5.5), plot.title.position = "plot")

# ---------------- a: ligand breadth
B <- rd("ligand_breadth.csv")
ord <- B %>% group_by(population) %>% summarise(m = median(ligand_breadth), .groups = "drop") %>% arrange(m)
pa <- ggplot(B %>% mutate(population = factor(relabel_pop(population), levels = relabel_pop(ord$population))),
             aes(ligand_breadth, population, colour = species)) +
  geom_jitter(height = 0.15, size = 0.6, alpha = 0.7) +
  stat_summary(fun = median, geom = "point", shape = 124, size = 2, colour = "black") +
  scale_colour_manual(values = COL_SPECIES, name = NULL) +
  th + theme(legend.position = "bottom") +
  labs(title = "a  Ligand breadth confounds interaction-strength scoring",
       subtitle = paste("% of 40 secreted/ECM ligands detected per cell. Methods that sum over ligand-receptor",
                        "pairs rank the broadest expressor first irrespective of biology. The breadth shown for the",
                        "PI16+ progenitor-like stroma is that of the fibroblast lineage until the source table is regenerated.", sep = "\n"),
       x = "% of surveyed ligands detected per cell", y = NULL)

# ---------------- b: random-panel null
r1 <- rd("random_panel_null.csv") %>%
  mutate(lig = factor(ligand_panel,
                      levels = c("Inflammatory", "Matrix", "MSC", "Neurotrophic_full", "Neurotrophic_core", "MDK_TNFSF12"),
                      labels = c("Inflammatory", "Matrix", "MSC", "Neurotrophic full", "Neurotrophic core", "MDK/TNFSF12")),
         rec = factor(sub("_receptor$", "", receptor_panel),
                      levels = c("MDK_PTN", "TNFSF12", "NRG", "SLIT", "SEMA", "IGF1", "CXCL12", "Neurotrophin"),
                      labels = c("MDK/PTN", "TNFSF12", "NRG", "SLIT", "SEMA", "IGF1", "CXCL12", "NTRK/NGFR")),
         sig = ifelse(fdr_empirical < 0.1, "*", ""))
pb <- ggplot(r1, aes(rec, lig, fill = z_vs_null)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%s", z_vs_null, sig)), size = 1.7) +
  facet_wrap(~ species) +
  scale_fill_gradient2(low = BLUE, mid = "white", high = RED, midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish, name = "z vs null") +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
             legend.position = "right", legend.key.width = grid::unit(0.18, "cm")) +
  labs(title = "b  Ligand-receptor coupling against a matched random-panel null",
       subtitle = paste("z of the observed partial rho relative to 200 random panel pairs matched for gene count and expression",
                        "bin. * FDR < 0.1. No coupling survives correction for 144 tests in either species.", sep = "\n"),
       x = "epithelial receptor panel", y = "sender ligand panel")

fig <- pa / pb + plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Extended Data Figure 3 | Specificity controls for stromal-epithelial inference",
    caption = paste(
      "a, Fraction of a 40-gene secreted and extracellular-matrix ligand panel detected per cell, by population, across all 16 single-cell datasets. Interaction-strength methods sum expression products over",
      "ligand-receptor pairs, so the population with the broadest ligand repertoire ranks highest by construction; fibroblast-lineage cells detect approximately 45% of the panel per cell compared with 7% for luminal",
      "epithelium. This is why ligand-receptor interaction strength was not used to nominate a dominant sender. b, For each sender ligand panel by epithelial receptor panel, the z of the observed dataset-blocked",
      "partial Spearman correlation relative to 200 random panel pairs matched for gene count and expression bin, scored through the identical control-gene-corrected path. Asterisks would denote FDR < 0.1; no coupling",
      "survives correction for the 144 tests performed in either species, so correlation between stromal ligand and epithelial receptor scores is not by itself evidence of a specific axis.",
      sep = "\n"),
    theme = theme(plot.title = element_text(face = "bold", size = 8.5, family = "Helvetica", margin = margin(b = 6)),
                  plot.caption = element_text(size = 4.6, colour = "grey40", hjust = 0, lineheight = 1.35)))
ggsave(file.path(OUT, "Extended Data Figure 3.pdf"), fig, width = 6.4, height = 6.6, device = cairo_pdf)
cat("written", file.path(OUT, "Extended Data Figure 3.pdf"), "\n")
