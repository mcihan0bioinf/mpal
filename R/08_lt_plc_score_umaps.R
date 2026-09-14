## 08_lt_plc_score_umaps.R
## Descriptive UMAPs colored by the continuous LT-PLC UCell score
## (up_all signature), one panel per atlas that has its own embedding
## (Granja/MPAL, Caron/B-ALL+T-ALL, van Galen/AML). This is the per-cell,
## visual companion to the patient-level statistics in
## 06_cross_atlas_comparison.R -- it shows *where* within each atlas the
## program is enriched, not whether the difference is significant.
##
## Input:  data/granja_seurat.rds, data/caron_seurat.rds,
##         data/vangalen_seurat.rds (all must carry $LT_PLC_up_UCell --
##         run 04_score_signature_all_atlases.R first)
## Output: results/lt_plc_score_umaps.png / .pdf

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(patchwork)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")

score_umap_panel <- function(path, title) {
  seu <- readRDS(path)
  stopifnot("LT_PLC_up_UCell" %in% colnames(seu@meta.data))
  emb <- Embeddings(seu, "umap")[, 1:2]
  df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2], score = seu$LT_PLC_up_UCell) %>%
    arrange(score)  # plot high scores last so they aren't buried under low ones

  ggplot(df, aes(UMAP_1, UMAP_2, color = score)) +
    geom_point(size = 0.35, alpha = 0.8) +
    scale_color_viridis_c(name = "UCell\nscore", option = "magma") +
    labs(title = title) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          panel.grid.minor = element_blank(),
          legend.key.height = unit(0.5, "cm"))
}

pA <- score_umap_panel(file.path(dir_data, "granja_seurat.rds"),   "Granja 2019 (MPAL + healthy)")
pB <- score_umap_panel(file.path(dir_data, "caron_seurat.rds"),    "Caron 2020 (B-ALL / T-ALL + healthy)")
pC <- score_umap_panel(file.path(dir_data, "vangalen_seurat.rds"), "van Galen 2019 (AML + healthy)")

p <- (pA | pB | pC) +
  plot_annotation(title = "LT-PLC signature score (UCell, up_all) across leukemia atlases",
                   subtitle = "each panel is that atlas's own UMAP embedding; color = per-cell UCell score",
                   theme = theme(plot.title = element_text(face = "bold", size = 13),
                                 plot.subtitle = element_text(size = 9.5, color = "grey30")))

ggsave(file.path(dir_out, "lt_plc_score_umaps.png"), p, width = 16, height = 5.5, dpi = 300)
ggsave(file.path(dir_out, "lt_plc_score_umaps.pdf"), p, width = 16, height = 5.5)
cat("Saved results/lt_plc_score_umaps.png / .pdf\n")
