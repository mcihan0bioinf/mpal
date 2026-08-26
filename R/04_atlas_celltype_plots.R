## 04_atlas_celltype_plots.R
## UMAP of each atlas colored and labeled by author-provided cell type.
##
## Input:  data/granja_seurat.rds
##         data/bonemarrowmap_seurat.rds
##         data/caron_seurat.rds
##         data/vangalen_seurat.rds
## Output: results/atlas_celltypes_<atlas>.png (one per atlas)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

plot_atlas <- function(rds_path, title, out_png) {
  if (!file.exists(rds_path)) {
    cat(sprintf("Skipping %s (not found: %s)\n", title, rds_path))
    return(invisible(NULL))
  }
  seu <- readRDS(rds_path)
  emb <- Embeddings(seu, "umap")[, 1:2]
  df <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                    celltype = as.character(seu$celltype))
  df <- df %>% filter(!is.na(celltype), celltype != "")

  # label at the median position of each cell type cluster
  label_df <- df %>% group_by(celltype) %>%
    summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), n = n(), .groups = "drop")

  n_types <- length(unique(df$celltype))
  pal <- if (n_types <= 12) scales::hue_pal()(n_types) else
    colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_types)

  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = celltype)) +
    geom_point(size = 0.3, alpha = 0.5) +
    ggrepel::geom_text_repel(data = label_df, aes(label = celltype), inherit.aes = TRUE,
                              color = "black", size = 3.2, fontface = "bold",
                              bg.color = "white", bg.r = 0.15, max.overlaps = 20) +
    scale_color_manual(values = pal) +
    guides(color = "none") +
    labs(title = title, subtitle = sprintf("%d cells, %d cell types", nrow(df), n_types),
         x = "UMAP 1", y = "UMAP 2") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 9.5, color = "grey30"),
          panel.grid.minor = element_blank())

  ggsave(out_png, p, width = 8, height = 7, dpi = 300)
  cat(sprintf("Saved %s (%d cells, %d cell types)\n", out_png, nrow(df), n_types))
  rm(seu); gc(verbose = FALSE)
}

plot_atlas(file.path(dir_data, "granja_seurat.rds"),
           "Granja et al. 2019 -- MPAL + healthy hematopoiesis",
           file.path(dir_out, "atlas_celltypes_granja.png"))

plot_atlas(file.path(dir_data, "bonemarrowmap_seurat.rds"),
           "BoneMarrowMap -- healthy bone marrow reference",
           file.path(dir_out, "atlas_celltypes_bonemarrowmap.png"))

plot_atlas(file.path(dir_data, "caron_seurat.rds"),
           "Caron et al. 2020 -- B-ALL / T-ALL / healthy (PBMMC)",
           file.path(dir_out, "atlas_celltypes_caron.png"))

plot_atlas(file.path(dir_data, "vangalen_seurat.rds"),
           "van Galen et al. 2019 -- AML + healthy BM",
           file.path(dir_out, "atlas_celltypes_vangalen.png"))

cat("\nDone.\n")
