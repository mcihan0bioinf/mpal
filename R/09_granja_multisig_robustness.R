## 09_granja_multisig_robustness.R
## Robustness check: does the Granja MPAL patient x cell-state pattern (and
## its UMAP footprint) depend on which LT-PLC gene-set variant is used to
## score it? Scores all human ortholog signature variants already built by
## 01_ortholog_mapping.R (maximal mapping) on the Granja atlas with UCell,
## then produces one multi-panel figure per view (heatmap grid, UMAP grid)
## so the seven variants can be compared side by side.
##
## Variants scored:
##   up_all       - all significant upregulated genes (padj<0.05), no LFC cutoff
##   up_lfc1      - upregulated, LFC>=1
##   up_lfc2      - upregulated, LFC>=2
##   up_top50     - top 50 upregulated by LFC
##   up_top100    - top 100 upregulated by LFC
##   up_top500    - top 500 upregulated by LFC
##   up_and_down  - up_all as positive markers + down_lfc1 as negative
##                  markers (UCell "-" suffix convention, w_neg=1)
##
## Input:  data/granja_seurat.rds (celltype, celltype_fine, sample, umap)
##         data/human_signature_sets.rds (maximal mapping)
## Output: results/granja_patient_heatmap_multisig.png / .pdf
##         results/granja_score_umap_multisig.png / .pdf
##         results/granja_score_violin_multisig.png / .pdf
##         data/granja_seurat.rds (updated in place: adds <variant>_UCell columns)

suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(ggplot2); library(dplyr); library(tidyr)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")

human_sets <- readRDS(file.path(dir_data, "human_signature_sets.rds"))

variant_order <- c("up_all", "up_lfc1", "up_lfc2", "up_top50", "up_top100", "up_top500", "up_and_down")
variant_labels <- c(
  up_all      = "Up: all (padj<0.05)",
  up_lfc1     = "Up: LFC>=1",
  up_lfc2     = "Up: LFC>=2",
  up_top50    = "Up: top 50",
  up_top100   = "Up: top 100",
  up_top500   = "Up: top 500",
  up_and_down = "Up (all) + Down (LFC<=-1)"
)

variants <- list(
  up_all      = human_sets$up_all,
  up_lfc1     = human_sets$up_lfc1,
  up_lfc2     = human_sets$up_lfc2,
  up_top50    = human_sets$up_top50,
  up_top100   = human_sets$up_top100,
  up_top500   = human_sets$up_top500,
  up_and_down = c(human_sets$up_all, paste0(human_sets$down_lfc1, "-"))
)
for (nm in names(variants)) cat(sprintf("  %-12s n_features = %d\n", nm, length(variants[[nm]])))

## ---- score Granja atlas for all variants at once, in memory-bounded batches ----
## (skipped if a previous run already saved these columns to data/granja_seurat.rds)
seu <- readRDS(file.path(dir_data, "granja_seurat.rds"))
DefaultAssay(seu) <- "RNA"

variant_cols <- paste0(names(variants), "_UCell")
if (all(variant_cols %in% colnames(seu@meta.data))) {
  cat("All variant scores already present in data/granja_seurat.rds -- skipping rescoring.\n\n")
  scores <- as.matrix(seu@meta.data[, variant_cols])
  colnames(scores) <- sub("_UCell$", "", colnames(scores))
} else {
  batch_size <- 10000
  cells <- colnames(seu)
  batches <- split(cells, ceiling(seq_along(cells) / batch_size))
  scores_list <- vector("list", length(batches))
  for (i in seq_along(batches)) {
    sub <- subset(seu, cells = batches[[i]])
    scores_list[[i]] <- ScoreSignatures_UCell(GetAssayData(sub, layer = "counts"),
                                               features = variants, maxRank = 4500)
    rm(sub); gc(verbose = FALSE)
    cat(sprintf("  batch %d/%d scored\n", i, length(batches)))
  }
  scores <- do.call(rbind, scores_list)[cells, , drop = FALSE]
  colnames(scores) <- sub("_UCell$", "", colnames(scores))
  for (nm in colnames(scores)) seu[[paste0(nm, "_UCell")]] <- scores[, nm]
  saveRDS(seu, file.path(dir_data, "granja_seurat.rds"))
  cat("Scored and saved all variants to data/granja_seurat.rds\n\n")
}

mpal_labels <- c("Progenitor_Like","Lymphoid_Like","Myeloid_Like","Healthy_Like","Erythroid_Like","TNK_Like")
meta <- seu@meta.data %>%
  mutate(cell = rownames(.),
         mpal_class = as.character(celltype),
         patient = sub("_T[0-9]+$", "", sample)) %>%
  filter(mpal_class %in% mpal_labels, grepl("^MPAL", patient))

## ================= Panel 1: patient x cell-state heatmap grid =================
cat("=== Patient x cell-state heatmap, per signature variant ===\n")

df_long <- meta %>%
  select(patient, mpal_class, all_of(paste0(variant_order, "_UCell"))) %>%
  pivot_longer(cols = ends_with("_UCell"), names_to = "variant", values_to = "score") %>%
  mutate(variant = sub("_UCell$", "", variant)) %>%
  group_by(variant, patient, mpal_class) %>%
  summarise(median_score = median(score), n = n(), .groups = "drop") %>%
  filter(n >= 20)

## fix patient/cell-state ordering from the up_all variant so all panels are comparable
ref <- df_long %>% filter(variant == "up_all")
state_order <- ref %>% group_by(mpal_class) %>% summarise(m = median(median_score)) %>% arrange(m) %>% pull(mpal_class)
patient_order <- ref %>% group_by(patient) %>% summarise(m = median(median_score)) %>% arrange(m) %>% pull(patient)

df_long <- df_long %>% mutate(
  mpal_class = factor(mpal_class, levels = state_order),
  patient = factor(patient, levels = patient_order),
  variant = factor(variant, levels = variant_order, labels = variant_labels[variant_order])
)

p_heat <- ggplot(df_long, aes(x = mpal_class, y = patient, fill = median_score)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", median_score)), size = 2.2) +
  facet_wrap(~variant, ncol = 4) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", name = "Median\nUCell score") +
  labs(title = "Granja et al. 2019: LT-PLC patient x cell-state reproducibility across signature variants",
       subtitle = "Median UCell score per patient x MPAL cell-state (min 20 cells); rows/columns ordered by the up_all variant",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        panel.grid = element_blank())

ggsave(file.path(dir_out, "granja_patient_heatmap_multisig.png"), p_heat, width = 15, height = 8, dpi = 300)
ggsave(file.path(dir_out, "granja_patient_heatmap_multisig.pdf"), p_heat, width = 15, height = 8)
cat("Saved results/granja_patient_heatmap_multisig.png / .pdf\n\n")

## ================= Panel 2: UMAP grid, one panel per signature variant =================
cat("=== UMAP colored by score, per signature variant ===\n")

emb <- Embeddings(seu, "umap")[, 1:2]
df_umap <- data.frame(UMAP_1 = emb[, 1], UMAP_2 = emb[, 2],
                       as.data.frame(scores[colnames(seu), , drop = FALSE])) %>%
  pivot_longer(cols = all_of(variant_order), names_to = "variant", values_to = "score") %>%
  mutate(variant = factor(variant, levels = variant_order, labels = variant_labels[variant_order])) %>%
  group_by(variant) %>%
  ## min-max normalize within each variant: signature size drives the raw UCell
  ## score's magnitude (up_all vs. top50), so a shared color scale would bury
  ## the spatial pattern of the smaller gene sets -- this puts every panel on
  ## the same 0-1 visual scale so patterns are comparable across variants.
  mutate(score_norm = (score - min(score)) / (max(score) - min(score))) %>%
  ungroup() %>%
  arrange(variant, score_norm)  # plot high scores last within each facet

p_umap <- ggplot(df_umap, aes(UMAP_1, UMAP_2, color = score_norm)) +
  geom_point(size = 0.25, alpha = 0.8) +
  facet_wrap(~variant, ncol = 4) +
  scale_color_viridis_c(name = "UCell score\n(min-max\nper panel)", option = "magma") +
  labs(title = "Granja et al. 2019: LT-PLC UMAP score across signature variants",
       subtitle = "each panel = same UMAP embedding; color = that variant's per-cell UCell score, min-max normalized within panel (raw scale differs by signature size)") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "granja_score_umap_multisig.png"), p_umap, width = 15, height = 8, dpi = 300)
ggsave(file.path(dir_out, "granja_score_umap_multisig.pdf"), p_umap, width = 15, height = 8)
cat("Saved results/granja_score_umap_multisig.png / .pdf\n")

## ================= Panel 3: violin by MPAL cell-state, per signature variant =================
cat("\n=== Violin by MPAL cell-state, per signature variant ===\n")

## state_order (from the heatmap panel) excludes TNK_Like -- no single patient
## has >=20 TNK_Like cells, so it never entered the patient x cell-state pool;
## drop it here too rather than let it fall through as a misleading "NA" bar.
df_violin <- meta %>%
  filter(mpal_class %in% state_order) %>%
  select(mpal_class, all_of(paste0(variant_order, "_UCell"))) %>%
  pivot_longer(cols = ends_with("_UCell"), names_to = "variant", values_to = "score") %>%
  mutate(variant = sub("_UCell$", "", variant),
         mpal_class = factor(mpal_class, levels = state_order),  # same order as the heatmap
         variant = factor(variant, levels = variant_order, labels = variant_labels[variant_order]))

pal_state <- scales::hue_pal()(length(state_order)); names(pal_state) <- state_order

p_violin <- ggplot(df_violin, aes(x = mpal_class, y = score, fill = mpal_class)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
  facet_wrap(~variant, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = pal_state, guide = "none") +
  labs(title = "Granja et al. 2019: LT-PLC score by MPAL cell-state, across signature variants",
       subtitle = "cells pooled across patients; TNK_Like excluded (no patient has >=20 cells, same as the heatmap panel); y-axis free per panel; x order = up_all median",
       x = NULL, y = "UCell score") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        strip.text = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "granja_score_violin_multisig.png"), p_violin, width = 15, height = 8, dpi = 300)
ggsave(file.path(dir_out, "granja_score_violin_multisig.pdf"), p_violin, width = 15, height = 8)
cat("Saved results/granja_score_violin_multisig.png / .pdf\n")
