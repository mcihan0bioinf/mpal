## 05_granja_robustness.R
## Three robustness/specificity follow-ups on the Granja et al. 2019 MPAL
## atlas, using the LT-PLC "up_all" UCell score computed in
## 04_score_signature_all_atlases.R (data/granja_seurat.rds$LT_PLC_up_UCell):
##
## Analysis 1 -- per-patient reproducibility: median LT-PLC score per MPAL
##   patient x MPAL cell-state, to check the enrichment isn't driven by one
##   patient's cells dominating the pool.
##
## Analysis 2 -- lineage-matched comparison: for each MPAL disease-classified
##   compartment, compare its LT-PLC score to the healthy reference cell
##   states of the *same* broad lineage, to test whether the malignant state
##   is elevated relative to its normal counterpart (not just "myeloid genes
##   score higher because myeloid cells express more genes generally").
##
## Analysis 3 -- driver genes: for the two compartments that DO show
##   elevation over their matched healthy counterpart (Erythroid_Like,
##   Lymphoid_Like), identify which signature genes drive the difference.
##
## Input:  data/granja_seurat.rds (celltype, celltype_fine, sample, LT_PLC_up_UCell)
##         data/human_signature_sets.rds (maximal mapping; up_all set)
## Output: results/granja_patient_heatmap.png / .pdf
##         results/granja_lineage_matched_comparison.png / .pdf
##         results/granja_lymphoid_driver_genes.csv
##         results/granja_erythroid_driver_genes.csv

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "results")

seu <- readRDS(file.path(dir_data, "granja_seurat.rds"))
stopifnot("LT_PLC_up_UCell" %in% colnames(seu@meta.data))

mpal_labels <- c("Progenitor_Like","Lymphoid_Like","Myeloid_Like","Healthy_Like","Erythroid_Like","TNK_Like")

meta <- seu@meta.data %>%
  mutate(cell = rownames(.),
         mpal_class = as.character(celltype),
         fine_class = as.character(celltype_fine),
         patient = sub("_T[0-9]+$", "", sample))  # e.g. MPAL4_T1/MPAL4_T2 -> MPAL4

## ================= Analysis 1: per-patient reproducibility =================
cat("=== Analysis 1: per-patient reproducibility ===\n")

dfA2 <- meta %>%
  filter(mpal_class %in% mpal_labels, grepl("^MPAL", patient)) %>%
  group_by(patient, mpal_class) %>%
  summarise(median_score = median(LT_PLC_up_UCell), n = n(), .groups = "drop")

cat(sprintf("Patients included: %s\n", paste(sort(unique(dfA2$patient)), collapse = ", ")))
print(dfA2 %>% arrange(mpal_class, desc(median_score)))

## keep only patient x cell-state combos with >=20 cells (avoid noisy small cells)
dfA2_f <- dfA2 %>% filter(n >= 20)

state_order <- dfA2_f %>% group_by(mpal_class) %>% summarise(m = median(median_score)) %>%
  arrange(m) %>% pull(mpal_class)
patient_order <- dfA2_f %>% group_by(patient) %>% summarise(m = median(median_score)) %>%
  arrange(m) %>% pull(patient)

dfA2_f <- dfA2_f %>% mutate(mpal_class = factor(mpal_class, levels = state_order),
                             patient = factor(patient, levels = patient_order))

pA2 <- ggplot(dfA2_f, aes(x = mpal_class, y = patient, fill = median_score)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f\n(n=%d)", median_score, n)), size = 2.8, lineheight = 0.85) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", name = "Median\nUCell score") +
  labs(title = "Granja et al. 2019: LT-PLC score reproducibility across MPAL patients",
       subtitle = "Median UCell score (up_all signature) per patient x MPAL cell-state; cells pooled across timepoints; min 20 cells/cell",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        panel.grid = element_blank())

ggsave(file.path(dir_out, "granja_patient_heatmap.png"), pA2, width = 9, height = 6, dpi = 300)
ggsave(file.path(dir_out, "granja_patient_heatmap.pdf"), pA2, width = 9, height = 6)
cat("Saved results/granja_patient_heatmap.png / .pdf\n\n")

## ================= Analysis 2: lineage-matched comparison =================
cat("=== Analysis 2: lineage-matched MPAL vs healthy comparison ===\n")

## define matched healthy counterparts for each MPAL compartment, by broad lineage
lineage_map <- list(
  Progenitor_Like = c("HSC", "CMP.LMPP"),
  Myeloid_Like    = c("GMP", "GMP.Neut", "CD14.Mono.1", "CD14.Mono.2", "CD16.Mono"),
  Erythroid_Like  = c("Early.Eryth", "Late.Eryth", "Early.Baso"),
  Lymphoid_Like   = c("CLP.1", "CLP.2", "Pre.B", "B"),
  TNK_Like        = c("CD4.N1", "CD4.N2", "CD4.M", "CD8.N", "CD8.CM", "CD8.EM", "NK"),
  Healthy_Like    = c("CD4.N1", "CD4.N2", "CD4.M", "CD8.N", "CD8.CM", "CD8.EM", "NK")
)

comparison_list <- list()
for (mp in names(lineage_map)) {
  mpal_scores <- meta %>% filter(mpal_class == mp) %>% pull(LT_PLC_up_UCell)
  healthy_scores <- meta %>% filter(fine_class %in% lineage_map[[mp]]) %>% pull(LT_PLC_up_UCell)
  if (length(mpal_scores) < 5 || length(healthy_scores) < 5) next
  wt <- wilcox.test(mpal_scores, healthy_scores)
  comparison_list[[mp]] <- data.frame(
    compartment = mp,
    group = c(rep("MPAL (malignant)", length(mpal_scores)), rep("Healthy (matched lineage)", length(healthy_scores))),
    score = c(mpal_scores, healthy_scores)
  )
  cat(sprintf("%-16s: MPAL median=%.4f (n=%d) vs Healthy median=%.4f (n=%d), matched to {%s}, Wilcoxon p=%.2e\n",
              mp, median(mpal_scores), length(mpal_scores), median(healthy_scores), length(healthy_scores),
              paste(lineage_map[[mp]], collapse=","), wt$p.value))
}
dfA3 <- bind_rows(comparison_list)

comp_order <- c("Progenitor_Like","Myeloid_Like","Erythroid_Like","Lymphoid_Like","TNK_Like","Healthy_Like")
dfA3 <- dfA3 %>% mutate(compartment = factor(compartment, levels = intersect(comp_order, unique(compartment))),
                         group = factor(group, levels = c("Healthy (matched lineage)", "MPAL (malignant)")))

pA3 <- ggplot(dfA3, aes(x = compartment, y = score, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.85, position = position_dodge(0.8), width = 0.75) +
  geom_boxplot(width = 0.15, outlier.shape = NA, position = position_dodge(0.8), fill = "white", alpha = 0.6) +
  scale_fill_manual(values = c("Healthy (matched lineage)" = "grey65", "MPAL (malignant)" = "#e6550d"), name = NULL) +
  labs(title = "Granja et al. 2019: LT-PLC score, MPAL compartment vs lineage-matched healthy cells",
       subtitle = "up_all signature (UCell); healthy counterparts: Progenitor->HSC/CMP.LMPP, Myeloid->GMP/monocytes,\nErythroid->Eryth/Baso, Lymphoid->CLP/B, TNK & Healthy_Like->CD4/CD8/NK",
       x = NULL, y = "UCell score") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "granja_lineage_matched_comparison.png"), pA3, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(dir_out, "granja_lineage_matched_comparison.pdf"), pA3, width = 10, height = 6.5)
cat("\nSaved results/granja_lineage_matched_comparison.png / .pdf\n")

## ================= Analysis 3: driver genes =================
## Erythroid_Like and Lymphoid_Like are the two compartments elevated over
## their matched healthy counterpart (Analysis 2) -- identify which genes
## drive that difference, and whether it's a few outlier genes or broad.
cat("\n=== Analysis 3: driver genes for Lymphoid_Like / Erythroid_Like elevation ===\n")

human_sets <- readRDS(file.path(dir_data, "human_signature_sets.rds"))
sig_genes <- intersect(human_sets$up_all, rownames(seu))
expr <- GetAssayData(seu, layer = "data")[sig_genes, ]

driver_analysis <- function(mpal_label, healthy_labels, top_n = 25) {
  mpal_cells <- meta$cell[meta$mpal_class == mpal_label]
  healthy_cells <- meta$cell[meta$fine_class %in% healthy_labels]
  mean_mpal <- Matrix::rowMeans(expr[, mpal_cells, drop = FALSE])
  mean_healthy <- Matrix::rowMeans(expr[, healthy_cells, drop = FALSE])
  delta <- mean_mpal - mean_healthy
  list(
    top = data.frame(gene = names(delta), mean_mpal = mean_mpal, mean_healthy = mean_healthy, delta = delta) %>%
      arrange(desc(delta)) %>% slice_head(n = top_n) %>% mutate(compartment = mpal_label),
    all_delta = delta
  )
}

report_driver <- function(mpal_label, healthy_labels) {
  res <- driver_analysis(mpal_label, healthy_labels, top_n = 25)
  cat(sprintf("\nTop genes driving %s elevation over {%s}:\n", mpal_label, paste(healthy_labels, collapse = ",")))
  print(res$top %>% select(gene, mean_mpal, mean_healthy, delta), row.names = FALSE)
  pos <- res$all_delta[res$all_delta > 0]
  top10sum <- sum(sort(pos, decreasing = TRUE)[1:10])
  cat(sprintf("%s: top10 genes account for %.1f%% of total positive delta (n=%d genes with positive delta)\n",
              mpal_label, 100 * top10sum / sum(pos), length(pos)))
  res$top
}

d_lymph <- report_driver("Lymphoid_Like", lineage_map[["Lymphoid_Like"]])
d_eryth <- report_driver("Erythroid_Like", lineage_map[["Erythroid_Like"]])

write.csv(d_lymph, file.path(dir_out, "granja_lymphoid_driver_genes.csv"), row.names = FALSE)
write.csv(d_eryth, file.path(dir_out, "granja_erythroid_driver_genes.csv"), row.names = FALSE)
cat("\nSaved results/granja_lymphoid_driver_genes.csv, granja_erythroid_driver_genes.csv\n")
