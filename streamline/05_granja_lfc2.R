## 05_granja_lfc2.R
## Robustness check on 02_granja_mpal_patient_level.R's Section B (compartment
## heterogeneity) and Section C (malignant vs own Healthy_Like, paired) --
## same Granja cells, same method, just the stricter up_lfc2 signature
## (padj<0.05 & LFC>=2, strict 1:1 orthologs, 797 genes vs up_lfc1's 1575).
##
## Input:  data/granja_seurat.rds (celltype, sample)
##         data/human_signature_sets_strict.rds (1:1 orthologs; up_lfc2)
## Output: streamline/results/granja_compartment_heatmap_lfc2.png / .pdf
##         streamline/results/granja_malignant_vs_healthy_paired_lfc2.png / .pdf
##         streamline/results/granja_lt_plc_lfc2_strict_percell.csv

suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

## exact permutation Friedman test for k=3 treatments -- see 02_granja_mpal_patient_level.R
## for rationale (n=4 is too small to trust the asymptotic chi-squared approximation).
friedman_exact_k3 <- function(mat) {
  stopifnot(ncol(mat) == 3)
  n <- nrow(mat); k <- 3
  rank_mat <- t(apply(mat, 1, rank))
  stat <- function(R) (12 / (n * k * (k + 1))) * sum(R^2) - 3 * n * (k + 1)
  Q_obs <- stat(colSums(rank_mat))
  perms3 <- matrix(c(1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1), byrow = TRUE, ncol = 3)
  idx_grid <- as.matrix(expand.grid(rep(list(1:6), n)))
  Q_null <- apply(idx_grid, 1, function(idx) stat(colSums(perms3[idx, , drop = FALSE])))
  list(statistic = Q_obs, p.value = mean(Q_null >= Q_obs - 1e-8), n_permutations = nrow(idx_grid))
}

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data")
dir_out   <- file.path(proj_root, "streamline", "results")

human_strict <- readRDS(file.path(dir_data, "human_signature_sets_strict.rds"))
sig <- list(LT_PLC_up_lfc2 = human_strict$up_lfc2)
cat(sprintf("Signature: up_lfc2, strict 1:1 human orthologs, %d genes\n", length(sig$LT_PLC_up_lfc2)))

seu <- readRDS(file.path(dir_data, "granja_seurat.rds"))
DefaultAssay(seu) <- "RNA"
n_present <- sum(sig$LT_PLC_up_lfc2 %in% rownames(seu))
cat(sprintf("gene coverage: %d / %d (%.0f%%) signature genes detected in Granja\n\n",
            n_present, length(sig$LT_PLC_up_lfc2), 100 * n_present / length(sig$LT_PLC_up_lfc2)))

batch_size <- 10000
cells <- colnames(seu)
batches <- split(cells, ceiling(seq_along(cells) / batch_size))
scores_list <- vector("list", length(batches))
for (i in seq_along(batches)) {
  sub <- subset(seu, cells = batches[[i]])
  scores_list[[i]] <- ScoreSignatures_UCell(GetAssayData(sub, layer = "counts"), features = sig, maxRank = 3000)
  rm(sub); gc(verbose = FALSE)
  cat(sprintf("  batch %d/%d scored\n", i, length(batches)))
}
scores <- do.call(rbind, scores_list)[cells, , drop = FALSE]

mpal_labels <- c("Progenitor_Like","Lymphoid_Like","Myeloid_Like","Healthy_Like","Erythroid_Like","TNK_Like")
meta <- data.frame(
  cell = cells,
  score = as.numeric(scores[, "LT_PLC_up_lfc2_UCell"]),
  celltype = as.character(seu$celltype),
  sample = as.character(seu$sample)
) %>%
  mutate(patient = sub("_?T[0-9]+$", "", sample)) %>%
  filter(celltype %in% mpal_labels, grepl("^MPAL", patient)) %>%
  mutate(is_relapse = patient == "MPAL5R",
         patient_primary = ifelse(is_relapse, "MPAL5", patient))

write.csv(meta, file.path(dir_out, "granja_lt_plc_lfc2_strict_percell.csv"), row.names = FALSE)
rm(seu); gc(verbose = FALSE)

primary_patients <- c("MPAL1","MPAL2","MPAL3","MPAL4","MPAL5")

## ================= Section B: heterogeneity across MPAL compartments =================
cat("\n=== Section B: patient x compartment heterogeneity (up_lfc2) ===\n")

dfB_all <- meta %>%
  group_by(patient, patient_primary, is_relapse, celltype) %>%
  summarise(median_score = median(score), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 20)

full_compartments <- c("Healthy_Like", "Progenitor_Like", "Myeloid_Like")

state_order <- dfB_all %>% filter(!is_relapse) %>% group_by(celltype) %>%
  summarise(m = median(median_score)) %>% arrange(m) %>% pull(celltype)
patient_order <- c(primary_patients, "MPAL5R")

dfB_plot <- dfB_all %>% mutate(
  celltype = factor(celltype, levels = state_order),
  patient = factor(patient, levels = patient_order)
)

pB1 <- ggplot(dfB_plot, aes(x = celltype, y = patient, fill = median_score)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f\n(n=%d)", median_score, n_cells)), size = 2.6, lineheight = 0.85) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", name = "Median\nUCell score") +
  labs(title = "B1. LT-PLC score by MPAL patient x compartment (up_lfc2 robustness check)",
       subtitle = "min 20 cells per patient x compartment; MPAL5R (relapse) shown as its own row, excluded from n=5 tests below",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid = element_blank())

wide_full <- dfB_all %>% filter(!is_relapse, celltype %in% full_compartments) %>%
  select(patient, celltype, median_score) %>%
  pivot_wider(names_from = celltype, values_from = median_score)
cat("\n"); print(wide_full)

if (all(full_compartments %in% colnames(wide_full))) {
  complete_rows <- stats::complete.cases(wide_full[, full_compartments])
  wide_complete <- wide_full[complete_rows, ]
  cat(sprintf("\nComplete-case patients for Friedman test: %d / %d (%s)\n",
              nrow(wide_complete), nrow(wide_full), paste(wide_complete$patient, collapse = ", ")))
  if (nrow(wide_complete) >= 3) {
    mat <- as.matrix(wide_complete[, full_compartments])
    print(friedman.test(mat))
    fr_exact <- friedman_exact_k3(mat)
    cat(sprintf("\nExact permutation Friedman test (enumerated over all %d row-label permutations, k=3):\n  statistic = %.3f, exact p = %.4f\n",
                fr_exact$n_permutations, fr_exact$statistic, fr_exact$p.value))
    cat(sprintf("\nPaired Wilcoxon post-hoc (n=%d complete-case patients):\n", nrow(mat)))
    for (p in combn(full_compartments, 2, simplify = FALSE)) {
      wt <- wilcox.test(mat[, p[1]], mat[, p[2]], paired = TRUE)
      cat(sprintf("%s vs %s: median %.4f vs %.4f, paired Wilcoxon p=%.4f\n",
                  p[1], p[2], median(mat[, p[1]]), median(mat[, p[2]]), wt$p.value))
    }
  }
}

dfB2 <- dfB_all %>% filter(!is_relapse, celltype %in% full_compartments) %>%
  mutate(celltype = factor(celltype, levels = full_compartments))
pB2 <- ggplot(dfB2, aes(x = celltype, y = median_score, group = patient, color = patient)) +
  geom_line(alpha = 0.7, linewidth = 0.8) + geom_point(size = 2.5) +
  labs(title = "B2. Per-patient trajectories across the 3 full-coverage compartments (up_lfc2)",
       x = NULL, y = "Median UCell score", color = "Patient") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.5), panel.grid.minor = element_blank())

p_secB <- pB1 / pB2
ggsave(file.path(dir_out, "granja_compartment_heatmap_lfc2.png"), p_secB, width = 9, height = 10, dpi = 300)
ggsave(file.path(dir_out, "granja_compartment_heatmap_lfc2.pdf"), p_secB, width = 9, height = 10)
cat("\nSaved streamline/results/granja_compartment_heatmap_lfc2.png / .pdf\n")

## ================= Section C: malignant vs own healthy-like, paired =================
cat("\n=== Section C: malignant vs own Healthy_Like, paired (up_lfc2) ===\n")

paired_compartments <- c("Progenitor_Like", "Myeloid_Like")
dfC_wide <- dfB_all %>% filter(!is_relapse, celltype %in% c(paired_compartments, "Healthy_Like")) %>%
  select(patient, celltype, median_score) %>%
  pivot_wider(names_from = celltype, values_from = median_score)

results_C <- list()
for (comp in paired_compartments) {
  if (all(c(comp, "Healthy_Like") %in% colnames(dfC_wide))) {
    ok <- !is.na(dfC_wide[[comp]]) & !is.na(dfC_wide$Healthy_Like)
    a <- dfC_wide[[comp]][ok]; b <- dfC_wide$Healthy_Like[ok]
    if (sum(ok) >= 2) {
      wt <- wilcox.test(a, b, paired = TRUE)
      cat(sprintf("\n%s vs own Healthy_Like (n=%d patients): median %.4f vs %.4f, paired Wilcoxon p=%.4f\n",
                  comp, sum(ok), median(a), median(b), wt$p.value))
      results_C[[comp]] <- data.frame(patient = dfC_wide$patient[ok], compartment = comp,
                                       malignant = a, healthy_like = b)
    }
  }
}
dfC_long <- bind_rows(results_C) %>%
  pivot_longer(cols = c(malignant, healthy_like), names_to = "group", values_to = "score") %>%
  mutate(group = factor(group, levels = c("healthy_like", "malignant"), labels = c("Own Healthy_Like", "Malignant")),
         compartment = factor(compartment, levels = paired_compartments))

pC <- ggplot(dfC_long, aes(x = group, y = score)) +
  geom_line(aes(group = patient), color = "grey60", alpha = 0.7) +
  geom_point(aes(color = patient), size = 2.8) +
  facet_wrap(~compartment) +
  labs(title = "C. LT-PLC score: malignant compartment vs own Healthy_Like (up_lfc2 robustness check)",
       subtitle = "n=5 patients (MPAL1-5), paired Wilcoxon signed-rank per compartment -- see console output for p-values",
       x = NULL, y = "Median UCell score", color = "Patient") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        strip.text = element_text(face = "bold", size = 10), panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "granja_malignant_vs_healthy_paired_lfc2.png"), pC, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(dir_out, "granja_malignant_vs_healthy_paired_lfc2.pdf"), pC, width = 9, height = 5.5)
cat("\nSaved streamline/results/granja_malignant_vs_healthy_paired_lfc2.png / .pdf\n")

## ---- compare with up_lfc1 per-cell scores (rank correlation) ----
lfc1_path <- file.path(dir_out, "granja_lt_plc_lfc1_strict_percell.csv")
if (file.exists(lfc1_path)) {
  lfc1 <- read.csv(lfc1_path) %>% select(cell, score_lfc1 = score)
  merged <- meta %>% select(cell, score_lfc2 = score) %>% inner_join(lfc1, by = "cell")
  cat(sprintf("\nSpearman correlation between up_lfc1 and up_lfc2 per-cell scores (n=%d cells):\n", nrow(merged)))
  print(cor.test(merged$score_lfc1, merged$score_lfc2, method = "spearman"))
}

cat("\nDone.\n")
