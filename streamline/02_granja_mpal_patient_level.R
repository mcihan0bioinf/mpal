## 02_granja_mpal_patient_level.R
## Comment 1's Granja piece ("show MPAL as an independent cohort") plus
## Comments 2 and 3, all built on ONE UCell scoring pass of the Granja et al.
## 2019 atlas with the strict 1:1-ortholog LT-PLC up_lfc1 signature (same
## signature and same per-cell-UCell-then-per-patient-median approach as
## 01_entity_comparison.R, so Section A's numbers are directly comparable to
## that script's MPAL entity).
##
## MPAL5R is the relapse sample of patient MPAL5, not a 6th independent
## patient. All formal statistics use the primary n=5 (MPAL1-5); MPAL5R is
## shown separately as a sensitivity/longitudinal point.
##
## Section A (Comment 1): per-patient LT-PLC score, MPAL as its own cohort.
## Section B (Comment 2): patient x MPAL-compartment heterogeneity. The three
##   compartments with full coverage in every patient (Healthy_Like,
##   Progenitor_Like, Myeloid_Like) are repeated measures on the same 5
##   patients -> Friedman test + paired Wilcoxon post-hoc, not Kruskal-Wallis.
##   Sparser compartments (Lymphoid_Like, Erythroid_Like, TNK_Like) are shown
##   descriptively only.
## Section C (Comment 3): within-patient paired comparison, malignant
##   (Progenitor_Like / Myeloid_Like) vs that same patient's own Healthy_Like
##   cells -- paired Wilcoxon signed-rank, n=5, plus paired-line plots.
##
## Input:  data/granja_seurat.rds (celltype, sample)
##         data/human_signature_sets_strict.rds (1:1 orthologs; up_lfc1)
## Output: streamline/results/granja_lt_plc_lfc1_strict_percell.csv
##         streamline/results/granja_patient_cohort.png / .pdf         (Section A)
##         streamline/results/granja_compartment_heatmap.png / .pdf    (Section B)
##         streamline/results/granja_malignant_vs_healthy_paired.png / .pdf (Section C)

suppressPackageStartupMessages({
  library(Seurat); library(UCell); library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
})

## ---- exact permutation Friedman test for k=3 treatments ----
## At n=4 complete-case patients, the usual chi-squared asymptotic
## approximation behind friedman.test() is questionable; with only 3
## treatments the exact null distribution is small enough (6^n row-label
## permutations) to enumerate directly rather than approximate, so this is
## used instead of/alongside the asymptotic p-value for the headline result.
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
if (!dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)

human_strict <- readRDS(file.path(dir_data, "human_signature_sets_strict.rds"))
sig <- list(LT_PLC_up_lfc1 = human_strict$up_lfc1)
cat(sprintf("Signature: up_lfc1, strict 1:1 human orthologs, %d genes\n", length(sig$LT_PLC_up_lfc1)))

## ================= score Granja, per cell, batched =================
seu <- readRDS(file.path(dir_data, "granja_seurat.rds"))
DefaultAssay(seu) <- "RNA"
n_present <- sum(sig$LT_PLC_up_lfc1 %in% rownames(seu))
cat(sprintf("gene coverage: %d / %d (%.0f%%) signature genes detected in Granja\n\n",
            n_present, length(sig$LT_PLC_up_lfc1), 100 * n_present / length(sig$LT_PLC_up_lfc1)))

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
  score = as.numeric(scores[, "LT_PLC_up_lfc1_UCell"]),
  celltype = as.character(seu$celltype),
  sample = as.character(seu$sample)
) %>%
  mutate(patient = sub("_?T[0-9]+$", "", sample)) %>%
  filter(celltype %in% mpal_labels, grepl("^MPAL", patient)) %>%
  mutate(is_relapse = patient == "MPAL5R",
         patient_primary = ifelse(is_relapse, "MPAL5", patient))

write.csv(meta, file.path(dir_out, "granja_lt_plc_lfc1_strict_percell.csv"), row.names = FALSE)
cat("\nSaved streamline/results/granja_lt_plc_lfc1_strict_percell.csv\n\n")

rm(seu); gc(verbose = FALSE)

## primary patients = MPAL1-5 (relapse MPAL5R excluded from n)
primary_patients <- c("MPAL1","MPAL2","MPAL3","MPAL4","MPAL5")

## ================= Section A: MPAL as an independent cohort =================
cat("=== Section A: MPAL as an independent cohort (Comment 1) ===\n")

dfA <- meta %>% filter(celltype != "Healthy_Like") %>%
  group_by(patient, is_relapse) %>%
  summarise(median_score = median(score), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 20)

print(dfA)

dfA_primary <- dfA %>% filter(!is_relapse) %>% mutate(patient = factor(patient, levels = primary_patients))

pA <- ggplot(dfA_primary, aes(x = patient, y = median_score)) +
  geom_col(fill = "#c51b8a", alpha = 0.85, width = 0.6) +
  geom_point(data = dfA %>% filter(is_relapse), aes(x = "MPAL5", y = median_score),
             shape = 23, size = 3, fill = "white", color = "#c51b8a", stroke = 1.2) +
  labs(title = "A. LT-PLC score per MPAL patient (malignant/disease-like compartments pooled)",
       subtitle = "bars = primary cohort (MPAL1-5, median UCell score, disease-like compartments, ≥20 cells);\nopen diamond = MPAL5R (relapse sample of MPAL5, shown for reference, not a 6th independent patient)",
       x = NULL, y = "Median UCell score (up_lfc1, strict orthologs)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "granja_patient_cohort.png"), pA, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(dir_out, "granja_patient_cohort.pdf"), pA, width = 7, height = 5.5)
cat("Saved streamline/results/granja_patient_cohort.png / .pdf\n\n")

## ================= Section B: heterogeneity across MPAL compartments =================
cat("=== Section B: patient x compartment heterogeneity (Comment 2) ===\n")

dfB_all <- meta %>%
  group_by(patient, patient_primary, is_relapse, celltype) %>%
  summarise(median_score = median(score), n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 20)

full_compartments <- c("Healthy_Like", "Progenitor_Like", "Myeloid_Like")

cat("\nPatient x compartment table (n_cells >= 20 only; * = full-coverage compartment used in Friedman test):\n")
print(dfB_all %>% mutate(full_cov = celltype %in% full_compartments) %>% arrange(celltype, patient))

## heatmap: all patients (primary labeled; relapse shown as its own row) x all compartments meeting threshold
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
  labs(title = "B1. LT-PLC score by MPAL patient x compartment",
       subtitle = "min 20 cells per patient x compartment; MPAL5R (relapse) shown as its own row, excluded from n=5 tests below",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid = element_blank())

## ---- Friedman test across the 3 full-coverage compartments, primary n=5 ----
cat("\n--- Repeated-measures test across Healthy_Like / Progenitor_Like / Myeloid_Like (MPAL1-5) ---\n")
wide_full <- dfB_all %>% filter(!is_relapse, celltype %in% full_compartments) %>%
  select(patient, celltype, median_score) %>%
  pivot_wider(names_from = celltype, values_from = median_score)
print(wide_full)

## friedman.test() does not handle missing cells safely: rank() defaults to
## na.last=TRUE, which silently ranks a missing value as if it were real
## (confirmed: it does not error or drop the row) -- so restrict to patients
## with a complete row across all 3 compartments, not just nrow(wide_full)==5.
if (all(full_compartments %in% colnames(wide_full))) {
  complete_rows <- stats::complete.cases(wide_full[, full_compartments])
  wide_complete <- wide_full[complete_rows, ]
  cat(sprintf("\nComplete-case patients for Friedman test: %d / %d (%s)\n",
              nrow(wide_complete), nrow(wide_full), paste(wide_complete$patient, collapse = ", ")))
  if (nrow(wide_complete) >= 3) {
    mat <- as.matrix(wide_complete[, full_compartments])
    fr <- friedman.test(mat)
    print(fr)

    fr_exact <- friedman_exact_k3(mat)
    cat(sprintf("\nExact permutation Friedman test (enumerated over all %d row-label permutations, k=3):\n  statistic = %.3f, exact p = %.4f\n",
                fr_exact$n_permutations, fr_exact$statistic, fr_exact$p.value))

    cat(sprintf("\nPaired Wilcoxon post-hoc (n=%d complete-case patients):\n", nrow(mat)))
    pairs <- combn(full_compartments, 2, simplify = FALSE)
    for (p in pairs) {
      wt <- wilcox.test(mat[, p[1]], mat[, p[2]], paired = TRUE)
      cat(sprintf("%s vs %s: median %.4f vs %.4f, paired Wilcoxon p=%.4f\n",
                  p[1], p[2], median(mat[, p[1]]), median(mat[, p[2]]), wt$p.value))
    }
  } else {
    cat("Fewer than 3 complete-case patients -- skipping Friedman test.\n")
  }
} else {
  cat("Not all 3 full compartments present -- skipping Friedman test.\n")
}

## paired-line plot across the 3 full-coverage compartments
dfB2 <- dfB_all %>% filter(!is_relapse, celltype %in% full_compartments) %>%
  mutate(celltype = factor(celltype, levels = full_compartments))

pB2 <- ggplot(dfB2, aes(x = celltype, y = median_score, group = patient, color = patient)) +
  geom_line(alpha = 0.7, linewidth = 0.8) +
  geom_point(size = 2.5) +
  labs(title = "B2. Per-patient trajectories across the 3 full-coverage compartments (MPAL1-5)",
       subtitle = "Friedman test (repeated measures, same patients across compartments) -- see console output",
       x = NULL, y = "Median UCell score", color = "Patient") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid.minor = element_blank())

p_secB <- pB1 / pB2 + plot_annotation(
  title = "Granja et al. 2019: LT-PLC heterogeneity across MPAL malignant states",
  theme = theme(plot.title = element_text(face = "bold", size = 13)))

ggsave(file.path(dir_out, "granja_compartment_heatmap.png"), p_secB, width = 9, height = 10, dpi = 300)
ggsave(file.path(dir_out, "granja_compartment_heatmap.pdf"), p_secB, width = 9, height = 10)
cat("\nSaved streamline/results/granja_compartment_heatmap.png / .pdf\n\n")

## ================= Section C: malignant vs own healthy-like, paired =================
cat("=== Section C: malignant vs own Healthy_Like, paired within patient (Comment 3) ===\n")

paired_compartments <- c("Progenitor_Like", "Myeloid_Like")
dfC_wide <- dfB_all %>% filter(!is_relapse, celltype %in% c(paired_compartments, "Healthy_Like")) %>%
  select(patient, celltype, median_score) %>%
  pivot_wider(names_from = celltype, values_from = median_score)
print(dfC_wide)

results_C <- list()
for (comp in paired_compartments) {
  if (all(c(comp, "Healthy_Like") %in% colnames(dfC_wide)) &&
      sum(!is.na(dfC_wide[[comp]]) & !is.na(dfC_wide$Healthy_Like)) >= 2) {
    ok <- !is.na(dfC_wide[[comp]]) & !is.na(dfC_wide$Healthy_Like)
    a <- dfC_wide[[comp]][ok]; b <- dfC_wide$Healthy_Like[ok]
    wt <- wilcox.test(a, b, paired = TRUE)
    cat(sprintf("\n%s vs own Healthy_Like (n=%d patients): median %.4f vs %.4f, paired Wilcoxon p=%.4f\n",
                comp, sum(ok), median(a), median(b), wt$p.value))
    results_C[[comp]] <- data.frame(patient = dfC_wide$patient[ok], compartment = comp,
                                     malignant = a, healthy_like = b)
  }
}
dfC_long <- bind_rows(results_C) %>%
  pivot_longer(cols = c(malignant, healthy_like), names_to = "group", values_to = "score") %>%
  mutate(group = factor(group, levels = c("healthy_like", "malignant"),
                         labels = c("Own Healthy_Like", "Malignant")),
         compartment = factor(compartment, levels = paired_compartments))

## MPAL5R sensitivity: same paired comparison including the relapse sample
cat("\n--- MPAL5R (relapse) sensitivity check ---\n")
dfC_relapse <- dfB_all %>% filter(is_relapse, celltype %in% c(paired_compartments, "Healthy_Like")) %>%
  select(celltype, median_score, n_cells)
print(dfC_relapse)

pC <- ggplot(dfC_long, aes(x = group, y = score)) +
  geom_line(aes(group = patient), color = "grey60", alpha = 0.7) +
  geom_point(aes(color = patient), size = 2.8) +
  facet_wrap(~compartment) +
  labs(title = "C. LT-PLC score: malignant compartment vs that patient's own Healthy_Like cells",
       subtitle = "n=5 patients (MPAL1-5), paired Wilcoxon signed-rank per compartment -- see console output for p-values",
       x = NULL, y = "Median UCell score", color = "Patient") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "granja_malignant_vs_healthy_paired.png"), pC, width = 9, height = 5.5, dpi = 300)
ggsave(file.path(dir_out, "granja_malignant_vs_healthy_paired.pdf"), pC, width = 9, height = 5.5)
cat("\nSaved streamline/results/granja_malignant_vs_healthy_paired.png / .pdf\n")

cat("\nDone.\n")
