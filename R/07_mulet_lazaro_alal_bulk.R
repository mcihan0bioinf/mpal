## 07_mulet_lazaro_alal_bulk.R
## Score the LT-PLC upregulated human signature per patient in Mulet-Lazaro
## et al. 2025's 30 ALAL (acute leukemia of ambiguous lineage) cohort, then
## group patients by the paper's own E-CAL transcriptional-classifier
## lineage assignment and WHO 2022 classification.
##
## Scoring method: mean of per-gene z-scores (z-scored across the 30
## patients on log2(TPM+1)) over the signature genes present in the data.
## This is a standard, dependency-light bulk single-sample scoring approach
## (avoids GSVA, whose current Bioconductor build has a broken
## SpatialExperiment/magick dependency chain on this machine).
##
## NOTE: this dataset's AML (n=145) / T-ALL (n=85) comparator cohorts from
## the same paper are EGA controlled-access and are NOT included here -- this
## is an ALAL-internal analysis only (patients grouped by their own
## expression-derived lineage call, not compared to an external cohort).
##
## Input:  data/mulet_lazaro_2025/TPM_salmon_ALAL.txt
##         data/mulet_lazaro_2025/hem370195-sup-0002-...supplementaltables.xlsx
##         data/human_signature_sets.rds (maximal mapping; up_all / top500)
## Output: results/mulet_lazaro_alal_signature_score.png / .pdf
##         results/mulet_lazaro_alal_patient_scores.csv

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(readxl)
})

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data", "mulet_lazaro_2025")
dir_out   <- file.path(proj_root, "results")

## ---- 1. expression matrix ----
tpm <- read.delim(file.path(dir_data, "TPM_salmon_ALAL.txt"), check.names = FALSE)
stopifnot(all(c("gene_name", "gene_id") %in% colnames(tpm)))
expr <- tpm %>% filter(!duplicated(gene_name)) %>% tibble::column_to_rownames("gene_name") %>%
  select(-gene_id) %>% as.matrix()
mode(expr) <- "numeric"
cat(sprintf("Expression matrix: %d genes x %d patients\n", nrow(expr), ncol(expr)))

## ---- 2. patient metadata (E-CAL lineage call, WHO classification) ----
xlsx <- file.path(dir_data, "hem370195-sup-0002-20250514_supplementaltables.xlsx")
s3 <- read_excel(xlsx, sheet = "S3 - Patient summary")[-1, ] %>%
  transmute(ID = as.character(ID), ECAL_S3 = ECAL, WHO_S3 = WHO)
s13 <- read_excel(xlsx, sheet = "S13 - Predictions")[-1, ] %>%
  transmute(ID = as.character(ID),
            ecal_label = `Label assigned by E-CAL`,
            who_class  = `WHO classification`)

meta <- s3 %>% left_join(s13, by = "ID") %>% filter(ID %in% colnames(expr))  # drops trailing footnote row
cat(sprintf("\nPatient metadata: %d patients\n", nrow(meta)))
cat("\nE-CAL label distribution:\n"); print(table(meta$ecal_label, useNA = "ifany"))
cat("\nWHO classification distribution:\n"); print(table(meta$who_class, useNA = "ifany"))

stopifnot(all(meta$ID %in% colnames(expr)))
expr <- expr[, meta$ID]  # order columns to match metadata

## ---- 3. score each patient: mean of per-gene z-scores on log2(TPM+1) ----
human_sets <- readRDS(file.path(proj_root, "data", "human_signature_sets.rds"))
sig <- list(LT_PLC_up_all = human_sets$up_all, LT_PLC_up_top500 = human_sets$up_top500)
for (nm in names(sig)) {
  n_present <- sum(sig[[nm]] %in% rownames(expr))
  cat(sprintf("%s: %d genes, %d present in ALAL expression matrix\n", nm, length(sig[[nm]]), n_present))
}

logexpr <- log2(expr + 1)
## z-score each gene across the 30 patients (only genes with nonzero variance)
gene_var <- apply(logexpr, 1, var)
z <- (logexpr[gene_var > 0, ] - rowMeans(logexpr[gene_var > 0, ])) / apply(logexpr[gene_var > 0, ], 1, sd)

score_signature <- function(genes) {
  genes <- intersect(genes, rownames(z))
  colMeans(z[genes, , drop = FALSE])
}
scores <- sapply(sig, score_signature)

df <- meta %>% bind_cols(as.data.frame(scores[meta$ID, ]))

## simplify E-CAL label: single-lineage calls vs multi-label ("AML;B-ALL" etc = ambiguous)
df <- df %>% mutate(
  ecal_simple = case_when(
    is.na(ecal_label) ~ NA_character_,
    grepl(";", ecal_label) ~ "Ambiguous (multi-label)",
    TRUE ~ ecal_label
  ),
  is_who_mpal = grepl("B", who_class) & grepl("T", who_class) | grepl("B", who_class) & grepl("M", who_class) | grepl("T", who_class) & grepl("M", who_class)
)

write.csv(df, file.path(dir_out, "mulet_lazaro_alal_patient_scores.csv"), row.names = FALSE)
cat("\nSaved results/mulet_lazaro_alal_patient_scores.csv\n")

## ---- 4. plot: per-patient score, ordered, colored by E-CAL lineage ----
plot_df <- df %>% filter(!is.na(ecal_simple)) %>%
  mutate(ID = reorder(ID, LT_PLC_up_all))

pal <- c("AML" = "#e6550d", "B-ALL" = "#3182bd", "T-ALL" = "#31a354",
         "Ambiguous (multi-label)" = "grey50")

p <- ggplot(plot_df, aes(x = ID, y = LT_PLC_up_all, fill = ecal_simple)) +
  geom_col(width = 0.7) +
  geom_point(aes(y = LT_PLC_up_top500), shape = 21, fill = "white", size = 1.8, stroke = 0.6) +
  scale_fill_manual(values = pal, name = "E-CAL lineage call") +
  labs(title = "Mulet-Lazaro et al. 2025: LT-PLC signature score across 30 ALAL patients",
       subtitle = "mean gene-wise z-score on bulk log2(TPM+1); bars = up_all signature (2380 genes), open circles = up_top500 signature (378 genes);\npatients ordered by up_all score, colored by the paper's own E-CAL transcriptional lineage classifier",
       x = "Patient", y = "Mean z-score (LT-PLC up)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "mulet_lazaro_alal_signature_score.png"), p, width = 11, height = 6.5, dpi = 300)
ggsave(file.path(dir_out, "mulet_lazaro_alal_signature_score.pdf"), p, width = 11, height = 6.5)
cat("Saved results/mulet_lazaro_alal_signature_score.png / .pdf\n")

## ---- 5. grouped comparison by lineage call ----
cat("\nMedian LT-PLC (up_all) score by E-CAL lineage call:\n")
print(df %>% filter(!is.na(ecal_simple)) %>% group_by(ecal_simple) %>%
        summarise(median = median(LT_PLC_up_all), n = n()) %>% arrange(desc(median)))

cat("\nMedian LT-PLC (up_all) score by WHO classification:\n")
print(df %>% filter(!is.na(who_class)) %>% group_by(who_class) %>%
        summarise(median = median(LT_PLC_up_all), n = n()) %>% arrange(desc(median)))

## ---- 6. significance tests (both groupings are small / likely underpowered) ----
cat("\nKruskal-Wallis, E-CAL lineage call (AML/B-ALL/T-ALL only):\n")
df_ecal <- df %>% filter(!is.na(ecal_simple), ecal_simple != "Ambiguous (multi-label)")
print(kruskal.test(LT_PLC_up_all ~ ecal_simple, data = df_ecal))

cat("\nKruskal-Wallis, WHO classification:\n")
df_who <- df %>% filter(!is.na(who_class))
print(kruskal.test(LT_PLC_up_all ~ who_class, data = df_who))

## ---- 7. plot: per-patient score, colored by WHO classification ----
plot_df_who <- df %>% filter(!is.na(who_class)) %>% mutate(ID = reorder(ID, LT_PLC_up_all))
pal_who <- c("M" = "#e6550d", "B/M" = "#9ecae1", "T" = "#31a354", "T/M" = "#a1d99b", "B/T" = "#756bb1")

p_who <- ggplot(plot_df_who, aes(x = ID, y = LT_PLC_up_all, fill = who_class)) +
  geom_col(width = 0.7) +
  geom_point(aes(y = LT_PLC_up_top500), shape = 21, fill = "white", size = 1.8, stroke = 0.6) +
  scale_fill_manual(values = pal_who, name = "WHO classification") +
  labs(title = "Mulet-Lazaro et al. 2025: LT-PLC signature score across 30 ALAL patients (WHO grouping)",
       subtitle = sprintf("same scores as previous figure, colored by WHO classification instead of E-CAL call; Kruskal-Wallis p = %.2f",
                           kruskal.test(LT_PLC_up_all ~ who_class, data = df_who)$p.value),
       x = "Patient", y = "Mean z-score (LT-PLC up)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        legend.position = "top",
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "mulet_lazaro_alal_signature_score_who.png"), p_who, width = 11, height = 6.5, dpi = 300)
ggsave(file.path(dir_out, "mulet_lazaro_alal_signature_score_who.pdf"), p_who, width = 11, height = 6.5)
cat("\nSaved results/mulet_lazaro_alal_signature_score_who.png / .pdf\n")
