## 04_target_bulk_lfc2.R
## Robustness check on 03_target_bulk_entity_comparison.R: same TARGET bulk
## cohort, same patients, same pipeline -- just the stricter up_lfc2 signature
## (padj<0.05 & LFC>=2, strict 1:1 orthologs, 797 genes vs up_lfc1's 1575)
## instead of up_lfc1. Reuses the cached expression matrix from script 03
## (data/target_bulk/target_bulk_expr.rds) -- no re-download.
##
## Score: leukemia-only z-score reference throughout (the version 03 determined
## was not confounded by pooling Normal with variable-purity tumor samples).
##
## Input:  data/target_bulk/target_bulk_expr.rds, target_bulk_meta.rds (from 03)
##         data/human_signature_sets_strict.rds (1:1 orthologs; up_lfc2)
## Output: streamline/results/target_entity_comparison_lfc2.png / .pdf
##         streamline/results/target_entity_scores_lfc2.csv

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

proj_root <- getwd()
dir_data  <- file.path(proj_root, "data", "target_bulk")
dir_out   <- file.path(proj_root, "streamline", "results")

stopifnot(file.exists(file.path(dir_data, "target_bulk_expr.rds")))
expr <- readRDS(file.path(dir_data, "target_bulk_expr.rds"))
meta <- readRDS(file.path(dir_data, "target_bulk_meta.rds"))

human_strict <- readRDS(file.path(proj_root, "data", "human_signature_sets_strict.rds"))
sig_genes <- human_strict$up_lfc2
n_present <- sum(sig_genes %in% rownames(expr))
cat(sprintf("Signature: up_lfc2, strict 1:1 human orthologs, %d genes (%d present in TARGET expression matrix, %.0f%%)\n",
            length(sig_genes), n_present, 100 * n_present / length(sig_genes)))

logexpr <- log2(expr + 1)

## leukemia-only z-score reference (Normal excluded) -- same convention as
## 03's confirmed-robust version
leuk_files <- meta$file_id[meta$entity != "Normal"]
logexpr_leuk <- logexpr[, leuk_files, drop = FALSE]
gene_var <- apply(logexpr_leuk, 1, var, na.rm = TRUE)
z <- (logexpr_leuk[gene_var > 0, ] - rowMeans(logexpr_leuk[gene_var > 0, ], na.rm = TRUE)) /
  apply(logexpr_leuk[gene_var > 0, ], 1, sd, na.rm = TRUE)
score_genes <- intersect(sig_genes, rownames(z))
patient_score <- colMeans(z[score_genes, , drop = FALSE], na.rm = TRUE)

df <- meta %>% mutate(score = patient_score[file_id]) %>% filter(!is.na(score), entity != "Normal")
write.csv(df, file.path(dir_out, "target_entity_scores_lfc2.csv"), row.names = FALSE)
cat(sprintf("\nSaved streamline/results/target_entity_scores_lfc2.csv (%d patients)\n", nrow(df)))

entity_order <- c("AML", "B-ALL", "T-My-MPAL", "B-My-MPAL", "T-ALL")
entity_pal <- c(AML = "#e6550d", `B-ALL` = "#3182bd", `T-My-MPAL` = "#756bb1",
                 `B-My-MPAL` = "#c51b8a", `T-ALL` = "#31a354")
df_plot <- df %>% mutate(entity = factor(entity, levels = entity_order))

cat("\nMedian score by entity (up_lfc2, leukemia-only reference):\n")
print(df_plot %>% group_by(entity) %>% summarise(n = n(), median_score = median(score)) %>% arrange(desc(median_score)))

p <- ggplot(df_plot, aes(x = entity, y = score, fill = entity)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85, width = 0.6) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
  scale_fill_manual(values = entity_pal, guide = "none") +
  labs(title = "TARGET bulk RNA-seq: LT-PLC score across leukemia entities (up_lfc2 robustness check)",
       subtitle = sprintf("mean gene-wise z-score (log2 TPM+1), up_lfc2 signature (strict 1:1 orthologs, %d genes, %d present); leukemia-only z-score reference; one point = one patient",
                           length(sig_genes), n_present),
       x = NULL, y = "Mean z-score (LT-PLC up_lfc2)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold", size = 12.5),
        plot.subtitle = element_text(size = 8.5, color = "grey30"),
        panel.grid.minor = element_blank())

ggsave(file.path(dir_out, "target_entity_comparison_lfc2.png"), p, width = 9, height = 6, dpi = 300)
ggsave(file.path(dir_out, "target_entity_comparison_lfc2.pdf"), p, width = 9, height = 6)
cat("Saved streamline/results/target_entity_comparison_lfc2.png / .pdf\n")

cat("\nKruskal-Wallis across leukemia entities (up_lfc2, leukemia-only reference):\n")
print(kruskal.test(score ~ entity, data = df_plot))

cat("\nT-My-MPAL vs T-ALL (up_lfc2):\n")
print(wilcox.test(df_plot$score[df_plot$entity == "T-My-MPAL"], df_plot$score[df_plot$entity == "T-ALL"]))

cat("\nT-My-MPAL vs AML (up_lfc2):\n")
print(wilcox.test(df_plot$score[df_plot$entity == "T-My-MPAL"], df_plot$score[df_plot$entity == "AML"]))

cat("\nB-ALL vs T-ALL (up_lfc2):\n")
print(wilcox.test(df_plot$score[df_plot$entity == "B-ALL"], df_plot$score[df_plot$entity == "T-ALL"]))

cat("\nT-My-MPAL vs B-My-MPAL (up_lfc2):\n")
print(wilcox.test(df_plot$score[df_plot$entity == "T-My-MPAL"], df_plot$score[df_plot$entity == "B-My-MPAL"]))

## ---- compare ranking with the up_lfc1 result (rank correlation of per-patient scores) ----
lfc1_scores <- read.csv(file.path(dir_out, "target_entity_scores.csv")) %>%
  filter(entity != "Normal") %>% select(file_id, score_lfc1 = score_leuk_only)
df_compare <- df %>% select(file_id, entity, score_lfc2 = score) %>% inner_join(lfc1_scores, by = "file_id")
cat(sprintf("\nSpearman correlation between up_lfc1 and up_lfc2 per-patient scores (n=%d, matched patients):\n", nrow(df_compare)))
print(cor.test(df_compare$score_lfc1, df_compare$score_lfc2, method = "spearman"))

cat("\nDone.\n")
