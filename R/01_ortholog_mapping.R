## 01_ortholog_mapping.R
## Build the LT-PLC mouse signature gene set variants from the DESeq2
## diff-expression table (LT-PLC vs non-LT-PLC, LSK primary leukemia, n=5),
## then map each set to human orthologs in two modes:
##   - strict:  1:1 orthologs only (conservative, unambiguous)
##   - maximal: every valid ortholog hit kept, incl. one-to-many/paralogs
##              (maximizes gene coverage; downstream scoring uses this
##              unless a gene set specifically needs the strict version)
##
## Input:  Diff_genes_annotated.csv (project root; not tracked in git)
##         columns: mgi_symbol, log2FoldChange, padj (+ others)
## Output: data/mouse_signature_sets.rds
##         data/human_signature_sets.rds         (maximal mapping; use this downstream)
##         data/human_signature_sets_strict.rds  (strict 1:1 mapping)
##         results/ortholog_mapping_summary.csv
##         results/ortholog_mapping_figure.png / .pdf

set.seed(1234)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gprofiler2)
})

proj_root <- getwd()
stopifnot("Diff_genes_annotated.csv" %in% list.files(proj_root))
dir_data <- file.path(proj_root, "data")
dir_out  <- file.path(proj_root, "results")
for (d in c(dir_data, dir_out)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

## ---- 1. Load DEG table and build mouse gene set variants ----

deg <- read.csv2(file.path(proj_root, "Diff_genes_annotated.csv"),
                  row.names = 1, stringsAsFactors = FALSE)
# read.csv2: ';' sep, ',' decimal (matches the file's EU-locale export)

stopifnot(all(c("log2FoldChange", "padj", "mgi_symbol") %in% names(deg)))

deg <- deg %>%
  filter(!is.na(mgi_symbol), mgi_symbol != "", !is.na(padj), !is.na(log2FoldChange)) %>%
  distinct(mgi_symbol, .keep_all = TRUE)

cat(sprintf("Loaded Diff_genes_annotated.csv: %d genes with usable mgi_symbol/padj/log2FC\n", nrow(deg)))

sig_sets <- list(
  up_all    = deg %>% filter(log2FoldChange > 0, padj < 0.05) %>% pull(mgi_symbol),
  up_lfc1   = deg %>% filter(log2FoldChange >= 1, padj < 0.05) %>% pull(mgi_symbol),
  up_lfc2   = deg %>% filter(log2FoldChange >= 2, padj < 0.05) %>% pull(mgi_symbol),
  up_top50  = deg %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange)) %>% slice_head(n = 50)  %>% pull(mgi_symbol),
  up_top100 = deg %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange)) %>% slice_head(n = 100) %>% pull(mgi_symbol),
  up_top500 = deg %>% filter(padj < 0.05) %>% arrange(desc(log2FoldChange)) %>% slice_head(n = 500) %>% pull(mgi_symbol),
  down_lfc1 = deg %>% filter(log2FoldChange <= -1, padj < 0.05) %>% pull(mgi_symbol)
)
for (nm in names(sig_sets)) cat(sprintf("  %-10s n = %d\n", nm, length(sig_sets[[nm]])))

saveRDS(sig_sets, file.path(dir_data, "mouse_signature_sets.rds"))
saveRDS(deg,       file.path(dir_data, "deg_table_clean.rds"))

## ---- 2. Map each set to human orthologs via gprofiler2::gorth ----
## Two modes computed from the same gorth() call:
##   strict  = keep only mouse genes with exactly one human ortholog, and
##             only human genes claimed by exactly one mouse gene
##   maximal = keep every valid ortholog hit (all paralogs / one-to-many
##             mappings included); maximizes the number of genes retained

clean_ortholog_names <- function(x) {
  x <- unique(x)
  x[x != "N/A" & !is.na(x) & x != ""]
}

map_orthologs <- function(genes, set_name) {
  if (length(genes) == 0) return(list(strict = character(0), maximal = character(0)))
  orth <- gorth(query = genes, source_organism = "mmusculus", target_organism = "hsapiens",
                mthreshold = Inf, filter_na = TRUE)

  maximal <- clean_ortholog_names(orth$ortholog_name)

  n_per_mouse <- table(orth$input)
  one_to_one_mouse <- names(n_per_mouse[n_per_mouse == 1])
  orth1 <- orth %>% filter(input %in% one_to_one_mouse)
  n_per_human <- table(orth1$ortholog_name)
  ambiguous_human <- names(n_per_human[n_per_human > 1])
  orth1 <- orth1 %>% filter(!ortholog_name %in% ambiguous_human)
  strict <- clean_ortholog_names(orth1$ortholog_name)

  cat(sprintf("  %-10s: %d mouse genes -> %d human orthologs (strict 1:1) / %d human orthologs (maximal)\n",
              set_name, length(genes), length(strict), length(maximal)))
  list(strict = strict, maximal = maximal)
}

mapped <- lapply(names(sig_sets), function(nm) map_orthologs(sig_sets[[nm]], nm))
names(mapped) <- names(sig_sets)

human_sets_strict  <- lapply(mapped, `[[`, "strict")
human_sets_maximal <- lapply(mapped, `[[`, "maximal")

saveRDS(human_sets_maximal, file.path(dir_data, "human_signature_sets.rds"))
saveRDS(human_sets_strict,  file.path(dir_data, "human_signature_sets_strict.rds"))

summary_df <- data.frame(set = names(sig_sets),
                          n_mouse = sapply(sig_sets, length),
                          n_human_strict  = sapply(human_sets_strict, length),
                          n_human_maximal = sapply(human_sets_maximal, length))
write.csv(summary_df, file.path(dir_out, "ortholog_mapping_summary.csv"), row.names = FALSE)
print(summary_df)

## ---- 3. mapping-yield figure: strict vs maximal, per gene-set variant ----

set_labels <- c(
  up_all    = "Up: all\n(padj<0.05)",
  up_lfc1   = "Up: LFC>=1",
  up_lfc2   = "Up: LFC>=2",
  up_top50  = "Up: top 50",
  up_top100 = "Up: top 100",
  up_top500 = "Up: top 500",
  down_lfc1 = "Down: LFC<=-1"
)
set_order <- c("up_all", "up_lfc1", "up_lfc2", "up_top50", "up_top100", "up_top500", "down_lfc1")

plot_df <- summary_df %>%
  mutate(set = factor(set, levels = set_order, labels = set_labels[set_order])) %>%
  select(set, n_mouse, strict = n_human_strict, maximal = n_human_maximal) %>%
  pivot_longer(cols = c(strict, maximal), names_to = "mode", values_to = "n_human") %>%
  mutate(
    mode  = factor(mode, levels = c("strict", "maximal"),
                    labels = c("Strict (1:1 only)", "Maximal (all orthologs)")),
    pct   = 100 * n_human / n_mouse,
    label = sprintf("%d (%.0f%%)", n_human, pct)
  )

p <- ggplot(plot_df, aes(x = set, y = n_human, fill = mode)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "white", linewidth = 0.3) +
  geom_text(aes(label = label), position = position_dodge(width = 0.7),
            vjust = -0.4, size = 2.9, fontface = "bold") +
  scale_fill_manual(values = c("Strict (1:1 only)" = "#2c7fb8",
                                "Maximal (all orthologs)" = "#e6550d"),
                     name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Mouse-to-human ortholog mapping of the LT-PLC signature",
    subtitle = "gprofiler2::gorth; % = fraction of mouse genes in that set retained",
    x = NULL, y = "Number of human orthologs"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9.5, color = "grey30"),
    axis.text.x = element_text(size = 9)
  )

ggsave(file.path(dir_out, "ortholog_mapping_figure.pdf"), p, width = 8, height = 5.5)
ggsave(file.path(dir_out, "ortholog_mapping_figure.png"), p, width = 8, height = 5.5, dpi = 300)
cat("Saved results/ortholog_mapping_figure.pdf and .png\n")
