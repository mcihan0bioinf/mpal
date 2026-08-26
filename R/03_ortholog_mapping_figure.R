## 03_ortholog_mapping_figure.R
## Figure: mouse -> human ortholog mapping yield per LT-PLC signature
## gene-set variant, comparing strict 1:1 mapping vs. maximal mapping
## (all valid ortholog hits, incl. paralogs / one-to-many).
##
## Input:  results/ortholog_mapping_summary.csv (produced by 01_ortholog_mapping.R)
## Output: results/ortholog_mapping_figure.pdf
##         results/ortholog_mapping_figure.png

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

proj_root <- getwd()
dir_out   <- file.path(proj_root, "results")
summary_path <- file.path(dir_out, "ortholog_mapping_summary.csv")
stopifnot(file.exists(summary_path))

df <- read.csv(summary_path, stringsAsFactors = FALSE)
stopifnot(all(c("n_human_strict", "n_human_maximal") %in% names(df)))

set_labels <- c(
  up_all    = "Up: all\n(padj<0.05)",
  up_lfc1   = "Up: LFC>=1",
  up_lfc2   = "Up: LFC>=2",
  up_top50  = "Up: top 50",
  up_top100 = "Up: top 100",
  down_lfc1 = "Down: LFC<=-1"
)
set_order <- c("up_all", "up_lfc1", "up_lfc2", "up_top50", "up_top100", "down_lfc1")

df <- df %>%
  mutate(set = factor(set, levels = set_order, labels = set_labels[set_order]))

plot_df <- df %>%
  select(set, n_mouse, strict = n_human_strict, maximal = n_human_maximal) %>%
  pivot_longer(cols = c(strict, maximal), names_to = "mode", values_to = "n_human") %>%
  mutate(
    mode    = factor(mode, levels = c("strict", "maximal"),
                      labels = c("Strict (1:1 only)", "Maximal (all orthologs)")),
    pct     = 100 * n_human / n_mouse,
    label   = sprintf("%d (%.0f%%)", n_human, pct)
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
print(df %>% select(set, n_mouse, n_human_strict, n_human_maximal))
