## 03_ortholog_mapping_figure.R
## Figure: mouse -> human 1:1 ortholog mapping yield per LT-PLC signature
## gene-set variant.
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
  mutate(
    set      = factor(set, levels = set_order, labels = set_labels[set_order]),
    n_lost   = n_mouse - n_human,
    pct_kept = 100 * n_human / n_mouse
  )

plot_df <- df %>%
  select(set, n_human, n_lost) %>%
  pivot_longer(cols = c(n_human, n_lost), names_to = "category", values_to = "n") %>%
  mutate(category = factor(category, levels = c("n_lost", "n_human"),
                            labels = c("Lost (no 1:1 ortholog)", "Mapped to human (1:1)")))

pct_labels <- df %>% mutate(label = sprintf("%.0f%%", pct_kept))

p <- ggplot(plot_df, aes(x = set, y = n, fill = category)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(
    data = pct_labels,
    aes(x = set, y = n_mouse, label = label),
    inherit.aes = FALSE, vjust = -0.5, size = 3.4, fontface = "bold"
  ) +
  geom_text(
    data = df,
    aes(x = set, y = n_human / 2, label = n_human),
    inherit.aes = FALSE, color = "white", size = 3.2, fontface = "bold"
  ) +
  scale_fill_manual(values = c("Lost (no 1:1 ortholog)" = "grey80",
                                "Mapped to human (1:1)" = "#2c7fb8"),
                     name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Mouse-to-human ortholog mapping of the LT-PLC signature",
    subtitle = "1:1 orthologs only (gprofiler2::gorth); % = fraction of mouse genes retained",
    x = NULL, y = "Number of genes"
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

ggsave(file.path(dir_out, "ortholog_mapping_figure.pdf"), p, width = 7.5, height = 5.5)
ggsave(file.path(dir_out, "ortholog_mapping_figure.png"), p, width = 7.5, height = 5.5, dpi = 300)

cat("Saved results/ortholog_mapping_figure.pdf and .png\n")
print(df %>% select(set, n_mouse, n_human, n_lost, pct_kept))
