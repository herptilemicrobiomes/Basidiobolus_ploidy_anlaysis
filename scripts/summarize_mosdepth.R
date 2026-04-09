#!/usr/bin/env Rscript
# Summarize mosdepth coverage data across strains
# Reads all results/mosdepth/*.mosdepth.summary.txt files
# Produces histogram and boxplot of per-scaffold coverage depth,
# distinguishing _region scaffolds from whole-scaffold entries.

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# ---- locate files ----
mosdepth_dir <- file.path("results", "mosdepth")
summary_files <- list.files(mosdepth_dir,
                            pattern = "\\.mosdepth\\.summary\\.txt$",
                            full.names = TRUE)

if (length(summary_files) == 0) {
  stop("No *.mosdepth.summary.txt files found in ", mosdepth_dir)
}

# ---- read and combine ----
read_summary <- function(path) {
  strain <- str_remove(basename(path), "\\.mosdepth\\.summary\\.txt$")
  df <- read_tsv(path, col_types = cols(
    chrom  = col_character(),
    length = col_double(),
    bases  = col_double(),
    mean   = col_double(),
    min    = col_double(),
    max    = col_double()
  ))
  df$strain <- strain
  df
}

dat <- bind_rows(lapply(summary_files, read_summary))

# ---- filter and classify ----
dat <- dat |>
  filter(!chrom %in% c("total", "total_region")) |>
  mutate(type = if_else(str_ends(chrom, "_region"), "region", "whole_scaffold"))

# column 4 (mean) is the depth of coverage
# dat$mean is already that column

# ---- output directory ----
out_dir <- file.path("results", "mosdepth_plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- colour palette ----
type_colors <- c(whole_scaffold = "#2166ac", region = "#d6604d")
type_labels  <- c(whole_scaffold = "Whole scaffold", region = "Region (BED)")

# ===========================================================================
# 1. Histogram – one panel per strain, two colours for scaffold type
# ===========================================================================
p_hist <- ggplot(dat, aes(x = mean, fill = type)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.65, colour = NA) +
  facet_wrap(~strain, scales = "free_y") +
  scale_fill_manual(values = type_colors, labels = type_labels,
                    name = "Scaffold type") +
  labs(
    title = "Distribution of per-scaffold mean coverage depth by strain",
    x     = "Mean depth of coverage",
    y     = "Count"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom"
  )

hist_out <- file.path(out_dir, "coverage_histogram.pdf")
ggsave(hist_out, p_hist, width = 14, height = 10)
message("Saved histogram -> ", hist_out)

# ===========================================================================
# 2. Box-and-whisker – one box per scaffold type, panels per strain
# ===========================================================================
p_box <- ggplot(dat, aes(x = type, y = mean, fill = type)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.5) +
  facet_wrap(~strain, scales = "free_y") +
  scale_fill_manual(values = type_colors, labels = type_labels,
                    name = "Scaffold type") +
  scale_x_discrete(labels = type_labels) +
  labs(
    title = "Per-scaffold mean coverage depth by strain and scaffold type",
    x     = "Scaffold type",
    y     = "Mean depth of coverage"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 20, hjust = 1),
    legend.position  = "none"
  )

box_out <- file.path(out_dir, "coverage_boxplot.pdf")
ggsave(box_out, p_box, width = 14, height = 10)
message("Saved boxplot  -> ", box_out)
