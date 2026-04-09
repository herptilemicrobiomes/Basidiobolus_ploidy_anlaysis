#!/usr/bin/env Rscript
# summarize_coverage_T2T_quantized.R
#
# For each strain in samples.csv:
#   1. Read T2T scaffold list from lib/T2T_IDs_genomes/<strain>_t2t_ID.txt
#   2. Load results/mosdepth/<strain>.quantized.bed.gz (BED4: chrom, start, end, label)
#      Labels: NO_COVERAGE, LOW_COVERAGE, CALLABLE, HIGH_COVERAGE, VERY_HIGH_COVERAGE
#   3. Filter to T2T scaffolds only
#   4. Produce per-strain PDFs:
#        heatmap - geom_rect tiles coloured by coverage class, all scaffolds stacked
#        bar     - stacked bar showing proportion of each class per scaffold
#
# Output: results/mosdepth_T2T_quantized/<strain>/{heatmap,bar}.pdf

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(forcats)

# ---- paths ----
samples_file  <- "samples.csv"
t2t_dir       <- file.path("lib", "T2T_IDs_genomes")
mosdepth_dir  <- file.path("results", "mosdepth")
out_base      <- file.path("results", "mosdepth_T2T_quantized")
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

# ---- load strain IDs ----
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))
strains <- samples[[1]]   # first column = Sample_ID

# ---- quantized coverage levels: cool -> hot ----
cov_levels <- c("NO_COVERAGE", "LOW_COVERAGE", "CALLABLE",
                 "HIGH_COVERAGE", "VERY_HIGH_COVERAGE")

# cool-to-hot palette matching the five classes
cov_colours <- c(
  NO_COVERAGE       = "#313695",   # deep blue  (cool)
  LOW_COVERAGE      = "#74add1",   # light blue
  CALLABLE          = "#ffffbf",   # pale yellow (neutral)
  HIGH_COVERAGE     = "#f46d43",   # orange
  VERY_HIGH_COVERAGE = "#a50026"   # deep red   (hot)
)

# ---- helper: save with message ----
save_pdf <- function(plot, path, w = 14, h = NULL, nscaff = 1) {
  if (is.null(h)) h <- max(3, min(nscaff * 1.2 + 2, 24))
  ggsave(path, plot, width = w, height = h, limitsize = FALSE)
  message("  saved -> ", path)
}

# ---- shared theme ----
base_theme <- theme_bw(base_size = 10) +
  theme(
    strip.text       = element_text(face = "bold", size = 7),
    strip.background = element_rect(fill = "#e8eaf6"),
    panel.spacing    = unit(0.3, "lines"),
    legend.position  = "right",
    axis.title       = element_text(size = 9)
  )

# ===========================================================================
# Main loop
# ===========================================================================
for (strain in strains) {

  message("\n=== ", strain, " ===")

  # ---- locate T2T scaffold list ----
  t2t_file <- file.path(t2t_dir, paste0(strain, "_t2t_ID.txt"))
  if (!file.exists(t2t_file)) {
    warning("T2T ID file not found, skipping: ", t2t_file)
    next
  }
  t2t_scaffolds <- readLines(t2t_file)
  t2t_scaffolds <- t2t_scaffolds[nzchar(trimws(t2t_scaffolds))]

  # ---- load quantized BED ----
  bed_file <- file.path(mosdepth_dir, paste0(strain, ".quantized.bed.gz"))
  if (!file.exists(bed_file)) {
    warning("Quantized BED file not found, skipping: ", bed_file)
    next
  }

  bed <- read_tsv(bed_file,
                  col_names = c("chrom", "start", "end", "label"),
                  col_types = cols(
                    chrom = col_character(),
                    start = col_double(),
                    end   = col_double(),
                    label = col_character()
                  ))

  # ---- filter to T2T scaffolds and set factor levels ----
  dat <- bed |>
    filter(chrom %in% t2t_scaffolds) |>
    mutate(
      chrom = factor(chrom, levels = t2t_scaffolds),
      label = factor(label, levels = cov_levels),
      width = end - start
    )

  if (nrow(dat) == 0) {
    warning("No rows remain after T2T filter for ", strain)
    next
  }

  nscaff      <- length(unique(dat$chrom))
  out_dir     <- file.path(out_base, strain)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  strain_label <- str_replace_all(strain, "\\.", " ")

  # =========================================================================
  # 1. HEATMAP  (geom_rect tiles, one row per scaffold, colour = class)
  # =========================================================================
  p_heat <- ggplot(dat, aes(xmin = start, xmax = end,
                             ymin = 0,    ymax = 1,
                             fill = label)) +
    geom_rect() +
    scale_fill_manual(values = cov_colours,
                      drop   = FALSE,
                      name   = "Coverage\nclass") +
    facet_wrap(~chrom, ncol = 1, scales = "free_x",
               strip.position = "left") +
    labs(
      title    = paste0(strain_label, " - T2T scaffold coverage (quantized heatmap)"),
      subtitle = "Colour = mosdepth quantized coverage class",
      x        = "Genomic position (bp)",
      y        = NULL
    ) +
    base_theme +
    theme(
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      axis.text.y       = element_blank(),
      axis.ticks.y      = element_blank(),
      panel.grid        = element_blank()
    )

  save_pdf(p_heat, file.path(out_dir, "heatmap.pdf"), nscaff = nscaff)

  # =========================================================================
  # 2. STACKED BAR  (proportion of each class per scaffold)
  # =========================================================================
  bar_dat <- dat |>
    group_by(chrom, label) |>
    summarise(total_bp = sum(width), .groups = "drop") |>
    group_by(chrom) |>
    mutate(pct = total_bp / sum(total_bp) * 100) |>
    ungroup()

  p_bar <- ggplot(bar_dat, aes(x = pct, y = fct_rev(chrom), fill = label)) +
    geom_col(width = 0.8) +
    scale_fill_manual(values = cov_colours,
                      drop   = FALSE,
                      name   = "Coverage\nclass") +
    scale_x_continuous(expand = c(0, 0),
                       labels = function(x) paste0(x, "%")) +
    labs(
      title    = paste0(strain_label, " - T2T scaffold coverage class proportions"),
      subtitle = "Percentage of scaffold bases in each quantized class",
      x        = "Percentage of scaffold (bp)",
      y        = "Scaffold"
    ) +
    base_theme

  h_bar <- max(3, min(nscaff * 0.45 + 2, 20))
  ggsave(file.path(out_dir, "bar.pdf"), p_bar,
         width = 12, height = h_bar, limitsize = FALSE)
  message("  saved -> ", file.path(out_dir, "bar.pdf"))

  message("  done: ", nscaff, " T2T scaffolds plotted")
}

message("\nAll strains processed. Output in: ", out_base)
