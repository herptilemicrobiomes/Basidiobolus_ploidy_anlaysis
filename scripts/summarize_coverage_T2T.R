#!/usr/bin/env Rscript
# summarize_coverage_T2T.R
#
# For each strain in samples.csv:
#   1. Read T2T scaffold list from lib/T2T_IDs_genomes/<strain>_t2t_ID.txt
#   2. Load results/mosdepth/<strain>.regions.bed.gz (BED5: chrom, start, end, name, coverage)
#   3. Filter to T2T scaffolds only
#   4. Normalise coverage to genome-wide mean (norm_cov = coverage / mean(all coverage))
#   5. Produce per-strain PDFs:
#        line plot    -continuous coverage profile, faceted by scaffold
#        impulse plot -vertical segments from 0 to coverage, faceted by scaffold
#        heatmap      -geom_rect tiles coloured by coverage, all scaffolds stacked
#      Each plot type is produced twice: linear scale and log10 scale.
#
# Output: results/mosdepth_T2T/<strain>/{line,impulse,heatmap}_{linear,log10}.pdf

library(ggplot2)
library(dplyr)
library(readr)
library(stringr)

# ---- paths ----
samples_file  <- "samples.csv"
t2t_dir       <- file.path("lib", "T2T_IDs_genomes")
mosdepth_dir  <- file.path("results", "mosdepth")
out_base      <- file.path("results", "mosdepth_T2T")
dir.create(out_base, showWarnings = FALSE, recursive = TRUE)

# ---- load strain IDs ----
samples <- read_csv(samples_file, col_types = cols(.default = col_character()))
strains <- samples[[1]]   # first column = Sample_ID

# ---- colour scales (shared) ----
cov_palette <- scale_fill_gradientn(
  colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
              "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"),
  name    = "Norm.\ncoverage"
)
cov_palette_log <- scale_fill_gradientn(
  colours = c("#313695", "#4575b4", "#74add1", "#abd9e9",
              "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"),
  trans   = "log10",
  name    = "Norm. cov\n(log10)"
)

# ---- helper: save with message ----
save_pdf <- function(plot, path, w = 14, h = NULL, nscaff = 1) {
  if (is.null(h)) h <- max(3, min(nscaff * 2.2, 24))
  ggsave(path, plot, width = w, height = h, limitsize = FALSE)
  message("  saved -> ", path)
}

# ---- shared theme ----
base_theme <- theme_bw(base_size = 10) +
  theme(
    strip.text      = element_text(face = "bold", size = 7),
    strip.background = element_rect(fill = "#e8eaf6"),
    panel.spacing   = unit(0.3, "lines"),
    legend.position = "right",
    axis.title      = element_text(size = 9)
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

  # ---- load regions BED ----
  bed_file <- file.path(mosdepth_dir, paste0(strain, ".regions.bed.gz"))
  if (!file.exists(bed_file)) {
    warning("BED file not found, skipping: ", bed_file)
    next
  }

  bed <- read_tsv(bed_file,
                  col_names = c("chrom", "start", "end", "name", "coverage"),
                  col_types = cols(
                    chrom    = col_character(),
                    start    = col_double(),
                    end      = col_double(),
                    name     = col_character(),
                    coverage = col_double()
                  ))

  # ---- normalise to genome-wide mean ----
  genome_mean <- mean(bed$coverage, na.rm = TRUE)
  message("  genome mean coverage: ", round(genome_mean, 2))
  bed <- bed |> mutate(norm_cov = coverage / genome_mean)

  # log10 floor for impulse/heatmap log plots (avoids log(0))
  log10_floor <- 0.001

  # ---- filter to T2T scaffolds and keep order ----
  dat <- bed |>
    filter(chrom %in% t2t_scaffolds) |>
    mutate(
      chrom          = factor(chrom, levels = t2t_scaffolds),
      midpos         = (start + end) / 2,          # window midpoint for line/impulse
      norm_cov_floor = pmax(norm_cov, log10_floor) # log-safe version
    )

  if (nrow(dat) == 0) {
    warning("No rows remain after T2T filter for ", strain)
    next
  }

  nscaff  <- length(unique(dat$chrom))
  out_dir <- file.path(out_base, strain)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # replace dots in strain name for plot titles
  strain_label <- str_replace_all(strain, "\\.", " ")

  # =========================================================================
  # 1. LINE PLOT  (continuous coverage profile)
  # =========================================================================
  make_line <- function(y_var, y_lab, y_scale, suffix) {
    p <- ggplot(dat, aes(x = midpos, y = .data[[y_var]])) +
      geom_line(colour = "#2166ac", linewidth = 0.4) +
      facet_wrap(~chrom, ncol = 1, scales = "free_x",
                 strip.position = "left") +
      y_scale +
      labs(
        title    = paste0(strain_label, " -T2T scaffold coverage (line)"),
        subtitle = "Coverage normalised to genome-wide mean",
        x        = "Genomic position (bp)",
        y        = y_lab
      ) +
      base_theme +
      theme(strip.text.y.left = element_text(angle = 0, hjust = 1))
    save_pdf(p, file.path(out_dir, paste0("line_", suffix, ".pdf")),
             nscaff = nscaff)
  }

  make_line("norm_cov",
            "Normalised coverage (×genome mean)",
            scale_y_continuous(limits = c(0, NA)),
            "linear")

  make_line("norm_cov",
            "Normalised coverage -log10",
            scale_y_log10(),
            "log10")

  # =========================================================================
  # 2. IMPULSE PLOT  (vertical segments from baseline to coverage value)
  # =========================================================================
  make_impulse <- function(y_var, y_lab, y_scale, suffix) {
    p <- ggplot(dat, aes(x = midpos)) +
      geom_segment(aes(xend = midpos, y = 0, yend = .data[[y_var]]),
                   colour = "#4575b4", linewidth = 0.3, alpha = 0.7) +
      facet_wrap(~chrom, ncol = 1, scales = "free_x",
                 strip.position = "left") +
      y_scale +
      labs(
        title    = paste0(strain_label, " -T2T scaffold coverage (impulse)"),
        subtitle = "Coverage normalised to genome-wide mean",
        x        = "Genomic position (bp)",
        y        = y_lab
      ) +
      base_theme +
      theme(strip.text.y.left = element_text(angle = 0, hjust = 1))
    save_pdf(p, file.path(out_dir, paste0("impulse_", suffix, ".pdf")),
             nscaff = nscaff)
  }

  make_impulse("norm_cov",
               "Normalised coverage (×genome mean)",
               scale_y_continuous(limits = c(0, NA)),
               "linear")

  # log10 impulse: segments from baseline 0 don't work on log scale,
  # so we draw from a small positive floor instead
  p_imp_log <- ggplot(dat, aes(x = midpos)) +
    geom_segment(aes(xend = midpos, y = log10_floor, yend = norm_cov_floor),
                 colour = "#4575b4", linewidth = 0.3, alpha = 0.7) +
    facet_wrap(~chrom, ncol = 1, scales = "free_x",
               strip.position = "left") +
    scale_y_log10() +
    labs(
      title    = paste0(strain_label, " -T2T scaffold coverage (impulse, log10)"),
      subtitle = "Coverage normalised to genome-wide mean",
      x        = "Genomic position (bp)",
      y        = "Normalised coverage -log10"
    ) +
    base_theme +
    theme(strip.text.y.left = element_text(angle = 0, hjust = 1))
  save_pdf(p_imp_log, file.path(out_dir, "impulse_log10.pdf"), nscaff = nscaff)

  # =========================================================================
  # 3. HEATMAP  (geom_rect: x = position, y = scaffold row, fill = coverage)
  #    All T2T scaffolds stacked in one panel — like a genome browser track.
  # =========================================================================
  make_heatmap <- function(fill_var, fill_scale, suffix) {
    p <- ggplot(dat, aes(xmin = start, xmax = end,
                         ymin = 0,    ymax = 1,
                         fill = .data[[fill_var]])) +
      geom_rect() +
      fill_scale +
      facet_wrap(~chrom, ncol = 1, scales = "free_x",
                 strip.position = "left") +
      labs(
        title    = paste0(strain_label, " -T2T scaffold coverage (heatmap)"),
        subtitle = "Coverage normalised to genome-wide mean; colour = depth",
        x        = "Genomic position (bp)",
        y        = NULL
      ) +
      base_theme +
      theme(
        strip.text.y.left  = element_text(angle = 0, hjust = 1),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        panel.grid         = element_blank()
      )
    # heatmap panels can be short — fixed height per scaffold
    h <- max(3, min(nscaff * 1.2 + 2, 24))
    ggsave(file.path(out_dir, paste0("heatmap_", suffix, ".pdf")),
           p, width = 16, height = h, limitsize = FALSE)
    message("  saved -> ", file.path(out_dir, paste0("heatmap_", suffix, ".pdf")))
  }

  make_heatmap("norm_cov",      cov_palette,     "linear")
  make_heatmap("norm_cov_floor", cov_palette_log, "log10")

  message("  done: ", nscaff, " T2T scaffolds plotted")
}

message("\nAll strains processed. Output in: ", out_base)
