# ==============================================================================
# 90_style.R — Shared visual style for every thesis figure (and the PDF tables)
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# One place that defines the palette, the ggplot theme, and the save routine, so
# every figure looks like part of one system and matches the manuscript's tables
# and body font. Sourced by 91_figures.R; the manuscript sets the SAME font as
# its PDF mainfont so figures and tables share typography.
#
# COLOUR CHOICE — validated, not eyeballed, and re-validated whenever it changes.
# The separation floor is 8 (OKLab distance x100) under normal vision, deuteranopia
# and protanopia, simulated with Machado, Oliveira and Fernandes (2009). The check
# lives in scripts/92b_palette_check.R and is run by hand, not by run_all.R.
#
# History of the choice matters only for the floor it established: green/red, used
# in the earliest figures, separates by 5.5 under deuteranopia and was rejected.
#
# CURRENT PAIR (burgundy / light blue), chosen 10/08:
#   normal 36.3 | deuteranopia 31.7 | protanopia 40.5
# The categorical sequence adds pink and dark grey. Over those four colours all
# six pairwise combinations clear the floor; the tightest is light blue against
# pink under protanopia (8.3), which is why those two are NOT placed adjacent in
# PAL_CAT. Over the pairs actually drawn -- the subject pair, every adjacent pair
# of PAL_CAT, and the first three against one another -- the tightest is burgundy
# against dark grey under deuteranopia, at 10.7. 92b_palette_check.R now reads
# PAL_CAT from here rather than keeping its own copy, and stops below the floor.
# ============================================================================

suppressPackageStartupMessages({ library(ggplot2) })
for (p in c("ragg", "systemfonts"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")

# Body/figure font. Must match the manuscript PDF `mainfont`. The OTFs ship with
# TinyTeX (tex-gyre, tex-gyre-math) and are copied to ~/Library/Fonts so that
# ragg can see them too; change in ONE place (and mirror it in the manuscript
# YAML `mainfont`).
THESIS_FONT <- "TeX Gyre Termes"
if (!any(grepl(THESIS_FONT, systemfonts::system_fonts()$family, fixed = TRUE)))
  stop("font not installed: ", THESIS_FONT,
       " -- copy texgyretermes-*.otf from TinyTeX to ~/Library/Fonts")

# ---- Palette ----------------------------------------------------------------
PAL <- list(
  mathematics = "#8C2D04",   # burgundy   — carries the results, so the heavier colour
  language    = "#6BAED6",   # light blue
  ink         = "#1A1A1A",   # primary text
  muted       = "#6B6B6B",   # axes, secondary text
  grid        = "#E6E6E6",   # recessive gridlines
  zero        = "#9A9A9A",   # reference line at 0
  fill_alpha  = 0.16         # confidence-band opacity
)
# Fixed categorical order for any 2+ category chart (assigned in order, never
# cycled). Ordered so that no adjacent pair sits near the floor: light blue and
# pink are separated by dark grey.
PAL_CAT <- c("#8C2D04", "#6BAED6", "#636363", "#E7969C", "#08519C", "#BDBDBD")

# Make PAL_CAT the DEFAULT discrete palette, so descriptive plots that do not set
# an explicit colour scale (03, 05) inherit it. Plots with an explicit scale
# (91_figures.R) are unaffected; continuous fills (e.g. the provincial-timing
# heatmap) use a sequential scale and are untouched.
options(ggplot2.discrete.colour = PAL_CAT, ggplot2.discrete.fill = PAL_CAT)

# Colour follows the ENTITY, not its position: language is always light blue and
# mathematics always burgundy, wherever they appear. Secondary two-category splits
# reuse the same two hues in order.
SUBJECT_COLORS <- c(
  "Language" = PAL$language, "Mathematics" = PAL$mathematics,
  "Non-promotion" = PAL$language, "Dropout" = PAL$mathematics)

# ---- Theme ------------------------------------------------------------------
theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = THESIS_FONT) +
    theme(
      text             = element_text(colour = PAL$ink),
      # No in-panel title or subtitle: the ONLY title a figure carries is its
      # Quarto caption below the panel, so a title here would be a second one.
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      plot.caption     = element_text(colour = PAL$muted, size = rel(0.72), hjust = 0,
                                      margin = margin(t = 8)),
      axis.title       = element_text(colour = PAL$muted, size = rel(0.9)),
      axis.text        = element_text(colour = PAL$muted),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = PAL$grid, linewidth = 0.3),
      legend.position  = "top",
      legend.justification = "left",
      legend.title     = element_blank(),
      legend.text      = element_text(colour = PAL$ink),
      strip.text       = element_text(colour = PAL$ink, size = rel(0.92)),
      plot.margin      = margin(10, 14, 8, 10)
    )
}

# Scales that reuse the palette by name.
scale_subject_colour <- function(...) scale_colour_manual(values = SUBJECT_COLORS, ...)
scale_subject_fill   <- function(...) scale_fill_manual(values = SUBJECT_COLORS, ...)

# ---- Save -------------------------------------------------------------------
# One save routine: ragg PNG (renders the system font crisply), 300 dpi, white
# background, consistent sizing. All figures go to output/figures.
save_fig <- function(p, file, width = 8, height = 4.6) {
  ggsave(here::here("output", "figures", file), plot = p,
         device = ragg::agg_png, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
}

