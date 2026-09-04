# ==============================================================================
# 91_figures.R — All manuscript figures, from saved results, in one style
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Reads the plot-data CSVs written by the estimation scripts (output/tables)
# and renders every figure with the shared style (90_style.R). Figures live here,
# not inside the estimation scripts, so that:
#   - a figure and its companion table always come from the same run;
#   - the visual style can be changed without re-estimating anything;
#   - every figure is guaranteed to share palette, theme and font.
#
# Runs in seconds (no estimation). Writes to output/figures with the SAME
# filenames the manuscript references, so the PDF picks them up unchanged.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
source(here::here("scripts", "90_style.R"))
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

rd  <- function(f) fread(file.path(DIR_TABLES, f))
cap <- function(x) { x <- tolower(as.character(x))
  fifelse(x == "language", "Language", fifelse(x == "mathematics", "Mathematics", x)) }

# ---- Generic event-study figure ---------------------------------------------
# Expects columns k, est, lo, hi, subject. Reference line at 0; dashed marker
# just before adoption; band + line + points, coloured by subject.
# `support` is an optional named vector mapping event time to the number of
# switchers identifying that horizon; when supplied it is printed underneath the
# corresponding tick, because in a non-absorbing design the support collapses as
# the horizon grows and the reader cannot otherwise see it.
event_study <- function(dt, support = NULL,
                        ylab = "ATT (SD units)",
                        xlab = "Event time in outcome waves (0 = last untreated wave, reference)") {
  ks <- sort(unique(dt$k))
  labs_x <- if (is.null(support)) as.character(ks) else vapply(ks, function(k) {
    n <- support[as.character(k)]
    if (is.na(n)) as.character(k) else sprintf("%d\n(%s)", k, format(n, big.mark = ","))
  }, character(1))
  ggplot(dt, aes(k, est, colour = subject, fill = subject)) +
    geom_hline(yintercept = 0, colour = PAL$zero, linewidth = 0.4) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = PAL$zero, linewidth = 0.3) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = PAL$fill_alpha, colour = NA) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.9) +
    scale_subject_colour() + scale_subject_fill() +
    scale_x_continuous(breaks = ks, labels = labs_x) +
    labs(x = xlab, y = ylab) +
    theme_thesis()
}

# ============================================================================
# 1. Main event study (CDH): Language clean, mathematics placebo at k = -2
# ============================================================================
d <- rd("07_event_study_cdh.csv")[spec == "uncontrolled"]
d[, subject := cap(subject)]
# NB: figure text is kept ASCII-only. ragg did not render the en-dash or the "oe"
# ligature (they came out as dots), so the estimator name is spelled out. The
# restriction was established under an earlier THESIS_FONT and has not been
# retested against the current one; it is kept because ASCII costs nothing here.
# Switcher support per post-treatment horizon, read off the identifying-variation
# table rather than restated: k = 1 is Effect_1, k = 2 is Effect_2, k = 3 is
# Effect_3, and the counts collapse across them.
.sup <- rd("11d_switchers_by_horizon.csv")[sample == "full"]
sup_v <- setNames(.sup$switchers,
                  as.character(as.integer(sub("Effect_", "", trimws(.sup$horizon)))))
save_fig(event_study(d, support = sup_v), "07_event_study_cdh.png", height = 5.0)

# ============================================================================
# 2. Dynamic profile (post-treatment path is erratic, not accumulating)
# ============================================================================
d <- rd("07c_dynamic_profile.csv"); setnames(d, "att", "est"); d[, subject := cap(subject)]
save_fig(event_study(d), "07c_dynamic_profile.png")

# ============================================================================
# 3. Genuine 6h+ arm event study (five-hour schools separated)
# ============================================================================
d <- rd("08_event_study_genuine.csv")[spec == "stable (school-level)"]
d[, subject := cap(subject)]
save_fig(event_study(d), "08_event_study_genuine.png")

# ============================================================================
# 4. FEct dynamic (counterfactual-imputation estimator)
# ============================================================================
d <- rd("11_fect_dynamic.csv")
setnames(d, c("ATT", "CI.lower", "CI.upper", "spec"), c("est", "lo", "hi", "subject"))
d[, subject := cap(subject)]
save_fig(event_study(d), "11_fect_dynamic.png")

# ============================================================================
# 5. Grade-progression trajectories (12 pre-periods, composition check)
# ============================================================================
d <- rd("09_event_study_trajectories.csv")
d[, subject := fifelse(subject == "nonpromotion", "Non-promotion", "Dropout")]
save_fig(event_study(d,
  ylab = "Effect (percentage points)",
  xlab = "Event time (school years; 0 = first treated year)"),
  "09_event_study_trajectories.png", height = 4.8)

# ============================================================================
# 6. Adoption intensity (validates the binary indicator) — single series
# ============================================================================
d <- rd("06_adoption_intensity.csv")
p <- ggplot(d, aes(k, mean_share)) +
  geom_hline(yintercept = 0.5, linetype = "dotted", colour = PAL$zero, linewidth = 0.3) +
  geom_vline(xintercept = -0.5, linetype = "dashed", colour = PAL$zero, linewidth = 0.3) +
  geom_line(colour = PAL$language, linewidth = 0.8) +
  geom_point(colour = PAL$language, size = 1.9) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = scales::breaks_width(2)) +
  labs(x = "Years relative to first adoption", y = "Mean expanded-day share") +
  theme_thesis()
save_fig(p, "06_adoption_intensity.png")

# ============================================================================
# 7. Heterogeneity by subgroup (points + CI, coloured by subject, faceted)
# ============================================================================
d <- rd("10_heterogeneity_results.csv")[sample == "full" & !is.na(att)]
d[, `:=`(subject = cap(subject), group = as.character(group))]
p <- ggplot(d, aes(att, group, colour = subject)) +
  geom_vline(xintercept = 0, colour = PAL$zero, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = att - 1.96 * se, xmax = att + 1.96 * se),
                 height = 0.22, linewidth = 0.5,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 1.9, position = position_dodge(width = 0.5)) +
  facet_grid(dimension ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_subject_colour() +
  labs(x = "ATT (SD units)", y = NULL) +
  theme_thesis() +
  theme(panel.spacing = unit(10, "pt"), strip.placement = "outside")
save_fig(p, "10_heterogeneity.png", width = 8.5, height = 9)

# ============================================================================
# 8-9. Event studies comparing specifications (colour = spec, facet = subject)
# ============================================================================
# Uncontrolled vs controlled, for the appendix. Colour marks the SPECIFICATION
# here (not the subject), so the two panels share a legend; PAL_CAT keeps it
# CVD-safe.
event_study_by_spec <- function(dt) {
  dt[, subject := cap(subject)]
  # The dashed rule separates pre- from post-treatment, so it is drawn only when
  # there is something on both sides. On a pre-periods-only panel it would sit on
  # the right-hand edge separating nothing.
  onset <- if (any(dt$k >= 0))
    geom_vline(xintercept = -0.5, linetype = "dashed", colour = PAL$zero, linewidth = 0.3)
  ggplot(dt, aes(k, est, colour = spec, fill = spec)) +
    geom_hline(yintercept = 0, colour = PAL$zero, linewidth = 0.4) +
    onset +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.7) +
    facet_wrap(~ subject) +
    scale_colour_manual(values = PAL_CAT) + scale_fill_manual(values = PAL_CAT) +
    scale_x_continuous(breaks = scales::breaks_width(1)) +
    labs(x = "Event time (outcome waves; 0 = first treated wave)",
         y = "ATT (SD units)") +
    theme_thesis()
}

save_fig(event_study_by_spec(rd("07_event_study_cdh_controls.csv")),
  "07_event_study_cdh_controls.png", height = 4.4)

# PRE-TREATMENT PERIODS ONLY, and that is a decision rather than a crop. This
# event study exists to feed HonestDiD and to test whether BASELINE covariates
# absorb the pre-trend -- the one test that cannot be run on the primary estimator,
# whose `controls` option takes only time-varying covariates. Its POST-treatment
# coefficients are not reported: when the SES covariate moved to the predetermined
# construction, k=+2 changed sign in both subjects (mathematics +0.100 to -0.049,
# language +0.039 to -0.093) while the average effect moved 0.0047. A coefficient
# that flips sign under a decision that does not move the result cannot carry a
# reading, and plotting it invites one. Dynamics come from the primary estimator
# (07_event_study_cdh.png).
save_fig(event_study_by_spec(rd("07_event_study_cs_covariates.csv")[k < 0]),
  "07_event_study_cs_covariates.png", height = 4.4)


# ============================================================================
# 10. Adoption is not absorbing (design figure, body)
# ============================================================================
# The estimator choice rests on this and nothing else: the policy switches ON and
# OFF. Adoptions above the axis, contractions below, so the two directions are
# read at a glance rather than inferred from a reversal rate.
tr <- rd("05_transitions.csv")
trl <- rbindlist(list(
  tr[, .(year, n = adopt_0to1,  flow = "Adoptions (0 to 1)")],
  tr[, .(year, n = -revert_1to0, flow = "Contractions (1 to 0)")]))
trl[, flow := factor(flow, levels = c("Adoptions (0 to 1)", "Contractions (1 to 0)"))]
p <- ggplot(trl, aes(year, n, fill = flow)) +
  geom_col(width = 0.72) +
  geom_hline(yintercept = 0, colour = PAL$ink, linewidth = 0.4) +
  scale_fill_manual(values = c(PAL_CAT[1], PAL_CAT[2]), name = NULL) +
  scale_y_continuous(labels = function(x) format(abs(x), big.mark = ",")) +
  scale_x_continuous(breaks = scales::breaks_width(2)) +
  labs(x = NULL, y = "Schools (adoptions above the axis, contractions below)") +
  theme_thesis()
save_fig(p, "05_transitions.png", height = 4.4)

# ============================================================================
# 11. What the census counts as extension (measurement figure, body)
# ============================================================================
# The substantive measurement contribution: the 2023 expansion recorded by the
# census is mostly five-hour days, not genuine extension, while the pre-2023
# stock is almost entirely genuine. The federal programme funds BOTH schedules,
# so it names the cause here and never one of the two categories.
ct <- rd("08_contamination_2023.csv")
ctl <- melt(ct, id.vars = "when", measure.vars = c("genuine", "horamas", "none"),
            variable.name = "kind", value.name = "n")
ctl[, kind := factor(fifelse(kind == "genuine", "Genuine 6h+ extension",
                     fifelse(kind == "horamas", "Five-hour day", "Hours not recorded")),
                     levels = c("Genuine 6h+ extension", "Five-hour day",
                                "Hours not recorded"))]
ctl[, when := factor(fifelse(when %like% "2023 \\(flow", "Adopted in 2023 (flow)",
                             "Adopted before 2023 (stock)"),
                     levels = c("Adopted before 2023 (stock)", "Adopted in 2023 (flow)"))]
p <- ggplot(ctl, aes(when, n, fill = kind)) +
  geom_col(width = 0.6, position = position_stack(reverse = TRUE)) +
  # Only the segments wide enough to hold a label get one; the thin ones would
  # collide and are read off the axis instead.
  geom_text(aes(label = fifelse(n >= 300, format(n, big.mark = ","), "")),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            colour = "white", size = 3.1, family = THESIS_FONT) +
  scale_fill_manual(values = PAL_CAT[c(1, 2, 3)], name = NULL) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(x = NULL, y = "Schools") +
  theme_thesis()
save_fig(p, "08_contamination_2023.png", height = 4.6)


# ============================================================================
# 12. HonestDiD sensitivity, all four specifications in one panel
# ============================================================================
# Built here rather than in 07 (where the estimation lives) so it follows the
# project rule that figures are rendered in one place and in one style, and so
# the four small panels become one legible exhibit. The package's own plot also
# draws the "Original" confidence set, which 07 constructs but does not save; it
# is omitted rather than reconstructed, because the claim the text makes is about
# the robust sets and the bound at which they first cover zero.
hd <- rbindlist(lapply(
  c("language_uncontrolled", "language_baselinecovariates",
    "mathematics_uncontrolled", "mathematics_baselinecovariates"),
  function(k) {
    d <- rd(sprintf("07_honestdid_%s.csv", k))
    setnames(d, tolower(names(d)))
    parts <- strsplit(k, "_", fixed = TRUE)[[1]]
    d[, `:=`(subject = cap(parts[1]),
             spec = fifelse(parts[2] == "uncontrolled",
                            "Uncontrolled", "With baseline covariates"))]
  }))
hd[, spec := factor(spec, levels = c("Uncontrolled", "With baseline covariates"))]
p <- ggplot(hd, aes(factor(mbar), ymin = lb, ymax = ub, colour = subject)) +
  geom_hline(yintercept = 0, colour = PAL$zero, linewidth = 0.4) +
  geom_errorbar(width = 0.14, linewidth = 0.7) +
  facet_grid(subject ~ spec) +
  scale_subject_colour(guide = "none") +
  labs(x = expression(paste("Relative-magnitude bound  ", bar(M))),
       y = "Robust confidence set (SD units)") +
  theme_thesis()
save_fig(p, "07_honestdid_panel.png", width = 8, height = 5.6)
