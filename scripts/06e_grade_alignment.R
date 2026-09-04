# ==============================================================================
# 06e_grade_alignment.R — Treatment (whole-primary) vs outcome (grade 6)
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# The treatment is a WHOLE-PRIMARY count of extended/full-day students (we don't
# observe which grades); the outcome is grade-6 test scores. Two facts to quantify
# for the identification discussion:
#   (1) primary structure: 6- vs 7-year primary (grade 6 is last vs second-to-last);
#   (2) grade-6 alignment: an arithmetic (pigeonhole) bound -- grade 6 is CERTAINLY
#       treated when n_expanded + enroll_g6 > enroll_primary (the expanded students
#       cannot all fit in the other grades). We report how often this holds and
#       re-estimate the headline on the "grade-6-certain" subset (a la Hincapie).
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(did); library(DIDmultiplegtDYN); library(polars)})
# No seed: nothing here draws on the random number generator. The only estimator
# in this script is did_multiplegt_dyn with its analytic clustered variance.

panel <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
# Treatment comes from 02 (the single source of truth) instead of being rebuilt
# here. The persistence rule used to live in nine copies, and shift() silently
# depends on the panel being sorted -- an invariant no script declared.
panel <- merge(panel,
               readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))[
                 , .(school_id, year, main_treatment)],
               by = c("school_id", "year"), all.x = TRUE)
if (!"n_expanded" %in% names(panel)) panel[, n_expanded := n_extended_day + n_full_day]

# Grade-6 CERTAINLY treated (pigeonhole): expanded students exceed all non-g6 seats.
# `g6_computable` is kept separate so that "not certain" is never confused with
# "not measurable" when the frequency is reported.
panel[, g6_computable := !is.na(enroll_g6) & !is.na(enroll_primary) & !is.na(n_expanded)]
panel[, g6_certain := as.integer(g6_computable & (n_expanded + enroll_g6 > enroll_primary))]
# Primary structure: 7-year (grade 7 present) vs 6-year
panel[, has_g7 := as.integer(!is.na(enroll_g7) & enroll_g7 > 0)]

d <- panel[year %in% WAVES & !is.na(score_lang),
           .(school_id, year, y_lang = score_lang / 100, y_math = score_math / 100,
             main_treatment, g6_certain, g6_computable, has_g7)]
d <- d[!is.na(main_treatment)]
d[, `:=`(wave = match(year, WAVES), gid = as.integer(factor(school_id)))]

# (1) Primary structure (school-level: ever has grade 7?)
sch_g7 <- d[, .(has_g7 = as.integer(any(has_g7 == 1, na.rm = TRUE))), by = school_id]

# (2) THE BOUND: how often grade-6 exposure is arithmetically forced, among the
# treated school-years that the estimation actually uses. This is the quantity the
# empirical-strategy section claims -- that the gap between assignment and receipt
# is bounded rather than assumed -- and it used to be computed and discarded.
# It bounds WHETHER any grade-6 pupil is exposed, not WHAT SHARE of grade 6 is.
tr <- d[main_treatment == 1]
bound <- rbind(
  tr[, .(year = NA_integer_, school_years = .N, computable = sum(g6_computable),
         pct_certain = round(100 * mean(g6_certain[g6_computable]), 1))],
  tr[, .(school_years = .N, computable = sum(g6_computable),
         pct_certain = round(100 * mean(g6_certain[g6_computable]), 1)), by = year][order(year)],
  fill = TRUE)
if (nrow(bound) != length(WAVES) + 1L) stop("06e: expected one row per wave plus the total")
if (bound[1L, computable / school_years] < 0.95)
  stop("06e: the pigeonhole bound is not computable for enough treated school-years")
fwrite(bound, file.path(DIR_TABLES, "06e_g6_certain_bound.csv"))

# School is "grade-6-aligned" if grade 6 is certain in its FIRST treated wave
first_tr <- d[main_treatment == 1, .SD[which.min(wave)], by = gid][, .(gid, aligned = g6_certain)]
d <- merge(d, first_tr, by = "gid", all.x = TRUE)
ever_tr <- d[, .(ever = any(main_treatment == 1)), by = gid]
d <- merge(d, ever_tr, by = "gid")

# (3) Headline on the grade-6-CERTAIN subset (aligned treated + never-treated)
d_al <- d[(ever & aligned == 1) | (!ever)]
cdh <- function(dat, y) tryCatch(did_multiplegt_dyn(df = as.data.frame(dat), outcome = y, group = "gid",
            time = "wave", treatment = "main_treatment", effects = 2, placebo = 1, cluster = "gid",
            graph_off = TRUE), error = function(e) NULL)
ate <- function(r) if (is.null(r)) c(NA, NA) else r$results$ATE[1:2]
res <- rbindlist(list(
  data.table(sample = "full",              subject = "Language",    att = ate(cdh(d,   "y_lang"))[1]),
  data.table(sample = "full",              subject = "Mathematics", att = ate(cdh(d,   "y_math"))[1]),
  data.table(sample = "grade-6-certain",   subject = "Language",    att = ate(cdh(d_al,"y_lang"))[1]),
  data.table(sample = "grade-6-certain",   subject = "Mathematics", att = ate(cdh(d_al,"y_math"))[1])
))
res[, att := round(att, 4)]
fwrite(res, file.path(DIR_TABLES, "06e_grade_alignment.csv"))

# (4) By primary structure (6- vs 7-year). NOT REPORTED, and the reason was
# measured on 11/08 rather than assumed: this split is a PROVINCIAL split.
# 18 of the 24 jurisdictions are more than 90% homogeneous in structure and hold
# 89.9% of schools, so the two groups are two blocs of provinces. They differ
# sharply in the thing this dissertation already identifies as provincially
# structured -- the five-hour margin the census records as extended day: the
# weighted recording rate is 4.3% in the 6-year bloc and 47.6% in the 7-year one.
# The treated group of the 7-year bloc is therefore far more contaminated with
# five-hour schools, which predicts its smaller estimate on its own, with no
# grade-6 alignment story. The contrast cannot separate alignment from
# contamination and is not interpretable as evidence about either.
#
# The earlier comment here claimed the estimates were "similar" and reassuring.
# They are not similar (+0.046 against +0.011 in mathematics) and the difference
# is not about structure. The institutional fact that matters -- that each
# jurisdiction chooses a 6- or 7-year primary, which is why treatment is measured
# over grades 1 to 7 -- is stated in the manuscript's context section from the
# law itself (Ley 26.206, art. 134), not from this table.
#
# Kept because it is cheap and its absence would be a silent omission; the output
# carries no standard errors and must not be read as an estimate.
d <- merge(d, sch_g7, by = "school_id", suffixes = c("", "_sch"))
by_struct <- rbindlist(lapply(c(0L, 1L), function(g7) {
  ds <- d[has_g7_sch == g7]
  data.table(structure = ifelse(g7 == 1, "7-year (g6 second-to-last)", "6-year (g6 last)"),
             schools = uniqueN(ds$gid),
             cdh_lang = round(ate(cdh(ds, "y_lang"))[1], 4),
             cdh_math = round(ate(cdh(ds, "y_math"))[1], 4))
}), fill = TRUE)
fwrite(by_struct, file.path(DIR_TABLES, "06e_by_primary_structure.csv"))
