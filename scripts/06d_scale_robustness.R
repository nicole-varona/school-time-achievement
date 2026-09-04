# ==============================================================================
# 06d_scale_robustness.R — Is the outcome scale choice consequential?
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Headline outcome = anchored Aprender score / 100 (SD units of the 2016 anchored
# TRI scale, comparable across waves by design). The comparability report flags a
# scaling issue in the 2021 wave — the same wave driving the COVID/pre-trend
# concern. This checks the alternative: standardising WITHIN each wave (z-score of
# school means by wave), which is robust to any wave-to-wave scale/dispersion
# drift but changes the estimand to a within-wave relative position.
#   (1) descriptive: dispersion of school means by wave (does 2021 look off?);
#   (2) re-estimate the CDH headline + Maths k=-2 placebo under both scalings.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(DIDmultiplegtDYN); library(polars)})
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
d <- panel[year %in% WAVES & !is.na(score_lang),
           .(school_id, year, score_lang, score_math, main_treatment)]
d <- d[!is.na(main_treatment)]

# (1) Dispersion of SCHOOL MEANS by wave (raw TRI points) — does 2021 look off?

# Outcomes: anchored (/100) vs within-wave z-score (school-mean SD per wave)
d[, `:=`(ya_lang = score_lang / 100, ya_math = score_math / 100)]
d[, `:=`(yz_lang = (score_lang - mean(score_lang, na.rm = TRUE)) / sd(score_lang, na.rm = TRUE),
         yz_math = (score_math - mean(score_math, na.rm = TRUE)) / sd(score_math, na.rm = TRUE)), by = year]
d[, `:=`(wave = match(year, WAVES), gid = as.integer(factor(school_id)))]

cdh <- function(y) tryCatch(did_multiplegt_dyn(df = as.data.frame(d), outcome = y, group = "gid",
            time = "wave", treatment = "main_treatment", effects = 2, placebo = 2, cluster = "gid",
            graph_off = TRUE), error = function(e) NULL)
ate <- function(r) if (is.null(r)) c(NA, NA) else r$results$ATE[1:2]
k2  <- function(r) { if (is.null(r) || is.null(r$results$Placebos)) return(c(NA, NA))
  p <- as.data.table(r$results$Placebos); p[, k := -as.integer(sub("Placebo_", "", rownames(r$results$Placebos)))]
  row <- p[k == -2]; if (!nrow(row)) return(c(NA, NA)); c(row$Estimate, row$SE) }

specs <- rbindlist(lapply(c("ya_lang", "ya_math", "yz_lang", "yz_math"), function(y) {
  r <- cdh(y); a <- ate(r); p <- k2(r)
  data.table(outcome = y, scale = ifelse(grepl("^ya", y), "anchored /100", "within-wave z"),
             subject = ifelse(grepl("lang", y), "Language", "Mathematics"),
             att = round(a[1], 4), se = round(a[2], 4), k2 = round(p[1], 4),
             k2_sig = as.integer(abs(p[1]) > 1.96 * p[2]))
}), fill = TRUE)
fwrite(specs, file.path(DIR_TABLES, "06d_scale_robustness.csv"))
