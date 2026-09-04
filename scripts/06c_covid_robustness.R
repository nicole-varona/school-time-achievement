# ==============================================================================
# 06c_covid_robustness.R — Which waves enter: COVID waves and the sample wave
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# The 2021 Aprender wave is anomalous (post-COVID score trough), and the
# decomposition (07b) shows the Maths k=-2 placebo is driven by the 2023 cohort's
# 2021->2023 comparison, which straddles it. This checks how the score estimates
# and the Maths pre-trend respond to excluding the COVID-affected wave(s):
#   (a) full 6 waves (headline);  (b) drop 2021 (COVID trough);
#   (c) drop 2021 + 2022 (COVID-disrupted period);
#   (d) drop 2022 alone — the only SAMPLE wave, added 06/08: it tests external
#       validity (does the answer depend on including a wave that is not a
#       census?) rather than the pandemic, and it removes the fragile G=2022
#       cohort, which holds ~35 schools because only those drawn into the
#       sample can belong to it.
# The treatment definition is unchanged; we only drop outcome observations and
# re-index waves / recompute cohorts. NOTE that dropping a wave RE-DATES
# adoption for schools first observed treated in it, so these rows change the
# estimand and are diagnostics rather than preferred specifications.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(did); library(DIDmultiplegtDYN); library(polars)})
# Randomness: seeded per call via set_call_seed() (00_setup.R), not globally, so
# adding or reordering a spec cannot move another spec's bootstrapped SEs.

panel <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
# Treatment comes from 02 (the single source of truth) instead of being rebuilt
# here. The persistence rule used to live in nine copies, and shift() silently
# depends on the panel being sorted -- an invariant no script declared.
panel <- merge(panel,
               readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))[
                 , .(school_id, year, main_treatment)],
               by = c("school_id", "year"), all.x = TRUE)
base <- panel[year %in% WAVES & !is.na(score_lang),
              .(school_id, year, y_lang = score_lang / 100, y_math = score_math / 100, main_treatment)]
base <- base[!is.na(main_treatment)]

ate <- function(r) if (is.null(r)) c(NA_real_, NA_real_) else r$results$ATE[1:2]
placebo_k <- function(r, kk) {
  if (is.null(r) || is.null(r$results$Placebos)) return(c(NA_real_, NA_real_))
  p <- as.data.table(r$results$Placebos)
  p[, k := -as.integer(sub("Placebo_", "", rownames(r$results$Placebos)))]
  row <- p[k == kk]; if (!nrow(row)) return(c(NA_real_, NA_real_))
  c(row$Estimate, row$SE)
}

run_spec <- function(drop_years, drop_schools, label) {
  set_call_seed(label)
  ds <- copy(base[!(year %in% drop_years) & !(school_id %in% drop_schools)])
  wv <- sort(unique(ds$year))
  ds[, `:=`(wave = match(year, wv), gid = as.integer(factor(school_id)))]
  cdh <- function(y) tryCatch(did_multiplegt_dyn(df = as.data.frame(ds), outcome = y,
              group = "gid", time = "wave", treatment = "main_treatment", effects = 2,
              placebo = 2, cluster = "gid", graph_off = TRUE), error = function(e) NULL)
  cm <- cdh("y_math"); cl <- cdh("y_lang")
  dcs <- cs_sample(ds)
  cs <- function(y) { fit <- tryCatch(att_gt(yname = y, tname = "wave", idname = "gid",
            gname = "G", data = dcs, control_group = "notyettreated",
            allow_unbalanced_panel = TRUE, bstrap = TRUE,
            biters = BITERS, est_method = "dr"), error = function(e) NULL)
    if (is.null(fit)) return(c(NA_real_, NA_real_))
    a <- suppressWarnings(aggte(fit, type = "simple", na.rm = TRUE)); c(a$overall.att, a$overall.se) }
  csm <- cs("y_math"); csl <- cs("y_lang")
  k2m <- placebo_k(cm, -2); k1m <- placebo_k(cm, -1)
  data.table(spec = label, waves = paste(wv, collapse = ","),
             cdh_lang = ate(cl)[1], cdh_math = ate(cm)[1],
             cs_lang = csl[1], cs_lang_se = csl[2], cs_math = csm[1], cs_math_se = csm[2],
             math_k2 = k2m[1], math_k2_se = k2m[2], math_k2_sig = as.integer(abs(k2m[1]) > 1.96 * k2m[2]),
             math_k1 = k1m[1], math_k1_se = k1m[2])
}

# Cohort first observed treated at the 2023 outcome wave (the
# adoption cohort where the Maths placebo is concentrated; 07b).
coh2023 <- base[main_treatment == 1, .(g = min(year)), by = school_id][g == 2023, school_id]

res <- rbindlist(list(
  run_spec(integer(0),      character(0), "full (6 waves, headline)"),
  run_spec(2021L,           character(0), "drop 2021 (COVID trough)"),
  run_spec(c(2021L, 2022L), character(0), "drop 2021+2022 (COVID period)"),
  run_spec(2022L,           character(0), "drop 2022 (sample wave)"),
  # Advised diagnostic: remove the cohort GENERATING the concern, keep all waves.
  # Present as a diagnostic, NOT the preferred spec.
  run_spec(integer(0),      coh2023,      "exclude 2023 cohort (full waves)")
), fill = TRUE)
res_r <- copy(res); num <- setdiff(names(res_r), c("spec", "waves"))
res_r[, (num) := lapply(.SD, function(x) round(x, 4)), .SDcols = num]
fwrite(res_r, file.path(DIR_TABLES, "06c_covid_robustness.csv"))
