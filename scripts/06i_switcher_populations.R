# ==============================================================================
# 06i_switcher_populations.R — Who supports each horizon, and how they differ
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds
# Outputs : output/tables/06i_switcher_populations.csv   baseline profile by group
#           output/tables/06i_switcher_counts.csv        the nested sets
#           output/tables/06i_switcher_timing.csv        adoption wave by group
#           output/tables/06i_size_within_cohort.csv     size gap within cohort
#           output/tables/06i_size_lpm.csv               that gap, with cohort FE
#
# The aggregation window changes the answer, and `same_switchers` shows why: the
# switchers that support each horizon are different schools. Fixing the set at the
# two horizons the headline aggregates RAISES the ATT (+0.034 -> +0.069), while
# fixing it at three collapses it to zero on 337 schools. That is compositional
# heterogeneity, not mechanical bias, so the question this script answers is what
# distinguishes the three nested populations before treatment.
#
# The sets are reconstructed from the treatment path rather than read off the
# estimator, which does not expose them: F_g is the first wave whose treatment
# differs from the school's first observed wave, and a school supports horizon l
# if it is observed at F_g - 1 (the step) and at F_g - 1 + l. The reconstruction
# is validated against the switcher counts the estimator reports; it stops if they
# disagree, because everything downstream would then describe the wrong schools.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))

d <- as.data.table(readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds")))
setorder(d, school_id, wave)

# ---- The nested switcher sets -------------------------------------------------
first_switch <- d[, {
  ch <- which(main_treatment != main_treatment[1L])
  .(F = if (length(ch)) wave[ch[1L]] else NA_integer_)
}, by = school_id]

sets <- d[, .(waves = list(wave)), by = school_id]
sets <- merge(sets, first_switch, by = "school_id")
supports <- function(w, f, l) !is.na(f) && (f - 1L + l) %in% w
sets[, has_step := mapply(function(w, f) supports(w, f, 0L), waves, F)]
for (l in 1:3)
  sets[[paste0("k", l)]] <- mapply(function(w, f) supports(w, f, l), sets$waves, sets$F)

sw <- sets[has_step == TRUE & k1 == TRUE]
sw[, group := fifelse(!k2, "k=1 only",
              fifelse(!k3, "k=1 and 2", "k=1 to 3"))]

# Validation against the counts the estimator reports at each horizon, read from
# 06g rather than typed in, so the check survives a re-estimation and still catches
# a reconstruction that describes different schools. Compared on language, whose
# outcome has no missing values (mathematics has 306, so its counts run lower). The
# nesting is asserted separately because it holds by construction and a violation
# would mean the horizon logic itself broke.
if (!(nrow(sw) >= sum(sw$group != "k=1 only") &&
      sum(sw$group != "k=1 only") >= sum(sw$group == "k=1 to 3")))
  stop("switcher sets are no longer nested across horizons")

hzf <- file.path(DIR_TABLES, "06g_horizon_sensitivity.csv")
if (file.exists(hzf)) {
  want <- fread(hzf)[regime == "combined" & spec == "1 horizon(s)" &
                     subject == "language" & parameter == "Effect_1", switchers]
  if (length(want) == 1L && !identical(as.integer(nrow(sw)), as.integer(want)))
    stop("switcher reconstruction mismatch at the first horizon: got ", nrow(sw),
         ", the estimator reports ", want)
}

counts <- sw[, .(schools = .N), by = group]
counts[, pct := round(100 * schools / sum(schools), 1)]
fwrite(counts[order(group)], file.path(DIR_TABLES, "06i_switcher_counts.csv"))

# ---- Why a school supports only one horizon -----------------------------------
# The covariate profile below is a gradient, but the binding constraint is timing:
# a school that first switches in the last wave cannot have a second horizon, and
# one that switches in 2021 needs to appear in the 2022 SAMPLE wave to have one.
# Reported with the 2023-cohort share because the natural guess -- that these are
# the diagnostically awkward 2023 adopters -- turns out to be wrong.
attrs <- as.data.table(readRDS(file.path(DIR_PROCESSED, "school_attributes.rds")))
sw <- merge(sw, attrs[, .(school_id, cohort_2023)], by = "school_id", all.x = TRUE)
timing <- sw[, .(schools = .N,
                 cohort_2023 = sum(cohort_2023, na.rm = TRUE)), by = .(group, first_switch_wave = F)]
timing[, switch_year := WAVES[first_switch_wave]]
setorder(timing, group, first_switch_wave)
fwrite(timing, file.path(DIR_TABLES, "06i_switcher_timing.csv"))

# ---- Baseline profile of each set ---------------------------------------------
# One row per school: every variable here is baseline or time-invariant, so the
# comparison cannot be contaminated by post-treatment behaviour. `F` is the wave
# of first switch, which is a treatment characteristic rather than a covariate and
# is reported separately for that reason.
base <- unique(d[, .(school_id, sector, area, enroll_primary, enroll_g6, log_enroll,
                     bl_pct_female, bl_pct_repeater, bl_pct_preschool, bl_ses,
                     bl_pct_mother_sec, bl_score_lang, bl_score_math)], by = "school_id")
prof <- merge(sw[, .(school_id, group, first_switch_wave = F)], base, by = "school_id")
prof[, `:=`(is_public = as.integer(sector == "public"),
            is_rural  = as.integer(area == "rural"))]

VARS <- c("enroll_primary", "enroll_g6", "log_enroll", "is_public", "is_rural",
          "bl_ses", "bl_pct_mother_sec", "bl_pct_repeater", "bl_pct_preschool",
          "bl_pct_female", "bl_score_lang", "bl_score_math", "first_switch_wave")

prof[, (VARS) := lapply(.SD, as.numeric), .SDcols = VARS]
long <- melt(prof, id.vars = "group", measure.vars = VARS,
             variable.name = "variable", value.name = "v")
out <- long[!is.na(v), .(n = .N, mean = mean(v), sd = sd(v)), by = .(variable, group)]

# Standardised difference against the schools that support only the first horizon,
# pooling the two standard deviations. Reported instead of a p-value: with 2,297
# schools everything is significant, and the question is how big the gap is.
ref <- out[group == "k=1 only", .(variable, m0 = mean, s0 = sd)]
out <- merge(out, ref, by = "variable")
out[, std_diff := fifelse(group == "k=1 only", NA_real_,
                          (mean - m0) / sqrt((sd^2 + s0^2) / 2))]
out[, c("m0", "s0") := NULL]
setorder(out, variable, group)
fwrite(out, file.path(DIR_TABLES, "06i_switcher_populations.csv"))

# ---- Is the size gradient timing, or selection on top of timing? --------------
# The gradient above compares groups that differ in adoption calendar, and cohorts
# differ in size, so the two are confounded. The question a reader will ask is
# whether schools of the same adoption cohort still differ in size according to
# whether they support the second horizon. Reported descriptively -- within-cohort
# means and a linear probability model of support on log enrolment with cohort
# fixed effects -- because this diagnoses composition, it does not estimate an
# effect.
prof[, supports_k2 := as.integer(group != "k=1 only")]

by_cohort <- prof[, .(schools = .N,
                      supports_k2 = sum(supports_k2),
                      mean_log_enroll_yes = mean(log_enroll[supports_k2 == 1], na.rm = TRUE),
                      mean_log_enroll_no  = mean(log_enroll[supports_k2 == 0], na.rm = TRUE)),
                  by = first_switch_wave]
by_cohort[, `:=`(switch_year = WAVES[first_switch_wave],
                 gap = mean_log_enroll_yes - mean_log_enroll_no)]
setorder(by_cohort, first_switch_wave)
fwrite(by_cohort, file.path(DIR_TABLES, "06i_size_within_cohort.csv"))

# Pooled: the log-enrolment coefficient without and with cohort fixed effects. If
# it collapses once cohort is absorbed, the gradient was cohort composition.
m_raw <- lm(supports_k2 ~ log_enroll, data = prof)
m_fe  <- lm(supports_k2 ~ log_enroll + factor(first_switch_wave), data = prof)
lpm <- data.table(
  model = c("log enrolment only", "log enrolment + adoption-cohort fixed effects"),
  coef  = c(coef(m_raw)[["log_enroll"]], coef(m_fe)[["log_enroll"]]),
  se    = c(summary(m_raw)$coefficients["log_enroll", "Std. Error"],
            summary(m_fe)$coefficients["log_enroll", "Std. Error"]),
  n     = c(nobs(m_raw), nobs(m_fe)))
lpm[, t := coef / se]
fwrite(lpm, file.path(DIR_TABLES, "06i_size_lpm.csv"))
