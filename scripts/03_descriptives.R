# ==============================================================================
# 03_descriptives.R — Descriptive statistics, distributions, and missingness
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/panel.rds, ra_panel.rds  (from 01_build_panel.R)
# Outputs : output/tables/03_*.csv   descriptive tables
#           output/figures/03_*.png  distribution plots
#
# Contents: (1) sample sizes and match rates; (2) missingness by variable and
#           wave; (3) summary statistics of key variables; (4) distributions
#           (treatment intensity, test scores, enrolment); (5) treatment
#           descriptives for the identification strategy: intensity,
#           threshold sensitivity, extended vs full-day, adoption cohorts and
#           reversions, provincial
#           timing, assessment coverage.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))

panel    <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
ra_panel <- readRDS(file.path(DIR_PROCESSED, "ra_panel.rds"))

# Descriptive figures are drawn here (they plot the raw panel, not saved
# estimates), but they use the shared style AND the shared save routine, so they
# match the results figures produced by 91_figures.R.
#
# This script used to redefine save_fig() right after sourcing 90_style.R, which
# shadowed the module's copy: a different signature (w/h against width/height), a
# different default height, and the .png appended here instead of by the caller.
# The heights below are therefore passed explicitly -- they were the shadow
# function's defaults, and writing them down keeps the figures at the size they
# have always had rather than letting them inherit the module's 4.6.
source(here::here("scripts", "90_style.R"))
theme_set(theme_thesis())

# ============================================================================
# 1. Sample overview
# ============================================================================

overview <- panel[, .(
  schools        = .N,
  students_tested= sum(n_tested, na.rm = TRUE),
  matched_ra     = sum(!is.na(share_expanded)),
  match_rate_pct = round(100 * mean(!is.na(share_expanded)), 2),
  mean_school_size = round(mean(n_tested, na.rm = TRUE), 1)
), by = .(year, wave_type)][order(year)]
fwrite(overview, file.path(DIR_TABLES, "03_sample_overview.csv"))

# ============================================================================
# 2. Missingness by variable and wave
# ============================================================================

key_vars <- intersect(
  c("score_lang", "score_math", "pct_below_basic_lang", "pct_satisf_plus_lang",
    "share_expanded", "enroll_primary", "province", "sector", "area",
    "pct_female", "pct_preschool", "pct_repeater", "ses_index",
    "pct_mother_secondary", "climate_index", "overage_g6", "repeaters_g6",
    "n_meals", "nonpromotion_rate_g6", "dropout_rate_g6"),
  names(panel))

missing_tab <- panel[, c(list(n = .N),
                         lapply(.SD, function(x) round(100 * mean(is.na(x)), 1))),
                     by = year, .SDcols = key_vars][order(year)]
fwrite(missing_tab, file.path(DIR_TABLES, "03_missingness_by_wave.csv"))

# ============================================================================
# 3. Summary statistics of key variables (pooled and by wave)
# ============================================================================

num_vars <- intersect(
  c("score_lang", "score_math", "share_expanded", "enroll_primary", "n_tested",
    "pct_below_basic_lang", "pct_satisf_plus_lang", "pct_female",
    "pct_preschool", "pct_repeater", "overage_g6", "n_meals"),
  names(panel))

summarise_var <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(as.list(rep(NA_real_, 8)))
  list(n = length(x), mean = mean(x), sd = sd(x), min = min(x),
       p25 = quantile(x, .25), median = median(x), p75 = quantile(x, .75),
       max = max(x))
}
summary_tab <- rbindlist(lapply(num_vars, function(v) {
  s <- summarise_var(panel[[v]])
  data.table(variable = v, n = s$n,
             mean = round(s$mean, 3), sd = round(s$sd, 3),
             min = round(s$min, 3), p25 = round(s$p25, 3),
             median = round(s$median, 3), p75 = round(s$p75, 3),
             max = round(s$max, 3))
}))
fwrite(summary_tab, file.path(DIR_TABLES, "03_summary_stats_pooled.csv"))

# Score means by wave: sanity check against the anchored scale (2016 = 500)
scores_by_wave <- panel[, .(
  score_lang = round(mean(score_lang, na.rm = TRUE), 1),
  score_math = round(mean(score_math, na.rm = TRUE), 1),
  sd_between_schools_lang = round(sd(score_lang, na.rm = TRUE), 1),
  pct_below_basic_lang = round(100 * mean(pct_below_basic_lang, na.rm = TRUE), 1)
), by = year][order(year)]
fwrite(scores_by_wave, file.path(DIR_TABLES, "03_scores_by_wave.csv"))

# ============================================================================
# 4. Distributions
# ============================================================================

# theme_thesis() blanks plot.title and plot.subtitle by design -- a figure's only
# title is its Quarto caption -- so descriptive figures that never enter the
# manuscript carry their description in `caption`, which the theme does render.
# They used to set `title`/`subtitle`, which drew nothing at all.

# Treatment intensity (informs the operational threshold, Decision 1)
p <- ggplot(ra_panel[!is.na(share_expanded) & share_expanded <= 1.5],
            aes(share_expanded)) +
  geom_histogram(bins = 50, fill = "grey55", colour = "grey30") +
  labs(caption = paste("Treatment intensity: share of primary enrolment in extended/full-day.",
                       "All RA school-years, 2011-2025. Bimodality favours a binary treatment."),
       x = "share_expanded", y = "school-years")
save_fig(p, "03_treatment_intensity.png", height = 5)

p <- ggplot(panel[!is.na(score_lang)], aes(score_lang)) +
  geom_density(aes(colour = factor(year)), linewidth = .7) +
  labs(caption = "School-mean language scores by wave (anchored TRI scale).",
       x = "score_lang (school mean)", colour = "wave")
save_fig(p, "03_score_lang_density.png", height = 5)

p <- ggplot(panel[!is.na(score_math)], aes(score_math)) +
  geom_density(aes(colour = factor(year)), linewidth = .7) +
  labs(caption = "School-mean mathematics scores by wave (anchored TRI scale).",
       x = "score_math (school mean)", colour = "wave")
save_fig(p, "03_score_math_density.png", height = 5)

p <- ggplot(panel[!is.na(enroll_primary) & enroll_primary > 0],
            aes(enroll_primary)) +
  geom_histogram(bins = 60, fill = "grey55", colour = "grey30") +
  scale_x_log10() +
  labs(caption = "Primary enrolment of Aprender schools (log scale).",
       x = "enroll_primary (log10)", y = "school-waves")
save_fig(p, "03_enrolment_dist.png", height = 5)

# ============================================================================
# 5. Treatment descriptives for the identification strategy
# ============================================================================

# --- 5.1 Threshold sensitivity ------------------------------------------------
# Two natural denominators, both reported:
#   * % of school-years (out of school-years with a valid denominator) = how
#     common treatment is in the data (the panel unit).
#   * % of unique schools ever treated (out of schools with any valid share) =
#     how many distinct schools the definition captures.
thresholds <- c(0.01, 0.10, 0.25, 0.50, 0.90)
den_sy      <- sum(!is.na(ra_panel$share_expanded))            # valid school-years
den_schools <- ra_panel[!is.na(share_expanded), uniqueN(school_id)]
thr_tab <- rbindlist(lapply(thresholds, function(u) {
  sy   <- sum(ra_panel$share_expanded >= u, na.rm = TRUE)
  esc  <- ra_panel[share_expanded >= u, uniqueN(school_id)]
  data.table(threshold = u,
             treated_school_years = sy,
             pct_of_school_years  = round(100 * sy / den_sy, 1),
             schools_ever_treated = esc,
             pct_of_schools       = round(100 * esc / den_schools, 1))
}))
fwrite(thr_tab, file.path(DIR_TABLES, "03_threshold_sensitivity.csv"))

# --- 5.2 Extended vs full-day composition --------------------------------------
comp <- ra_panel[n_expanded >= 1, .(
  extended_only = sum(n_extended_day >= 1 & n_full_day == 0),
  full_only     = sum(n_extended_day == 0 & n_full_day >= 1),
  both          = sum(n_extended_day >= 1 & n_full_day >= 1)
)]
fwrite(comp, file.path(DIR_TABLES, "03_extended_vs_full.csv"))

# --- 5.3 Adoption cohorts and reversions --------------------------------------
# These used to be built from a LEVEL indicator (share >= 0.50 in a single year),
# which is NOT the treatment any estimator uses and produced a cohort table that
# disagreed with the rest of the analysis: 3,809 schools first treated in 2023
# against the 2,188 the estimation actually works with. The definitions now come
# from `treatment_annual` (built in 02), so the descriptive table and the
# estimation describe the same object.
treat_ann <- readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))
est_ids   <- unique(readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))$school_id)
setorder(treat_ann, school_id, year)
ann <- treat_ann[!is.na(main_treatment)]

# First year of adoption on the primary indicator. Schools already treated in the first census
# year are left-censored: the RA starts in 2011, so their true adoption year is
# unobserved and they must not be read as a 2011 cohort.
first_yr <- ann[main_treatment == 1, .(first_treated = min(year)), by = school_id]
min_year <- ann[, min(year)]
cohorts  <- merge(data.table(school_id = unique(ann$school_id)), first_yr,
                  by = "school_id", all.x = TRUE)
cohorts[, cohort := fifelse(is.na(first_treated), "never treated",
                    fifelse(first_treated == min_year,
                            paste0("already treated in ", min_year, " (left-censored)"),
                            as.character(first_treated)))]
cohorts[, in_est := school_id %in% est_ids]
cohort_dist <- cohorts[, .(schools = .N, in_estimation_sample = sum(in_est)),
                       by = .(cohort, first_treated)][order(first_treated, na.last = TRUE)]
cohort_dist[, first_treated := NULL]
fwrite(cohort_dist, file.path(DIR_TABLES, "03_adoption_cohorts.csv"))

# Reversions by province, on the primary definition.
revs <- ann[, .(reverts = .N >= 2 && any(diff(main_treatment) < 0),
                ever    = any(main_treatment == 1)), by = school_id]
rev_by_prov <- merge(revs[ever == TRUE & reverts == TRUE],
                     unique(ra_panel[!is.na(province)],
                            by = "school_id")[, .(school_id, province)],
                     by = "school_id", all.x = TRUE)[, .N, by = province][order(-N)]
fwrite(rev_by_prov, file.path(DIR_TABLES, "03_reversions_by_province.csv"))

# --- 5.3b Baseline balance: ever-treated vs never-treated ----------------------
# What an examiner asks for before any estimate: were the schools that eventually
# adopt comparable at baseline? They are not, and the design does not require
# them to be -- difference-in-differences identifies from CHANGES, not levels --
# but the gaps have to be visible, because they are what the conditional
# parallel-trends assumption is asked to carry. Covariates are the 2016 baseline
# vector used by the doubly-robust specification; the last column is the
# standardised difference, which is scale-free and the conventional balance
# metric (a school-count difference in means would not be comparable across rows).
est  <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
bal_ids <- unique(est[, .(school_id, sector, area, log_enroll, bl_pct_female,
                          bl_pct_repeater, bl_pct_preschool, bl_pct_mother_sec,
                          bl_ses, bl_score_lang, bl_score_math)], by = "school_id")
ever <- ann[, .(ever_treated = any(main_treatment == 1)), by = school_id]
bal  <- merge(bal_ids, ever, by = "school_id")
bal[, `:=`(pct_public = as.integer(sector == "public"),
           pct_rural  = as.integer(area == "rural"))]

BAL_VARS <- c("pct_public", "pct_rural", "log_enroll", "bl_pct_female",
              "bl_pct_repeater", "bl_pct_preschool", "bl_pct_mother_sec",
              "bl_ses", "bl_score_lang", "bl_score_math")
balance <- rbindlist(lapply(BAL_VARS, function(v) {
  t <- bal[ever_treated == TRUE][[v]]; c0 <- bal[ever_treated == FALSE][[v]]
  mt <- mean(t, na.rm = TRUE); mc <- mean(c0, na.rm = TRUE)
  sd_pool <- sqrt((var(t, na.rm = TRUE) + var(c0, na.rm = TRUE)) / 2)
  data.table(variable = v,
             ever_treated = round(mt, 3), never_treated = round(mc, 3),
             difference = round(mt - mc, 3),
             std_difference = round((mt - mc) / sd_pool, 3),
             n_treated = sum(!is.na(t)), n_control = sum(!is.na(c0)))
}))
fwrite(balance, file.path(DIR_TABLES, "03_baseline_balance.csv"))

# --- 5.3c How the treatment MOVES ---------------------------------------------
# The identification argument rests on the treatment switching on and off, so the
# switching itself has to be shown rather than asserted: how many transitions go
# up, how many go down, and how many schools switch more than once. Computed on
# the estimation sample, between consecutive outcome waves, which is the variation
# the estimator actually uses.
est_sw <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
setorder(est_sw, gid, wave)
est_sw[, dD := main_treatment - shift(main_treatment), by = gid]
moves <- est_sw[!is.na(dD), .N, by = dD][order(dD)]
moves[, `:=`(direction = fifelse(dD > 0, "increase", fifelse(dD < 0, "decrease", "no change")),
             pct = round(100 * N / sum(N), 2))]
fwrite(moves[, .(delta_D = dD, direction, transitions = N, pct)],
       file.path(DIR_TABLES, "03_treatment_transitions.csv"))

n_switch <- est_sw[!is.na(dD), .(switches = sum(dD != 0)), by = gid][
  , .(schools = .N), by = switches][order(switches)]
n_switch[, pct := round(100 * schools / sum(schools), 1)]
fwrite(n_switch, file.path(DIR_TABLES, "03_switches_per_school.csv"))

# --- 5.4 Provincial adoption timing (share of treated schools) -----------------
prov_timing <- merge(ann, unique(ra_panel[!is.na(province)], by = "school_id")[
                       , .(school_id, province)], by = "school_id")[
                     , .(pct_treated = round(100 * mean(main_treatment), 1)), by = .(province, year)]
# The legend has to name the variable, not its units: the shading is the share of
# a province's primary schools whose expanded share reaches the 50% threshold in
# that year, which is neither a share of enrolment nor a share of the country.
# No internal title -- body figures carry their caption from the chunk, as every
# figure built in 91_figures.R does.
p <- ggplot(prov_timing, aes(year, province, fill = pct_treated)) +
  geom_tile(colour = "white") +
  scale_fill_gradient(low = "white", high = "steelblue4",
                      name = "Primary schools above the treatment threshold (%)") +
  labs(x = NULL, y = NULL) +
  # The project theme blanks every legend title, which is right for the discrete
  # legends elsewhere (their categories are self-naming) and wrong here: a
  # continuous colour scale that does not say what it encodes is unreadable.
  # Overridden locally so no other figure changes.
  theme(legend.title = element_text(colour = PAL$ink))
save_fig(p, "03_provincial_timing.png", width = 9, height = 7)
fwrite(prov_timing, file.path(DIR_TABLES, "03_provincial_timing.csv"))

# --- 5.5 Schools by regime and year (extended vs full-day) --------------------
# Counts and shares of primary schools under each regime, by RA year. The
# denominator each year is the primary schools observed that year (enroll_primary
# > 0). Three cut-offs make explicit how "whole-school" each regime is:
#   has_regime   : the school has any pupils in the regime (n > 0);
#   share_ge_50  : the regime covers at least half of primary enrolment;
#   share_eq_100 : the whole primary is in the regime (share >= 1; a handful of
#                  school-years exceed 1 by an RA count/enrolment mismatch and
#                  are included here as fully converted).
# Extended-day and full-day are reported in separate tables. Full-day is close
# to all-or-nothing (has_regime ~ share_ge_50 ~ share_eq_100), whereas extended
# day is often partial (has_regime > share_ge_50), which is why the treatment
# threshold matters for the combined/extended measure but not for full-day.
regime_by_year <- function(count_col) {
  d <- ra_panel[!is.na(enroll_primary) & enroll_primary > 0]
  d[, sh := get(count_col) / enroll_primary]
  out <- d[, .(`primary schools` = .N,
               has_regime   = sum(get(count_col) > 0, na.rm = TRUE),
               share_ge_50  = sum(sh >= 0.50, na.rm = TRUE),
               share_eq_100 = sum(sh >= 1,    na.rm = TRUE)), by = year][order(year)]
  out[, `:=`(has_regime_pct   = round(100 * has_regime   / `primary schools`, 1),
             share_ge_50_pct  = round(100 * share_ge_50  / `primary schools`, 1),
             share_eq_100_pct = round(100 * share_eq_100 / `primary schools`, 1))]
  setcolorder(out, c("year", "primary schools",
                     "has_regime", "has_regime_pct",
                     "share_ge_50", "share_ge_50_pct",
                     "share_eq_100", "share_eq_100_pct"))
  out[]
}
fwrite(regime_by_year("n_extended_day"),
       file.path(DIR_TABLES, "03_extended_day_by_year.csv"))
fwrite(regime_by_year("n_full_day"),
       file.path(DIR_TABLES, "03_full_day_by_year.csv"))

# --- 5.6 Treated schools by outcome wave --------------------------------------
# How many schools are treated at each Aprender wave (the years that have an
# outcome), from the start of the panel (2016) to the end (2025). Built from the
# estimation sample (02_build_estimation_sample.R), which carries the primary
# treatment definitions; one row per school x wave, so the row count per wave is
# the number of schools tested that wave.
#
# Columns (the binary is contemporaneous: the share clears the cut-off in the year itself):
#   year                          Aprender wave
#   wave_type                     census / sample (2022 is a school-level sample,
#                                 so its LEVELS are not comparable to census waves;
#                                 compare the percentages instead)
#   schools_tested                schools with a grade-6 outcome that wave
#                                 (the denominator for the percentages below)
#   treated_ext_or_full_ge50      # treated under the HEADLINE definition: extended-
#                                 OR full-day at >=50% of primary enrolment (main_treatment)
#   pct_treated_ext_or_full_ge50  the same as % of schools_tested
#   treated_fullday_ge50          # treated under full-day only at >=50% (main_treatment_full)
#   pct_treated_fullday_ge50      the same as % of schools_tested
#   treated_any_expansion         # any expansion at all, share >0
#                                 the Nistal-style benchmark)
#   pct_treated_any_expansion     the same as % of schools_tested
est <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
treated_by_wave <- est[, .(
  schools_tested           = .N,
  treated_ext_or_full_ge50 = sum(main_treatment == 1, na.rm = TRUE),
  treated_fullday_ge50     = sum(main_treatment_full == 1, na.rm = TRUE),
  treated_any_expansion    = sum(share_expanded > 0, na.rm = TRUE)
), by = year][order(year)]
treated_by_wave[, wave_type := WAVE_TYPE[as.character(year)]]
treated_by_wave[, `:=`(
  pct_treated_ext_or_full_ge50 = round(100 * treated_ext_or_full_ge50 / schools_tested, 1),
  pct_treated_fullday_ge50     = round(100 * treated_fullday_ge50     / schools_tested, 1),
  pct_treated_any_expansion    = round(100 * treated_any_expansion    / schools_tested, 1))]
setcolorder(treated_by_wave, c(
  "year", "wave_type", "schools_tested",
  "treated_ext_or_full_ge50", "pct_treated_ext_or_full_ge50",
  "treated_fullday_ge50",     "pct_treated_fullday_ge50",
  "treated_any_expansion",    "pct_treated_any_expansion"))
fwrite(treated_by_wave, file.path(DIR_TABLES, "03_treated_by_wave.csv"))

# --- 5.7 Assessment coverage: who sits Aprender, and how they differ ----------
# Linkage (how many assessed schools match the census) is near-complete, but
# coverage (how many primary schools are assessed at all) is not: Aprender
# reaches roughly two thirds of them, and the ones it misses are much smaller.
# The manuscript reads these figures inline, so they must not be hard-coded.
census_waves <- names(WAVE_TYPE)[WAVE_TYPE == "census"]
cw <- panel[year %in% as.integer(census_waves)]
cw[, assessed := !is.na(score_lang)]
coverage <- rbind(
  data.table(
    metric = "primary schools in the RA census (ever)",
    value  = panel[, uniqueN(school_id)]),
  data.table(
    metric = "schools with at least one Aprender score",
    value  = panel[!is.na(score_lang), uniqueN(school_id)]),
  data.table(
    metric = "school-years in census waves",
    value  = nrow(cw)),
  data.table(
    metric = "school-years in census waves that are assessed",
    value  = sum(cw$assessed)),
  cw[, .(metric = paste0("median primary enrolment, ",
                         fifelse(assessed, "assessed", "not assessed")),
         value  = median(enroll_primary, na.rm = TRUE)), by = assessed][, -"assessed"],
  cw[, .(metric = paste0("% treated (share >= 50%), ",
                         fifelse(assessed, "assessed", "not assessed")),
         value  = round(100 * mean(share_expanded >= 0.5, na.rm = TRUE), 1)),
     by = assessed][, -"assessed"])
coverage[, value := round(value, 1)]
# Written once, below, with the panel-balance rows appended. It used to be written
# here as well, and that first copy was overwritten a few lines later.

# Panel balance: a school enters in any wave it is assessed, so the panel is
# unbalanced. The manuscript states the share present in all six waves.
# Computed on the ESTIMATION PANEL, not on every assessed school: the sentence
# it feeds describes how unbalanced the panel the estimates run on actually is,
# and the two sets differ (22,935 against 22,942) because a handful of assessed
# schools carry no treatment history.
bal <- est[, .N, by = school_id][
  , .(metric = "% of estimation-panel schools present in all six waves",
      value  = round(100 * mean(N == length(WAVE_TYPE)), 1))]
# The six-wave figure reads as a fragile panel and is not one: one of the six is
# a SAMPLE wave, so most schools are absent from it by design rather than by
# attrition. The five-census-wave figure is the one that describes how stable
# the panel actually is, and the gap between them is the 2022 design.
cens_yrs <- as.integer(names(WAVE_TYPE)[WAVE_TYPE == "census"])
bal5 <- est[year %in% cens_yrs, .N, by = school_id][
  , .(metric = "% of estimation-panel schools present in all five census waves",
      value  = round(100 * mean(N == length(cens_yrs)), 1))]
bal <- rbind(bal, bal5)
fwrite(rbind(coverage, bal), file.path(DIR_TABLES, "03_assessment_coverage.csv"))


# ---- 6. Jurisdiction profile ------------------------------------------------
# One row per jurisdiction, describing the system the programme landed on rather
# than the programme itself: how many primary schools it has, how they split by
# management sector and by area, how big they are, and where the treated share
# stood at each end of the census. The point of the table is that the 24
# jurisdictions are not comparable units -- they differ by two orders of
# magnitude in size and by their whole composition -- so a national average over
# them is a weighted statement about very unequal systems, and the provincial
# clustering question of 06f is a question about these cells.
#
# Structural columns are taken from the LAST census year, so the table describes
# the system as it now stands; the treated shares are given at both ends so the
# movement is visible. Schools are counted once, from their last appearance.
LAST <- max(panel$year)

prof_struct <- panel[year == LAST, .(
  schools      = .N,
  pct_state    = round(100 * mean(sector == "public", na.rm = TRUE), 1),
  pct_rural    = round(100 * mean(area == "rural",    na.rm = TRUE), 1),
  enrol        = sum(enroll_primary, na.rm = TRUE),
  median_size  = as.numeric(median(enroll_primary, na.rm = TRUE))), by = province]

# Treated share at each end, on the same indicator the estimates use.
tr <- function(y) panel[year == y & !is.na(share_expanded),
                        .(v = round(100 * mean(share_expanded >= 0.50), 1)), by = province]
prof <- Reduce(function(a, b) merge(a, b, by = "province", all.x = TRUE),
               list(prof_struct,
                    setnames(tr(min(panel$year)), "v", "pct_treated_first"),
                    setnames(tr(LAST),            "v", "pct_treated_last")))

# How much of each jurisdiction the assessment frame reaches, which is what makes
# the estimand local (02, and the coverage note in the Data section). Numerator
# and denominator are both EVER-counts: mixing a last-year denominator with an
# ever-assessed numerator returns coverage above 100% wherever schools closed
# between their assessment and the last census year.
ever <- panel[, .(schools_ever = uniqueN(school_id)), by = province]
assessed <- panel[!is.na(score_lang), .(schools_assessed = uniqueN(school_id)), by = province]
prof <- merge(merge(prof, ever, by = "province", all.x = TRUE),
              assessed, by = "province", all.x = TRUE)
prof[, pct_assessed := round(100 * schools_assessed / schools_ever, 1)]
setnames(prof, "schools", "schools_last")

setorder(prof, -schools_last)
fwrite(prof, file.path(DIR_TABLES, "03_jurisdiction_profile.csv"))
