# ==============================================================================
# 06l_exposure_validity.R — Does "years since first switch" equal "years treated"?
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds
#           data/processed/treatment_annual.rds
# Outputs : output/tables/06l_exposure_validity.csv
#
# The exposure classes in 06j are built from time elapsed since the first switch.
# In an ABSORBING design that is the same thing as time treated. This design is
# not absorbing -- 32% of treated schools revert and 17.7% re-adopt -- so the two
# can separate: a school can switch on, switch off, and be measured at a wave
# having accumulated fewer treated years than have elapsed.
#
# This checks how often they separate, by exposure class. If they coincide for
# nearly every switcher, the partition does what it claims. If they do not, the
# classes have to be built on years actually treated instead.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))

d  <- as.data.table(readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds")))
ta <- as.data.table(readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds")))
setorder(d, school_id, wave)

fs <- d[, {
  ch <- which(main_treatment != main_treatment[1L])
  .(F = if (length(ch)) wave[ch[1L]] else NA_integer_,
    d1 = main_treatment[1L], y0 = WAVES[wave[1L]])
}, by = school_id][!is.na(F)]

yr <- merge(ta[, .(school_id, year, tr = main_treatment)], fs[, .(school_id, d1, y0)],
            by = "school_id")
chg <- yr[year >= y0][, {
  ch <- which(tr != d1[1L])
  .(switch_year = if (length(ch)) year[ch[1L]] else NA_integer_)
}, by = school_id]

sw <- merge(fs, chg, by = "school_id")[!is.na(switch_year)]
sw[, years_since_adoption := WAVES[F] - switch_year]

# Years actually treated between the switch year and the wave at each horizon,
# counted on the annual series. `elapsed` is what the partition uses.
obs <- d[, .(waves = list(wave)), by = school_id]
sw  <- merge(sw, obs, by = "school_id")
ann <- ta[, .(school_id, year, tr = main_treatment)]

rows <- rbindlist(lapply(1:2, function(l) {
  ok <- mapply(function(w, f) (f - 1L + l) %in% w, sw$waves, sw$F)
  x  <- sw[ok]
  x[, target_year := WAVES[F - 1L + l]]
  z  <- merge(ann, x[, .(school_id, switch_year, target_year, years_since_adoption)],
              by = "school_id")
  z  <- z[year >= switch_year & year <= target_year]
  agg <- z[, .(treated_years = sum(tr == 1L, na.rm = TRUE),
               window_years  = .N), by = .(school_id, years_since_adoption)]
  agg[, `:=`(horizon = l, gap = window_years - treated_years)]
  agg[]
}), fill = TRUE)

out <- rows[, .(switchers = .N,
                mean_window = round(mean(window_years), 2),
                mean_treated = round(mean(treated_years), 2),
                pct_no_gap = round(100 * mean(gap == 0), 1),
                pct_gap_ge2 = round(100 * mean(gap >= 2), 1)),
            by = .(horizon, years_since_adoption)]
setorder(out, horizon, years_since_adoption)
fwrite(out, file.path(DIR_TABLES, "06l_exposure_validity.csv"))
