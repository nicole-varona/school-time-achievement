# ==============================================================================
# 06m_calendar_support.R — What population each candidate outcome calendar
#                          identifies, before any coefficient is looked at
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds
#           data/processed/treatment_annual.rds
#           data/processed/school_attributes.rds
# Outputs : data/processed/calendar_designs.rds        classification, for 06n
#           output/tables/06m_calendar_support.csv     support of each calendar
#           output/tables/06m_calendar_cohorts.csv     first-adoption cohorts
#           output/tables/06m_exposure_by_horizon.csv  calendar time per Effect_k
#           output/tables/06m_cohort_by_exposure.csv   cohorts within each class
#
# THE PROBLEM THIS ADDRESSES. The outcome is observed on an irregular calendar
# (2016, 2018, 2021, 2022, 2023, 2025: gaps of 2, 3, 1, 1, 2 years) while the
# treatment moves every year. `did_multiplegt_dyn` indexes event time by period,
# so a single `Effect_k` averages switchers whose calendar time since adoption
# differs by years.
# That damages the estimand (Effect_1 is not "one year after adoption") and the
# diagnostics alike: an aggregate placebo near zero can be cancellation between
# timing classes with pre-trends of opposite sign, which is what 06j found.
#
# THREE CANDIDATE CALENDARS are therefore compared as DESIGNS, not as robustness
# checks, because they identify different populations:
#   six waves      — every wave. Most information, most cohorts, most history;
#                    irregular spacing and one wave (2022) that is a sample.
#   five censuses  — drops 2022. Census only, four intervals of 2, 3, 2, 2 years.
#   three censuses — 2021, 2023, 2025. Regular two-year spacing, census only,
#                    the main treatment, but only three periods.
#
# This script measures WHO each calendar identifies. It reports no estimate, and
# the design is not to be chosen from coefficients: choosing after seeing them is
# specification search. The comparable estimates are produced by 06n, which reads
# the classification saved here rather than re-deriving it.
#
# TREATMENT: `main_treatment` throughout, the headline definition -- at least half
# of primary enrolment in extended or full day, dated on the annual census. It is
# used in three places and it is the SAME series
# in all three: the estimation in 06n, the switch year that dates adoption here,
# and `first_treated_year`, which dates the adoption cohorts. Nothing here varies
# the treatment; what varies is the calendar the outcome is observed on.
#
# TWO QUANTITIES THAT ARE NOT THE SAME THING, and the distinction is the point of
# this script. `years_since_adoption` is CALENDAR TIME between the annual switch
# and the wave at that horizon, which is what event time indexes. `treated_years`
# is ACCUMULATED EXPOSURE, the years actually treated inside that window. Under an
# absorbing treatment they coincide; here they do not, because a school can switch
# on, contract, and be measured having accumulated fewer treated years than have
# elapsed -- so 0,1,0,0 and 0,1,1,1 share a value of the first and differ in the
# second. Neither is called "years of exposure": the first is not exposure, and
# only the second is. Both come from the ANNUAL history, used ONLY to classify
# schools; the estimation in 06n stays at the outcome frequency.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))

DESIGNS <- list("six waves"      = c(2016, 2018, 2021, 2022, 2023, 2025),
                "five censuses"  = c(2016, 2018, 2021, 2023, 2025),
                "three censuses" = c(2021, 2023, 2025))
HORIZONS <- 1:3

d  <- as.data.table(readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds")))
ta <- as.data.table(readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds")))

# Both definitions come from 02, where they are built once and are NA on exactly
# the same school-years. Building the main treatment here from the share, as
# this script used to, reintroduced a difference in missing-data handling on top
# of the difference in rule, so the comparison confounded the two.
DEFS <- list(main = list(tv = "main_treatment",
                                    annual = ta[!is.na(main_treatment),
                                                .(school_id, year, tr = main_treatment)]))

if (!all(unlist(DESIGNS) %in% WAVES))
  stop("a design asks for a year with no outcome wave")
if (!all(WAVE_TYPE[as.character(unlist(DESIGNS[-1]))] == "census"))
  stop("the census-only designs must not contain the sample wave")

# ---- Classify one calendar ---------------------------------------------------
# Waves are re-indexed 1..K within the calendar, which is what the estimator sees.
# F is the first observed wave whose treatment differs from the school's first
# observed wave; the switch year is the first annual year, from the first observed
# wave onwards, at which the annual series differs. A school supports horizon l if
# it is observed at F-1 (the step) and at F-1+l.
build <- function(yrs, tv, ann) {
  x <- d[year %in% yrs & !is.na(get(tv))]
  x[, wv := match(year, yrs)]
  setorder(x, school_id, wv)

  s <- x[, {
    ch <- which(get(tv) != get(tv)[1L])
    .(n_waves = .N, waves = list(wv), ever = any(get(tv) == 1L),
      n_changes = length(ch), n_down = sum(diff(get(tv)) < 0),
      F = if (length(ch)) wv[ch[1L]] else NA_integer_,
      d1 = get(tv)[1L], y0 = yrs[wv[1L]])
  }, by = school_id]

  sw  <- s[!is.na(F)]
  yr  <- merge(ann, sw[, .(school_id, d1, y0)], by = "school_id")
  chg <- yr[year >= y0, {
    ch <- which(tr != d1[1L])
    .(switch_year = if (length(ch)) year[ch[1L]] else NA_integer_)
  }, by = school_id]
  sw <- merge(sw, chg, by = "school_id")[!is.na(switch_year)]

  hz <- rbindlist(lapply(HORIZONS, function(l) {
    ok <- mapply(function(w, f) (f - 1L) %in% w && (f - 1L + l) %in% w, sw$waves, sw$F)
    z  <- sw[ok, .(school_id, F, d1, switch_year)]
    z[, `:=`(horizon = l, target_year = yrs[F - 1L + l])]
    z[, years_since_adoption := target_year - switch_year]
    z[]
  }))

  # Years actually treated between the switch and the wave at that horizon.
  ty <- merge(ann, hz[, .(school_id, horizon, switch_year, target_year)],
              by = "school_id", allow.cartesian = TRUE)
  ty <- ty[year >= switch_year & year <= target_year,
           .(treated_years = sum(tr == 1L)), by = .(school_id, horizon)]
  hz <- merge(hz, ty, by = c("school_id", "horizon"), all.x = TRUE)
  hz[, treated_throughout := treated_years == years_since_adoption + 1L]

  list(schools = s, switchers = sw, horizons = hz)
}


# Every calendar under every definition. `each()` folds the two loops so that no
# downstream table can be built for one definition and not the other.
B <- lapply(DEFS, function(def)
       lapply(DESIGNS, function(yrs) build(yrs, def$tv, def$annual)))
each <- function(f) rbindlist(lapply(names(DEFS), function(dn)
          rbindlist(lapply(names(DESIGNS), function(nm) f(dn, nm)))))

# Adoption cohort and reversal are properties of the ANNUAL path, so they are
# re-derived per definition rather than read from school_attributes.rds, which
# only carries the older ones. Under the main treatment the 2023 cohort
# is larger and a 2025 cohort exists at all, which the retired persistence rule
# had removed by right-censoring.
ann_attr <- lapply(DEFS, function(def) def$annual[order(school_id, year), .(
  first_treated_year = if (any(tr == 1L)) min(year[tr == 1L]) else NA_integer_,
  reverts_annual     = any(diff(tr) < 0)), by = school_id])

# The six-wave classification is the one the rest of the project rests on, so it
# doubles as a check on this reconstruction: these are the switcher
# counts the estimator reports under each horizon (06g, language). They are
# support PER horizon, not the nested set: 06i counts 338 schools supporting all
# three at once, which is what `same_switchers` estimates on and a different thing.
k <- B$main[["six waves"]]$horizons[, .N, by = horizon][order(horizon), N]
if (length(k) != 3L || k[1] < k[2] || k[2] < k[3])
  stop("horizon support is no longer nested: got ", paste(k, collapse = "/"))

# ---- Support of each calendar ------------------------------------------------
support <- each(function(dn, nm) {
  b  <- B[[dn]][[nm]]
  s  <- merge(b$schools, ann_attr[[dn]], by = "school_id", all.x = TRUE)
  hz <- b$horizons
  lost <- setdiff(b$switchers$school_id, hz[horizon == 1L, school_id])
  data.table(
    definition          = dn,
    design              = nm,
    waves               = length(DESIGNS[[nm]]),
    schools             = nrow(s),
    single_wave         = s[n_waves == 1L, .N],
    ever_treated        = s[ever == TRUE, .N],
    never_treated       = s[ever == FALSE, .N],
    never_switchers     = s[is.na(F), .N],
    switchers           = s[!is.na(F), .N],
    switchers_dated     = nrow(b$switchers),
    switch_up           = b$switchers[d1 == 0L, .N],
    switch_down         = b$switchers[d1 == 1L, .N],
    multiple_switchers  = s[n_changes >= 2L, .N],
    reverters_in_design = s[n_down >= 1L, .N],
    reverters_annual    = s[reverts_annual == TRUE, .N],
    # A switcher identifies nothing unless the wave BEFORE its first change is
    # observed, so `supports_k1` is the count the estimator actually uses, not
    # `switchers`. The second column says who pays: a school first switching in
    # 2023 needs a 2022 record, and 2022 is the sample wave, so keeping it in the
    # calendar is what removes most of the 2023 cohort from the first horizon.
    no_horizon          = length(lost),
    no_horizon_2023     = s[school_id %in% lost & first_treated_year == 2023L, .N],
    supports_k1         = hz[horizon == 1L, .N],
    supports_k2         = hz[horizon == 2L, .N],
    supports_k3         = hz[horizon == 3L, .N])
})
fwrite(support, file.path(DIR_TABLES, "06m_calendar_support.csv"))

# ---- First-adoption cohorts --------------------------------------------------
# Dated on the annual history, so a school's cohort does not change with the
# calendar; what changes is which cohorts the calendar can observe and use.
cohorts <- each(function(dn, nm) {
  s <- merge(B[[dn]][[nm]]$schools[, .(school_id, F)],
             ann_attr[[dn]][, .(school_id, first_treated_year)], by = "school_id")
  s[, .(definition = dn, design = nm, schools = .N, switchers = sum(!is.na(F))),
    by = .(cohort = first_treated_year)]
})
cohorts[, pct_switchers := round(100 * switchers / sum(switchers), 1),
        by = .(definition, design)]
setorder(cohorts, definition, design, cohort, na.last = TRUE)
fwrite(cohorts, file.path(DIR_TABLES, "06m_calendar_cohorts.csv"))

# ---- Calendar time inside each Effect_k --------------------------------------
# The decisive table: how homogeneous the calendar time since adoption actually is
# within the horizon the estimator reports as one number. `modal_pct` is the share
# of the largest class, so a calendar solves the problem only if it is high. The
# accumulated-exposure columns sit alongside because a homogeneous horizon in
# calendar time need not be homogeneous in treatment received.
#
# Those two columns are computed over ADOPTERS only (0->1). A switcher in the
# other direction has no treated years by construction, so pooling the two would
# report a mechanical zero as if it were attrition of exposure.
expo <- each(function(dn, nm)
  B[[dn]][[nm]]$horizons[, .(definition = dn, design = nm, switchers = .N,
                       pct_adopters = round(100 * mean(d1 == 0L), 1),
                       mean_treated_years_cumulative = round(mean(treated_years[d1 == 0L]), 2),
                       pct_treated_throughout = round(100 * mean(treated_throughout[d1 == 0L]), 1)),
                   by = .(horizon, years_since_adoption)])
expo[, pct := round(100 * switchers / sum(switchers), 1),
     by = .(definition, design, horizon)]
expo[, modal_pct := max(pct), by = .(definition, design, horizon)]
setcolorder(expo, c("definition", "design", "horizon", "years_since_adoption", "switchers",
                    "pct", "modal_pct"))
setorder(expo, definition, design, horizon, years_since_adoption)
fwrite(expo, file.path(DIR_TABLES, "06m_exposure_by_horizon.csv"))

# ---- Cohort composition of each exposure class -------------------------------
# Calendar time since adoption and adoption cohort are not separable by design: a
# class that is one cohort measures that cohort, not that elapsed time. This is where the 2023
# cohort, which carries both the mathematics positive and the pre-trend problem,
# becomes visible per class.
comp <- each(function(dn, nm) {
  z <- merge(B[[dn]][[nm]]$horizons[, .(school_id, horizon, years_since_adoption)],
             ann_attr[[dn]][, .(school_id, first_treated_year)], by = "school_id")
  z[, .(definition = dn, design = nm, schools = .N),
    by = .(horizon, years_since_adoption, cohort = first_treated_year)]
})
comp[, pct := round(100 * schools / sum(schools), 1),
     by = .(definition, design, horizon, years_since_adoption)]
setorder(comp, definition, design, horizon, years_since_adoption, cohort, na.last = TRUE)
fwrite(comp, file.path(DIR_TABLES, "06m_cohort_by_exposure.csv"))

# ---- Save the classification -------------------------------------------------
# 06n and 06p read this instead of re-deriving the classes, so no two scripts can
# describe different schools. `exposure` is the class at the first horizon, which
# is what the partition is built on. It is calendar time, NOT accumulated exposure.
schools <- each(function(dn, nm) {
  s  <- B[[dn]][[nm]]$schools[, .(definition = dn, design = nm, school_id,
                                  first_switch_wave = F, n_waves)]
  e1 <- B[[dn]][[nm]]$horizons[horizon == 1L, .(school_id, switch_year,
                                                years_since_adoption)]
  merge(s, e1, by = "school_id", all.x = TRUE)
})
horizons <- each(function(dn, nm)
  cbind(definition = dn, design = nm,
        B[[dn]][[nm]]$horizons[, .(school_id, horizon, years_since_adoption, treated_years)]))
saveRDS(list(designs = DESIGNS, defs = vapply(DEFS, `[[`, character(1), "tv"),
             schools = schools, horizons = horizons,
             cohorts = rbindlist(lapply(names(DEFS), function(dn)
               cbind(definition = dn, ann_attr[[dn]])))),
        file.path(DIR_PROCESSED, "calendar_designs.rds"))
