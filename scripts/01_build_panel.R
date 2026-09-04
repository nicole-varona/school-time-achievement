# ==============================================================================
# 01_build_panel.R — Build the school × year analysis panel
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : Raw RA census files (2011-2025) and Aprender grade-6 microdata
#           (2016, 2018, 2021, 2022, 2023, 2025). See 00_setup.R for layout.
# Outputs : data/processed/ra_panel.rds        school × year, RA variables
#           data/processed/aprender_school.rds school × wave, Aprender vars
#           data/processed/panel.rds           merged analysis panel
#           data/processed/panel.csv           human-readable copy of panel.rds
#           data/processed/data_dictionary.csv variable documentation
#           output/tables/01_match_rate_by_wave.csv  Aprender->RA match rate
#
# Pipeline: (1) RA: treatment counts + school covariates + grade-6 trajectory,
#           one row per school-year, 2011-2025. Restricted to primary schools
#           (the thesis studies grade 6 only; see the filter below).
#           (2) Aprender: student microdata collapsed to school level
#           (unweighted means, per official methodology; label-based
#           aggregation for SPSS waves, code-based for the 2025 CSV).
#           (3) Merge on (school_id, year); save with English variable names.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))

# ============================================================================
# 1. RA panel: treatment + covariates + grade-6 trajectory (2011-2025)
# ============================================================================

build_ra_year <- function(year) {
  ydir <- ra_year_dir(year)
  if (is.na(ydir)) { return(NULL) }

  pop <- read_ra(find_file(ydir, RA_FILES$population))
  enr <- read_ra(find_file(ydir, RA_FILES$enrolment, exclude = "edad"))
  if (is.null(pop) || is.null(enr)) {
    return(NULL)
  }

  # --- Population base: treatment counts + school demographics ----------------
  n_ext <- num_matched(pop, RA_PATTERNS$extended_day)
  n_ful <- num_matched(pop, RA_PATTERNS$full_day)
  if (is.null(n_ext) || is.null(n_ful))
    stop(sprintf("[ra] %d: extended/full-day columns not found.", year))
  d_pop <- data.table(
    school_id      = pop$school_id,
    n_extended_day = fifelse(is.na(n_ext), 0, n_ext),
    n_full_day     = fifelse(is.na(n_ful), 0, n_ful)
  )
  opt <- function(v, n) if (is.null(v)) rep(NA_real_, n) else v
  d_pop[, `:=`(
    # "Turno Doble" is an "X"/"" marker in the RA, not a count -> 0/1 flag
    double_shift = opt(flag_matched(pop, RA_PATTERNS$double_shift), nrow(pop)),
    n_indigenous = opt(num_matched(pop, RA_PATTERNS$indigenous),   nrow(pop)),
    n_computer   = opt(num_matched(pop, RA_PATTERNS$computer),     nrow(pop)),
    n_disability = sum_matched(pop, RA_PATTERNS$disability),  # 3 cols summed
    n_meals      = sum_matched(pop, RA_PATTERNS$meals),       # 5 cols summed
    # Lunch on its own. n_meals sums all five services, so it cannot say whether a
    # school feeds pupils at midday -- which is the material signature of a day that
    # crosses noon, and the evidence that a five-hour 7:45-12:45 schedule is not a
    # genuine extension. Carried separately so that argument is
    # reproducible from the pipeline rather than computed by hand.
    n_lunch      = sum_matched(pop, RA_PATTERNS$lunch)
  )]
  d_pop[!is.na(n_meals) & n_meals < 0, n_meals := NA_real_]   # source error

  # --- Enrolment base: reference vars + denominator (grades 1-7) --------------
  enroll <- sum_matched(enr, RA_PATTERNS$enrolment_grades_1_7)
  if (all(is.na(enroll)))
    stop(sprintf("[ra] %d: grade enrolment columns not found.", year))
  d_enr <- data.table(
    school_id      = enr$school_id,
    enroll_primary = enroll,
    enroll_male    = sum_matched(enr, RA_PATTERNS$enrolment_male_1_7),
    repeaters_g6   = opt(num_matched(enr, RA_PATTERNS$repeaters_g6), nrow(enr)),
    overage_g6     = opt(num_matched(enr, RA_PATTERNS$overage_g6),   nrow(enr))
  )
  # Per-grade enrolment (kept so treatment counts can be matched against exact
  # grade combinations: which grades does the extended day cover?)
  for (g in 1:7)
    d_enr[, (sprintf("enroll_g%d", g)) :=
            opt(num_matched(enr, sprintf("^x[._]+%d$", g)), nrow(enr))]
  for (v in c("province", "sector", "area")) {
    src <- get_col(enr, c(province = "provincia", sector = "sector",
                          area = "ambito")[[v]])
    if (!is.null(src)) d_enr[, (v) := trimws(src)]
  }

  # --- Trajectory base (grade 6; year t records school year t-1) --------------
  trj <- read_ra(find_file(ydir, RA_FILES$trajectory))
  d_trj <- NULL
  if (!is.null(trj) && "school_id" %in% names(trj)) {
    d_trj <- data.table(
      school_id        = trj$school_id,
      promoted_g6      = opt(num_matched(trj, RA_PATTERNS$promoted_g6),      nrow(trj)),
      promoted_exam_g6 = opt(num_matched(trj, RA_PATTERNS$promoted_ex_g6),   nrow(trj)),
      not_promoted_g6  = opt(num_matched(trj, RA_PATTERNS$not_promoted_g6),  nrow(trj)),
      left_transfer_g6 = opt(num_matched(trj, RA_PATTERNS$left_transfer_g6), nrow(trj)),
      dropout_g6       = opt(num_matched(trj, RA_PATTERNS$dropout_g6),       nrow(trj))
    )
  }

  out <- merge(d_pop, d_enr, by = "school_id", all = TRUE)
  if (!is.null(d_trj)) out <- merge(out, d_trj, by = "school_id", all.x = TRUE)
  out[, year := year]
  out[]
}

ra_panel <- rbindlist(lapply(RA_YEARS, build_ra_year), fill = TRUE)

# ---- Restrict to primary schools (the thesis universe) ----------------------
# The dissertation studies primary grade 6 only. The RA population base also
# lists initial- and secondary-only schools (no grade 1-7 enrolment): these have
# no primary treatment and never sit the grade-6 Aprender assessment, so they
# would only dilute the descriptive statistics and add spurious untreated
# controls. Keep every school that is primary in at least one RA year, retaining
# its full annual history (intermittent zero-primary years are left in place so
# they do not fragment the adoption timeline); drop schools that are never
# primary.
primary_schools <- ra_panel[!is.na(enroll_primary) & enroll_primary > 0,
                            unique(school_id)]
ra_panel <- ra_panel[school_id %in% primary_schools]

# PROVINCE ONLY. The 2020 census file ships double-encoded at source (UTF-8
# bytes saved as Latin-1), so five provinces arrive under a second, corrupted
# spelling in that year alone. R treats the two spellings as different strings,
# which silently splits a grouping key: province is used for clustering, for the
# regional cuts and for the provincial tables, and a split category surfaces as
# a spurious extra jurisdiction. Resolved by MODE across the school's own years
# rather than by matching the corrupted strings — one bad year cannot outvote
# fourteen good ones, and the repair uses no string literals, so it does not
# depend on the session locale (a literal-matching version failed for exactly
# that reason).
#
# NOT applied to sector or area, deliberately: their within-school variation is
# REAL, not a defect. Checked on the pre-fix panel — 248 of the 258 schools whose
# area changes do so exactly once and the new value persists (a locality
# reclassified urban/rural), and 48 of 50 sector changes are single, persistent
# transitions (a school moving between private and state management). Collapsing
# those to a mode would erase genuine reclassification: it altered 462 school-
# waves of area before this was caught.
md <- ra_panel[!is.na(province), .N, by = .(school_id, province)][
  order(school_id, -N)][, .SD[1], by = school_id]
setnames(md, "province", "..mode")
ra_panel[md, province := i...mode, on = "school_id"]

# English recodes for sector / area (raw values kept if unmatched)
lo <- function(x) tolower(trimws(x))
ra_panel[, sector := fifelse(startsWith(lo(sector), "est"), "public",
                     fifelse(startsWith(lo(sector), "pri"), "private", sector))]
ra_panel[, area   := fifelse(startsWith(lo(area), "urb"), "urban",
                     fifelse(startsWith(lo(area), "rur"), "rural", area))]

# Derived treatment variables
ra_panel[, n_expanded := n_extended_day + n_full_day]
ra_panel[, share_raw := fifelse(!is.na(enroll_primary) & enroll_primary > 0,
                                n_expanded / enroll_primary, NA_real_)]

# The numerator comes from Base 6-Poblacion and the denominator from Base 2-
# Matricula. The two bases disagree by a few pupils for some fully-covered
# schools in 2011-2012, and a handful of records enter the same count under both
# regimes, so the raw ratio can exceed 1. The share is a proportion by
# definition, so it is capped; every affected record already clears the 0.50
# threshold, so no treatment status changes.
# The discrepancy is confined to 164 school-years in 2011-2012. Stop rather than
# cap silently if a later RA wave carries the same problem.
OVER_YEARS <- c(2011L, 2012L)
OVER_N     <- 164L
over <- ra_panel[!is.na(share_raw) & share_raw > 1]
if (!all(over$year %in% OVER_YEARS) || nrow(over) > OVER_N)
  stop("share_expanded exceeds 1 beyond the known records: ", nrow(over),
       " school-years in ", paste(sort(unique(over$year)), collapse = ", "),
       " (expected at most ", OVER_N, " in ",
       paste(OVER_YEARS, collapse = ", "),
       "). Reconcile Base 6-Poblacion against Base 2-Matricula before capping.")

ra_panel[, share_expanded := pmin(share_raw, 1)]
ra_panel[, share_raw := NULL]

# Grade-6 trajectory rates (denominator: all grade-6 outcomes recorded)
traj_cols <- c("promoted_g6", "promoted_exam_g6", "not_promoted_g6",
               "left_transfer_g6", "dropout_g6")
if (all(traj_cols %in% names(ra_panel))) {
  ra_panel[, traj_total_g6 := rowSums(.SD, na.rm = TRUE), .SDcols = traj_cols]
  ra_panel[, `:=`(
    nonpromotion_rate_g6 = fifelse(traj_total_g6 > 0,
                                   not_promoted_g6 / traj_total_g6, NA_real_),
    dropout_rate_g6      = fifelse(traj_total_g6 > 0,
                                   dropout_g6 / traj_total_g6, NA_real_)
  )]
}

# ============================================================================
# 2. Aprender: collapse student microdata to school level (per wave)
# ============================================================================
# School means are unweighted, per the official Aprender methodology (survey
# weights correct jurisdiction/sector totals, not within-school averages).

# Performance levels -> standardized code 1-4 (1 = below basic ... 4 = advanced)
level_code <- function(x, coded) {
  if (is.null(x)) return(NULL)
  if (coded) { v <- to_num(x); v[!is.na(v) & (v < 1 | v > 4)] <- NA_real_; return(v) }
  lab <- tolower(labels_of(x))
  out <- rep(NA_real_, length(lab))
  out[grepl("debajo", lab)] <- 1
  out[grepl("b.sico", lab) & !grepl("debajo", lab)] <- 2
  out[grepl("satisfact", lab)] <- 3
  out[grepl("avanz", lab)] <- 4
  out
}

# Label match preserving NA (grepl maps NA -> FALSE, which would bias shares)
match_label <- function(x, pattern, invert = FALSE) {
  if (is.null(x)) return(NULL)
  lab <- labels_of(x)
  res <- grepl(pattern, tolower(lab))
  res[is.na(lab) | lab == "NA"] <- NA
  if (invert) !res else res
}

build_aprender_wave <- function(year) {
  ydir <- aprender_dir(year)
  f <- find_file(ydir, include = c("estudiante", "primaria"),
                 exclude = "secundaria", prefer_not = "muestral")
  if (is.na(f)) { return(NULL) }
  ap <- read_aprender(f)
  if (!"school_id" %in% names(ap)) stop(sprintf("[aprender] %d: school ID not found.", year))

  map   <- APRENDER_MAP[[as.character(year)]]
  coded <- identical(map$format, "csv_coded")

  score_lang <- to_num(get_col(ap, APRENDER_OUTCOME$score_lang))
  score_math <- to_num(get_col(ap, APRENDER_OUTCOME$score_math))
  if (is.null(score_lang) || is.null(score_math))
    stop(sprintf("[aprender] %d: score columns not found.", year))

  stu <- data.table(
    school_id  = ap$school_id,
    score_lang = score_lang,
    score_math = score_math,
    lev_lang   = level_code(get_col(ap, APRENDER_OUTCOME$level_lang), coded),
    lev_math   = level_code(get_col(ap, APRENDER_OUTCOME$level_math), coded)
  )

  # --- Composition items (wave-specific coding) -------------------------------
  add <- function(name, v) if (!is.null(v)) stu[, (name) := v]
  if (coded) {
    # Standard Aprender codings (verified against the official 2025 dictionary).
    sexv <- to_num(get_col(ap, map$sex))
    # 1 = male, 2 = female, 3 = non-binary ("X" in the DNI): a distinct category,
    # kept as not-female (stays in the denominator), NOT treated as missing.
    if (!is.null(sexv)) add("is_female", sexv == 2)
    prev <- to_num(get_col(ap, map$preschool))
    if (!is.null(prev) && any(!is.na(prev))) {
      yes <- if (max(prev, na.rm = TRUE) <= 2) prev == 1 else prev < 4
      add("went_preschool", yes)
    }
    repv <- to_num(get_col(ap, map$repeated))
    if (!is.null(repv)) add("repeated_grade", repv >= 2)  # 1 = never repeated
    sesv <- get_col(ap, map$ses)                      # "Q1".."Q5" -> 1-5
    if (!is.null(sesv)) add("ses_index", to_num(sub("^\\s*Q", "", as.character(sesv))))
    # Parental education: constructed 1-7 scale (this wave has no self-report item
    # and no "No se" code). >= 4 = secondary-incomplete or higher, matching the
    # label patterns used for the .sav waves.
    med <- to_num(get_col(ap, map$mother_ed)); if (!is.null(med)) add("mother_ed_sec", med >= 4)
    fed <- to_num(get_col(ap, map$father_ed)); if (!is.null(fed)) add("father_ed_sec", fed >= 4)
  } else {
    add("is_female",      match_label(get_col(ap, map$sex), "mujer|femenin"))
    add("went_preschool", match_label(get_col(ap, map$preschool),
                                      "no fui|\\bno\\b", invert = TRUE))
    add("repeated_grade", match_label(get_col(ap, map$repeated),
                                      "nunca|\\bno\\b", invert = TRUE))
    add("mother_ed_sec",  match_label(get_col(ap, map$mother_ed),
                                      "secundari|superior|univ|terciari|posgrado"))
    add("father_ed_sec",  match_label(get_col(ap, map$father_ed),
                                      "secundari|superior|univ|terciari|posgrado"))
    ses <- get_col(ap, map$ses); if (!is.null(ses)) add("ses_index", to_num(ses))
  }
  add("age_code",     if (!is.na(map$age))     to_num(get_col(ap, map$age)))
  add("absence_code", if (!is.na(map$absence)) to_num(get_col(ap, map$absence)))
  add("climate_index",   if (!is.na(map$climate))   to_num(get_col(ap, map$climate)))
  add("context_quartile", if (!is.na(map$context_q)) to_num(get_col(ap, map$context_q)))

  # --- Collapse to school level (unweighted) ----------------------------------
  school <- stu[, {
    res <- list(
      n_tested            = .N,
      score_lang          = mean_safe(score_lang),
      score_math          = mean_safe(score_math),
      pct_below_basic_lang = mean_safe(lev_lang == 1),
      pct_basic_lang       = mean_safe(lev_lang == 2),
      pct_satisf_lang      = mean_safe(lev_lang == 3),
      pct_advanced_lang    = mean_safe(lev_lang == 4),
      pct_satisf_plus_lang = mean_safe(lev_lang >= 3),
      pct_below_basic_math = mean_safe(lev_math == 1),
      pct_basic_math       = mean_safe(lev_math == 2),
      pct_satisf_math      = mean_safe(lev_math == 3),
      pct_advanced_math    = mean_safe(lev_math == 4),
      pct_satisf_plus_math = mean_safe(lev_math >= 3)
    )
    for (v in c("is_female", "went_preschool", "repeated_grade",
                "mother_ed_sec", "father_ed_sec")) {
      out_name <- c(is_female = "pct_female", went_preschool = "pct_preschool",
                    repeated_grade = "pct_repeater",
                    mother_ed_sec = "pct_mother_secondary",
                    father_ed_sec = "pct_father_secondary")[[v]]
      if (v %in% names(.SD)) res[[out_name]] <- mean_safe(.SD[[v]])
    }
    for (v in c("age_code", "absence_code", "ses_index", "climate_index",
                "context_quartile"))
      if (v %in% names(.SD)) res[[v]] <- mean_safe(.SD[[v]])
    res
  }, by = school_id]

  school[, `:=`(year = year, wave_type = WAVE_TYPE[[as.character(year)]])]
  rm(ap, stu); gc(verbose = FALSE)
  school[]
}

aprender_school <- rbindlist(lapply(APRENDER_WAVES, build_aprender_wave),
                             fill = TRUE)

# ============================================================================
# 3. Row order (no treatment definition lives in this script)
# ============================================================================
# NO TREATMENT INDICATOR IS DERIVED HERE. This script used to build one
# (`treated`, at a 50% cut-off) plus a school-level trajectory summary
# (`ever_treated`, `first_treated_year`, `years_treated`, `reverts`), and carry
# them into the panel. They were a fourth copy of a definition that
# 02_build_estimation_sample.R owns, and the duplication was not benign: the
# column collided by name with the one 05 merges in, so a merge produced
# `treated.x`/`treated.y` and every line downstream would have found neither.
# Nothing read them -- every consumer derives its own from `treatment_annual.rds`
# -- so they are gone. This script's job is the panel: enrolment, shares,
# covariates and outcomes. `share_expanded` is what a treatment definition is
# built FROM, and it stays; the definitions themselves live in 02.
setorder(ra_panel, school_id, year)

# ============================================================================
# 4. Merge: analysis panel (Aprender outcome waves x RA treatment/covariates)
# ============================================================================

# Unified panel: one row per school x year over the full RA history (2011-2025).
# RA variables (treatment, grade-6 trajectory, structural covariates) are present
# every year; Aprender variables (scores, composition, performance levels) are
# merged in and are NA outside the six assessment waves. A single panel serves
# every analysis, which filters to what it needs: the score/achievement analyses
# keep rows with non-missing scores (the six waves), the trajectory analysis uses
# all annual rows.
panel <- merge(ra_panel, aprender_school,
               by = c("school_id", "year"), all.x = TRUE)

# ---- Hora+ identification variables (from Aprender Director bases) -----------
# The RA cannot distinguish the "Hora+/1h+" programme (adds 1 hour to a simple
# day) from a genuine extended day. Two director-reported sources help:
#   * grade-6 daily hours (indirect: ~5h => simple + 1h = Hora+): censal waves
#     2021 (dp23), 2022 (dp18), 2023 (dp47). Merged per (school_id, year).
#   * explicit participation flag (dp21, 2024 sample): school-level, merged to
#     all of a school's waves as a validation flag.
hours_items <- list("2021" = c(dir_hours_g6 = "dp23"),
                    "2022" = c(dir_hours_g6 = "dp18"),
                    "2023" = c(dir_hours_g6 = "dp47"))
dir_hours <- rbindlist(lapply(names(hours_items), function(y)
  read_director_items(as.integer(y), hours_items[[y]])), fill = TRUE)
if (!is.null(dir_hours) && "dir_hours_g6" %in% names(dir_hours))
  panel <- merge(panel, dir_hours, by = c("school_id", "year"), all.x = TRUE)

hm24 <- read_director_items(2024, c(hora_mas_2024 = "dp21"))
if (!is.null(hm24) && "hora_mas_2024" %in% names(hm24)) {
  hm24 <- unique(hm24[!is.na(hora_mas_2024), .(school_id, hora_mas_2024)])
  panel <- merge(panel, hm24, by = "school_id", all.x = TRUE)
}

# Aprender -> RA match, computed on the assessment-wave rows only (n_tested is
# the Aprender per-school tested count; NA outside the six waves).
match_diag <- panel[!is.na(n_tested), .(
  schools     = .N,
  matched_ra  = sum(!is.na(share_expanded)),
  match_rate  = round(100 * mean(!is.na(share_expanded)), 2)
), by = .(year, wave_type)][order(year)]

# ============================================================================
# 5. Save processed data + data dictionary
# ============================================================================

saveRDS(ra_panel,        file.path(DIR_PROCESSED, "ra_panel.rds"))
saveRDS(aprender_school, file.path(DIR_PROCESSED, "aprender_school.rds"))
saveRDS(panel,           file.path(DIR_PROCESSED, "panel.rds"))
fwrite(panel,            file.path(DIR_PROCESSED, "panel.csv"))  # human-readable copy
fwrite(match_diag,       file.path(DIR_TABLES, "01_match_rate_by_wave.csv"))

dictionary <- data.table(rbind(
  c("school_id",       "15-digit fictitious school ID (links RA and Aprender)"),
  c("year",            "Calendar year (RA wave / Aprender wave)"),
  c("province",        "Province (24 jurisdictions)"),
  c("sector",          "School sector: public / private"),
  c("area",            "urban / rural"),
  c("enroll_primary",  "Primary enrolment, grades 1-7 (treatment denominator)"),
  c("enroll_g1 ... enroll_g7", "Enrolment by grade (to match treatment counts against grade combinations)"),
  c("enroll_male",     "Male primary enrolment, grades 1-7"),
  c("n_extended_day",  "Students in extended-day schooling, primary (RA)"),
  c("n_full_day",      "Students in full-day schooling, primary (RA)"),
  c("n_expanded",      "n_extended_day + n_full_day"),
  c("share_expanded",  "n_expanded / enroll_primary (treatment intensity)"),
  c("double_shift",    "School offers double shift in primary (0/1 flag; RA marker column). Validates n_full_day"),
  c("n_indigenous",    "Indigenous students, primary"),
  c("n_computer",      "Students with computer at home, primary"),
  c("n_disability",    "Students with disability, primary (3 categories summed)"),
  c("n_meals",         "Free-meal beneficiaries (5 meal types summed; SES proxy)"),
  c("repeaters_g6",    "Grade-6 repeaters (RA enrolment base)"),
  c("overage_g6",      "Grade-6 over-age students (RA; SES proxy)"),
  c("promoted_g6",     "Grade-6 promoted (trajectory, refers to year t-1)"),
  c("promoted_exam_g6","Grade-6 promoted after exam"),
  c("not_promoted_g6", "Grade-6 not promoted"),
  c("left_transfer_g6","Grade-6 left with transfer pass"),
  c("dropout_g6",      "Grade-6 left without pass (dropout)"),
  c("nonpromotion_rate_g6", "not_promoted_g6 / all grade-6 outcomes"),
  c("dropout_rate_g6", "dropout_g6 / all grade-6 outcomes"),
  c("n_tested",        "Aprender: students tested in the school-wave"),
  c("score_lang",      "Aprender language score, school mean (TRI scale, 2016=500/100)"),
  c("score_math",      "Aprender mathematics score, school mean"),
  c("pct_below_basic_lang", "Share below basic level, language (0-1)"),
  c("pct_satisf_plus_lang", "Share satisfactory or advanced, language"),
  c("pct_female",      "Share female students (tested)"),
  c("pct_preschool",   "Share who attended preschool"),
  c("pct_repeater",    "Share who ever repeated a grade"),
  c("pct_mother_secondary", "Share with mother >= secondary education (all waves)"),
  c("pct_father_secondary", "Share with father >= secondary education (all waves)"),
  c("ses_index",       "Student SES index, school mean (wave-specific scale!)"),
  c("climate_index",   "School climate index (2016, 2023 only)"),
  c("dir_hours_g6",    "Director-reported daily hours for grade 6 (2021/2022/2023): indirect Hora+ signal (~5h)"),
  c("hora_mas_2024",   "School reports participating in the Hora+/1h+ programme (explicit item, 2024 sample)"),
  c("wave_type",       "census / sample (2022 is a school-level sample)")
))
setnames(dictionary, c("variable", "description"))
fwrite(dictionary, file.path(DIR_PROCESSED, "data_dictionary.csv"))
