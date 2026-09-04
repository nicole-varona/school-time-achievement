# ==============================================================================
# 92_estimate_registry.R — One registry of every estimate the project produced
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : the estimate tables already written by 04, 06*, 08, 09, 10, 11*
# Outputs : output/tables/92_estimate_registry.csv
#
# Estimates were spread across eighteen files with eighteen different shapes, so
# the same quantity appeared in several places and no single file answered "what
# was estimated, and what role does it play". This reads all of them and emits
# ONE long table with a fixed schema:
#
#   role       main | estimator robustness | specification robustness |
#              diagnostic | heterogeneity (exploratory) | secondary outcome
#              -- the taxonomy of the working document. Its sixth category,
#              SENSITIVITY, has no rows here on purpose: HonestDiD bounds and the
#              equivalence tests are not average effects, so they cannot be a row
#              in this schema. They stay in 07_honestdid_*.csv and 11_fect_tests.csv.
#   estimator  CDH | CS | FEct | TWFE
#   spec       the one thing that varies relative to the main specification
#   treatment  which treatment definition
#   sample     which schools/waves
#   outcome    language | mathematics | non-promotion | dropout
#   units      effects are NOT all in the same unit; see below
#   att, se    the average effect and its standard error (se is NA where the
#              source table does not carry one)
#   on_primary TRUE only where the estimator is the primary one (CDH). A row with
#              FALSE cannot substitute for the same check on the primary.
#   note       only where a row would otherwise be misread
#   source     the file the row was read from
#
# `role`, `treatment` and `sample` are JUDGEMENTS, not data: they are assigned by
# the explicit maps below rather than inferred, so the classification is auditable
# in one place instead of being implied by prose.
#
# UNITS. Score effects are in SD of the anchored 2016 scale (/100). Two groups are
# not in that unit and are labelled: the within-wave standardised rows (SD of the
# within-wave distribution of school means -- comparable in sign, not in level)
# and the grade-progression outcomes (percentage points).
#
# DELIBERATELY NOT HARVESTED, because these are not average effects and would turn
# a registry into a coefficient dump: event studies and placebos (07, 07b, 08,
# 09), the dynamic profile (07c), HonestDiD bounds (07), FEct's test battery
# (11_fect_tests), and the per-horizon effects and placebos inside 06g, 06j, 06n
# and 06p. Those stay in their own files; this table is the index of ATT-level
# estimates.
#
# Every table that carries an ATT is now read. Seven were neither harvested nor
# excluded until 13/08 -- they were simply absent, which is the failure mode a
# registry exists to prevent: 06j (the exposure classes of @tbl-main Panel A and
# all of @tbl-timing-placebo), 06k (the leave-one-province check quoted in the
# inference section), 06n, 06p, 06q, 11_fect_carryover and 11b_dynamic_matched.
# The completeness check at the foot of this script now asserts it rather than
# leaving it to be noticed.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))

# Read a source table, failing loudly: a missing file would silently shrink the
# registry, which is worse than not having one.
src <- function(f) {
  p <- file.path(DIR_TABLES, f)
  if (!file.exists(p)) stop("missing source table: ", f, " (run the pipeline first)")
  fread(p)
}

SD16 <- "SD (anchored 2016 scale)"

# Subject labels differ across files ("language", "Language", "y_lang", "ya_lang").
subj <- function(x) {
  s <- tolower(as.character(x))
  fifelse(s %like% "lang", "language",
  fifelse(s %like% "math", "mathematics", NA_character_))
}

# Estimator labels differ across files ("CS (robustness)", "Callaway-Sant'Anna",
# "CDH (did_multiplegt_dyn)"). Collapse to the four families plus the Stata command.
norm_est <- function(x) {
  s <- as.character(x)
  fifelse(s %like% "CDH|multiplegt_dyn", "CDH",
  fifelse(s %like% "^CS|Callaway",       "CS",
  fifelse(s %like% "FEct",               "FEct",
  fifelse(s %like% "TWFE",               "TWFE", NA_character_))))
}

reg <- function(role, estimator, spec, treatment, sample, outcome, att, se = NA_real_,
                units = SD16, note = NA_character_, source)
  data.table(role, estimator, spec, treatment, sample, outcome,
             units, att = as.numeric(att), se = as.numeric(se), note, source)

R <- list()

# ---- 06_estimation_main.csv --------------------------------------------------
# The main table. Role, treatment and sample are keyed on the three columns that
# identify a specification there.
m <- src("06_estimation_main.csv")
MAIN_MAP <- data.table(
  estimator = c("CDH (headline)", "CDH (headline)", "CS (absorbing)", "CS (absorbing)",
                "CS (absorbing)", "CS (absorbing)", "CS (absorbing)",
                "CDH (enrol-weighted)", "CDH (never-switchers)",
                "CS (enrol-weighted)", "TWFE (benchmark)"),
  regime = c("combined", "full-day", rep("combined", 5), "combined", "combined",
             "combined", "combined"),
  covariates = c("none (absorbed)", "none (absorbed)", "none", "structural",
                 "none (comp-sample)", "composition", "composition, no SES",
                 "none (absorbed)", "none (absorbed)",
                 "composition", "none (absorbed)"),
  out_estimator = c("CDH", "CDH", rep("CS", 5), "CDH", "CDH", "CS", "TWFE"),
  # Only ONE of the five CS rows is estimator robustness: the CS reference spec,
  # which is Callaway-Sant'Anna's answer to the main question. The other four vary
  # the COVARIATE SET inside the robustness estimator, so they are specification
  # robustness OF CS -- and by the working document's own rule they cannot stand as
  # robustness of the headline. `on_primary` below makes that readable at a glance.
  role = c("main", "specification robustness",
           "specification robustness",   # no covariates
           "specification robustness",   # structural only
           "specification robustness",   # no covariates, composition subsample
           "estimator robustness",       # composition + SES = the CS reference
           "specification robustness",   # composition without SES
           "specification robustness", "specification robustness",
           "specification robustness", "estimator robustness"),
  spec = c("reference specification", "full-day regime only",
           "no covariates", "structural covariates only",
           "no covariates, composition-complete subsample",
           "composition + SES covariates (CS reference)", "composition without SES",
           "enrolment-weighted (pupil-average estimand)",
           "never-switchers as comparison group",
           "enrolment-weighted (pupil-average estimand)", "static benchmark"),
  treatment = c("combined >=50%", "full-day >=50%",
                rep("combined >=50%", 5), "combined >=50%",
                "combined >=50%", "combined >=50%", "combined >=50%"),
  sample = c("estimation sample", "estimation sample (full-day defined)",
             "absorbing subsample", "absorbing subsample",
             "absorbing subsample, composition-complete",
             "absorbing subsample, composition-complete",
             "absorbing subsample, composition-complete",
             "estimation sample", "estimation sample",
             "absorbing subsample, composition-complete", "estimation sample"),
  # Two rows left this map on 13/08 with the specifications they described: the
  # continuous dose and the quasi value-added spec, both removed on 11/08. They
  # matched nothing, so they cost nothing -- but a map that still names retired
  # specifications is a map nobody can read against the pipeline.
  note = c(NA, NA, NA, NA, NA, NA,
           "cannot be fitted: singular (g,t) cell at the 2022 sample-wave cohort",
           NA, NA, NA,
           "biased under staggered timing with reversals; reference only"))
# Appended rather than inserted: the map is joined by key, not by position, and a
# single-row rbind cannot misalign the nine parallel vectors above.
MAIN_MAP <- rbind(MAIN_MAP, data.table(
  estimator  = c("CDH (grade-6-weighted)", "CS (grade-6-weighted)"),
  regime     = "combined",
  covariates = c("none (absorbed)", "composition"),
  out_estimator = c("CDH", "CS"),
  role       = "specification robustness",
  spec       = "grade-6-enrolment-weighted (pupil-average estimand)",
  treatment  = "combined >=50%",
  sample     = c("estimation sample", "absorbing subsample, composition-complete"),
  note       = NA_character_))

m <- merge(m, MAIN_MAP, by = c("estimator", "regime", "covariates"), all.x = TRUE)
if (anyNA(m$role)) stop("unmapped specification in 06_estimation_main.csv")
R$main <- reg(m$role, m$out_estimator, m$spec, m$treatment, m$sample, subj(m$subject),
              m$att, m$se, note = m$note, source = "06_estimation_main.csv")

# ---- 04_threshold_sensitivity_cdh.csv ---------------------------------------
t4 <- src("04_threshold_sensitivity_cdh.csv")
R$thr <- reg("specification robustness", "CDH", "treatment cut-off",
             t4$threshold, "estimation sample", subj(t4$subject), t4$att, t4$se,
             note = fifelse(t4$threshold %like% "headline", "same as main", NA_character_),
             source = "04_threshold_sensitivity_cdh.csv")

# ---- 06c_covid_robustness.csv (wide: two estimators x two subjects) ---------
cv <- src("06c_covid_robustness.csv")
cv_long <- rbindlist(list(
  data.table(est = "CDH", outcome = "language",    att = cv$cdh_lang, se = NA_real_),
  data.table(est = "CDH", outcome = "mathematics", att = cv$cdh_math, se = NA_real_),
  data.table(est = "CS",  outcome = "language",    att = cv$cs_lang,  se = cv$cs_lang_se),
  data.table(est = "CS",  outcome = "mathematics", att = cv$cs_math,  se = cv$cs_math_se)
))[, `:=`(spec = rep(cv$spec, 4), waves = rep(cv$waves, 4))]
R$covid <- reg("diagnostic", cv_long$est, cv_long$spec, "combined >=50%",
               paste("waves:", cv_long$waves), cv_long$outcome, cv_long$att, cv_long$se,
               note = fifelse(cv_long$spec %like% "headline", "same as main",
                       "changes the estimand: dropping waves re-dates cohorts and re-spaces event time"),
               source = "06c_covid_robustness.csv")

# ---- 06d_scale_robustness.csv -----------------------------------------------
sc <- src("06d_scale_robustness.csv")
R$scale <- reg("specification robustness", "CDH", paste("outcome scale:", sc$scale),
               "combined >=50%", "estimation sample", subj(sc$subject),
               sc$att, sc$se,
               units = fifelse(sc$scale %like% "anchored", SD16,
                               "SD (within-wave school-mean distribution)"),
               note = fifelse(sc$scale %like% "anchored", "same as main",
                              "not comparable in level with the anchored rows"),
               source = "06d_scale_robustness.csv")

# ---- 06e: grade alignment ---------------------------------------------------
ga <- src("06e_grade_alignment.csv")
R$grade <- reg("specification robustness", "CDH", paste("grade-6 exposure:", ga$sample),
               "combined >=50%", ga$sample, subj(ga$subject), ga$att,
               note = fifelse(ga$sample == "full", "same as main",
                              "source table carries no standard error"),
               source = "06e_grade_alignment.csv")

ps <- src("06e_by_primary_structure.csv")
R$struct <- reg("specification robustness", "CDH",
                paste("primary structure:", rep(ps$structure, 2)),
                "combined >=50%", paste("primary structure:", rep(ps$structure, 2)),
                rep(c("language", "mathematics"), each = nrow(ps)),
                c(ps$cdh_lang, ps$cdh_math),
                note = "source table carries no standard error",
                source = "06e_by_primary_structure.csv")

# ---- 06f: clustering level and province-specific trends ---------------------
cl <- src("06f_clustering.csv")
R$clus <- reg("specification robustness", norm_est(cl$estimator),
              paste("standard errors clustered at:", cl$cluster_level),
              "combined >=50%",
              fifelse(cl$estimator %like% "CS", "absorbing subsample, composition-complete",
                      "estimation sample"),
              subj(cl$subject), cl$att, cl$se,
              note = fifelse(cl$cluster_level == "school", "same as main",
                             "sensitivity only: 8.1 effective clusters at province, over 24"),
              source = "06f_clustering.csv")

pt <- src("06f_province_trends.csv")
R$ptrend <- reg("specification robustness", "CDH",
                fifelse(pt$spec == "main", "reference specification",
                        "province-specific year trends (trends_nonparam)"),
                "combined >=50%", "estimation sample", subj(pt$subject),
                pt$att, pt$se,
                note = fifelse(pt$spec == "main", "same as main",
                               "varies the identifying assumption: parallel trends within province-year"),
                source = "06f_province_trends.csv")

# ---- 06g: aggregation window and three estimator options -------------------
# ATE rows only; the per-horizon effects and placebos stay in the source file.
hw <- src("06g_horizon_sensitivity.csv")[parameter == "ATE"]
R$window <- reg("specification robustness", "CDH",
                paste("aggregation window:", hw$spec),
                fifelse(hw$regime == "full-day", "full-day >=50%",
                        "combined >=50%"),
                fifelse(hw$regime == "full-day", "estimation sample (full-day defined)",
                        "estimation sample"),
                subj(hw$subject), hw$estimate, hw$se,
                note = fifelse(hw$spec %like% "^2" & hw$regime == "combined", "same as main",
                               NA_character_),
                source = "06g_horizon_sensitivity.csv")

op <- src("06g_cdh_options.csv")[parameter == "ATE"]
OPT_ROLE <- c("3 horizons, all switchers" = "specification robustness",
              "3 horizons, same switchers" = "specification robustness",
              "2 horizons, same switchers" = "specification robustness",
              "2 horizons, same switchers, pupil-weighted" = "specification robustness",
              "more granular demeaning" = "specification robustness",
              "main (2 horizons)" = "specification robustness",
              "group-specific linear trends" = "specification robustness",
              "complete-case, uncontrolled" = "specification robustness",
              "complete-case, composition controls" = "specification robustness",
              "SES-stratum-specific year trends" = "specification robustness")
if (!all(op$spec %in% names(OPT_ROLE))) stop("unmapped spec in 06g_cdh_options.csv")
R$opts <- reg(OPT_ROLE[op$spec], "CDH", op$spec, "combined >=50%",
              fifelse(op$spec %like% "complete-case",
                      "estimation sample, composition-complete", "estimation sample"),
              subj(op$subject), op$estimate, op$se,
              note = fifelse(op$spec %like% "linear trends",
                      "average cumulative effect is not defined under trends_lin; per-horizon effects in the source file",
                      fifelse(op$spec %like% "^main", "same as main", NA_character_)),
              source = "06g_cdh_options.csv")

# ---- 06h_ses_covariate.csv (how the SES covariate is built) -----------------
# Only the estimate rows; the coverage rows in that file are counts, not effects.
sc <- src("06h_ses_covariate.csv")[!is.na(att)]
R$sescov <- reg("specification robustness", "CS",
                paste("SES covariate:", sc$ses_covariate), "combined >=50%",
                "absorbing subsample, composition-complete", sc$subject, sc$att, sc$se,
                note = fifelse(sc$ses_covariate %like% "2016",
                        "previous construction, kept so the choice is reportable",
                        "same as the CS reference spec"),
                source = "06h_ses_covariate.csv")

# ---- 08_arm_results.csv -----------------------------------------------------
ar <- src("08_arm_results.csv")
R$arms <- reg("specification robustness", ar$estimator,
              paste("treatment arm:", ar$arm),
              fifelse(ar$arm %like% "genuine", "genuine 6h+", "1h+ programme (5h)"),
              paste("arm subsample:", ar$arm), subj(ar$subject), ar$att, ar$se,
              source = "08_arm_results.csv")

# ---- 09_trajectory_results.csv (secondary outcome, percentage points) -------
tr <- src("09_trajectory_results.csv")
R$traj <- reg("secondary outcome", fifelse(tr$spec %like% "CS", "CS", "CDH"), tr$spec,
              "combined >=50%", "annual census panel",
              fifelse(tr$outcome == "nonpromotion", "non-promotion", "dropout"),
              tr$att_pp, tr$se_pp, units = "percentage points",
              source = "09_trajectory_results.csv")

# ---- 11: FEct, non-reverting subsample, like-for-like, sample grid ----------
fe <- src("11_fect_results.csv")
R$fect <- reg("estimator robustness", "FEct",
              sub("^FEct, ?", "", sub("FEct \\((.*)\\)", "\\1", fe$estimator)),
              "combined >=50%",
              "estimation sample", subj(fe$subject), fe$att, fe$se,
              source = "11_fect_results.csv")

nr <- src("11_nonreverter_results.csv")
R$nonrev <- reg("diagnostic", norm_est(nr$estimator), "non-reverting subsample",
                "combined >=50%", "non-reverting schools", subj(nr$subject),
                nr$att, nr$se,
                note = "conditions on post-treatment behaviour, i.e. selection",
                source = "11_nonreverter_results.csv")

lf <- src("11b_like_for_like.csv")
R$lfl <- reg("estimator robustness", norm_est(lf$estimator),
             paste("like-for-like,", lf$counterfactual_trend), "combined >=50%",
             "non-reverting schools (common sample)", subj(lf$subject), lf$att, lf$se,
             note = fifelse(lf$estimator %like% "multiplegt",
                     "aggregates THREE horizons, unlike the two of the main specification",
                     NA_character_),
             source = "11b_like_for_like.csv")

sg <- src("11c_sample_grid.csv")
R$grid <- reg("diagnostic", "FEct", paste("sample grid:", sg$sample),
              "combined >=50%", sg$sample, subj(sg$subject),
              sg$fect_att, sg$fect_se, source = "11c_sample_grid.csv")

# ---- 06j: the exposure classes of Panel A ----------------------------------
# ATE rows only, and only the two windows the manuscript reads. `spec` is the
# aggregation window, `exposure_years` the class, so both belong in `spec` here:
# the class IS the specification, and the sample is what distinguishes the rows.
pj <- src("06j_exposure_partition.csv")[parameter == "ATE"]
R$expo <- reg("specification robustness", "CDH",
              sprintf("time since first adoption = %s years (%s)",
                      pj$exposure_years, pj$spec),
              "combined >=50%",
              sprintf("class %s + all never-switching schools", pj$exposure_years),
              subj(pj$subject), pj$estimate, pj$se,
              note = "estimated against all never-switching schools; classes are not pooled",
              source = "06j_exposure_partition.csv")

# ---- 06k: leave one province out -------------------------------------------
lop <- src("06k_leave_one_province.csv")
R$lopo <- reg("diagnostic", "CDH",
              fifelse(lop$dropped %like% "^none", "reference specification",
                      paste("dropping", lop$dropped)),
              "combined >=50%",
              fifelse(lop$dropped %like% "^none", "estimation sample",
                      paste("estimation sample less", lop$dropped)),
              subj(lop$subject), lop$att, lop$se,
              note = fifelse(lop$dropped %like% "^none", "same as main",
                             "leverage check on point estimates, not an inference procedure"),
              source = "06k_leave_one_province.csv")

# ---- 06n: the same estimand under each candidate outcome calendar ----------
cal <- src("06n_calendar_estimates.csv")[parameter == "ATE"]
R$cal <- reg("diagnostic", "CDH",
             sprintf("outcome calendar: %s (%s)", cal$design, cal$spec),
             "combined >=50%", paste("waves:", cal$design),
             subj(cal$subject), cal$estimate, cal$se,
             note = paste("changes which waves enter, which re-dates adoption and",
                          "re-spaces event time; the sign of the ATT depends on it"),
             source = "06n_calendar_estimates.csv")

# ---- 06p: first-spell length ------------------------------------------------
sp <- src("06p_spell_length_estimates.csv")[parameter == "ATE"]
R$spell <- reg("diagnostic", "CDH",
               sprintf("first spell %s (%s, %s)", sp$first_spell, sp$spec, sp$design),
               "combined >=50%",
               sprintf("switchers whose first spell lasts %s", sp$first_spell),
               subj(sp$subject), sp$estimate, sp$se,
               note = "descriptive across classes, NOT a dose-response: spell length is realised after adoption",
               source = "06p_spell_length_estimates.csv")

# ---- 06q: the headline under three inference procedures --------------------
# One estimate, three variances. Harvested at the school-clustered variance, which
# is the one the headline reports; the other two live in the source file, because
# a registry keyed on (spec, sample, outcome) cannot hold three rows that differ
# only in how the same number is bounded.
qb <- src("06q_cdh_bootstrap.csv")[parameter == "ATE"]
R$boot <- reg("specification robustness", "CDH",
              "province-clustered inference (analytic and bootstrap)",
              "combined >=50%", "estimation sample", qb$subject, qb$att, qb$se_school,
              note = paste0("same estimate as the headline; province analytic s.e. ",
                            "and the pairs cluster bootstrap are in the source file"),
              source = "06q_cdh_bootstrap.csv")

# ---- 11: carryover, and matched event time ---------------------------------
co <- src("11_fect_carryover.csv")
R$carry <- reg("diagnostic", "FEct", co$spec, "combined >=50%",
               fifelse(co$spec %like% "removed",
                       "estimation sample less the waves following a contraction",
                       "estimation sample"),
               subj(co$subject), co$att, co$se,
               note = "tests whether carryover into the untreated pool explains the FEct/CDH gap",
               source = "11_fect_carryover.csv")

# 11b's matched event time is a dynamic profile, so only its aggregate-free rows
# would fit here; the per-k estimates are the point of that table and stay in it.
# Recorded as an explicit exclusion rather than an omission.
EXCLUDED_BY_DESIGN <- c("11b_dynamic_matched.csv (per-event-time, not an ATT)")

# ---- 10_heterogeneity_results.csv ------------------------------------------
he <- src("10_heterogeneity_results.csv")
R$het <- reg("heterogeneity (exploratory)", "CDH",
             sprintf("%s = %s", he$dimension, he$group), "combined >=50%",
             paste0(he$sample, " subsample"), subj(he$subject), he$att, he$se,
             note = "split-sample subgroup effect, reported as descriptive",
             source = "10_heterogeneity_results.csv")

# ---- Assemble and validate -------------------------------------------------
out <- rbindlist(R, use.names = TRUE)
# Whether the row is estimated with the PRIMARY estimator. Mechanical, not a
# judgement, and it encodes the rule that specification robustness belongs on the
# estimator the headline rests on: a row with on_primary = FALSE answers "what does
# this OTHER estimator depend on", which cannot substitute for the same check on
# the primary. Filter on it to get the checks that bear on the headline.
out[, on_primary := estimator == "CDH"]
setorder(out, role, estimator, spec, outcome)
out[, `:=`(att = round(att, 4), se = round(se, 4))]

if (anyNA(out$outcome)) stop("unrecognised subject label in a source table")
if (anyNA(out$estimator)) stop("unrecognised estimator label in a source table")

# COMPLETENESS. Every committed table carrying an `att`-like column is either read
# above or named as a deliberate exclusion. Asserted rather than trusted: the
# registry claims to index every ATT the project produces, and seven tables were
# missing from it for weeks precisely because nothing checked.
ATT_LIKE  <- c("att", "att_pp", "estimate", "fect_att", "cdh_lang", "cdh_math")
EXCLUDED  <- c(
  # Not ATTs: event studies, per-period profiles, test batteries, bounds.
  "06_placebos.csv", "06g_horizon_sensitivity.csv", "06g_max_placebos.csv",
  "06g_cdh_options.csv", "07_event_study_cdh.csv", "07_event_study_cdh_controls.csv",
  "07_event_study_cs_covariates.csv", "07b_cohort_decomposition.csv",
  "07b_kminus2_by_region.csv", "07c_dynamic_profile.csv",
  "08_event_study_genuine.csv", "09_event_study_cdh.csv",
  "09_event_study_trajectories.csv", "11_fect_dynamic.csv", "11_fect_tests.csv",
  "11c_pretrend_leads.csv", "11d_switchers_by_horizon.csv",
  "06m_exposure_by_horizon.csv", "06n_exposure_estimates.csv",
  "06n_placebo_cancellation.csv", "06n_cohort_diagnostic.csv",
  "06i_size_lpm.csv", "06i_size_within_cohort.csv", "06h_composition_by_group.csv",
  "10_quartile_direction.csv", "06q_bootstrap_draws.csv",
  "06q_cdh_bootstrap_partA.csv",
  # Per-event-time by design; the profile is the point of the table.
  "11b_dynamic_matched.csv",
  # Not an achievement effect: the outcome is the SES covariate itself, so these
  # rows test whether treatment moves school composition. Same numbers, different
  # dependent variable, and `outcome` here only takes the four this schema knows.
  "06h_composition_test.csv",
  # One estimate under several inference procedures, as in 06q. The registry is
  # keyed on (spec, sample, outcome), so it cannot hold the analytic and the
  # resampled rows of the same coefficient without a duplicate key. The point
  # estimates already enter through 06f_clustering.csv; what this table adds is
  # the procedure, which belongs beside the intervals it produces.
  "06f_resampling.csv",
  # This file.
  "92_estimate_registry.csv")
have_att <- Filter(function(f)
  any(ATT_LIKE %in% names(fread(file.path(DIR_TABLES, f), nrows = 0))),
  list.files(DIR_TABLES, pattern = "[.]csv$"))
missed <- setdiff(have_att, c(unique(out$source), EXCLUDED))
if (length(missed))
  stop("tables carrying estimates that the registry neither reads nor excludes: ",
       paste(missed, collapse = ", "),
       ". Add them above, or name them in EXCLUDED with the reason.")
if (out[role == "main", .N] != 2L)
  stop("expected exactly two rows with role 'main' (one per subject), got ",
       out[role == "main", .N])
if (nrow(unique(out[, .(role, estimator, spec, treatment, sample, outcome)])) != nrow(out))
  stop("duplicate specification keys in the registry")

fwrite(out, file.path(DIR_TABLES, "92_estimate_registry.csv"))
