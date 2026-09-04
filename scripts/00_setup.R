# ==============================================================================
# 00_setup.R — Configuration, paths, and helper functions
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Data (not included in the repository; see README for access):
#   - Relevamiento Anual (RA) 2011-2025: annual census of all schools (CSV,
#     UTF-8 BOM, ';' separator, ',' decimal). One folder per year inside
#     data/bases_ra_*.
#   - Aprender national assessment, grade 6 primary: waves 2016, 2018, 2021,
#     2022 (sample), 2023, 2025. Student-level microdata (.sav for 2016-2023,
#     .csv for 2025) in data/aprender_<year>/.
#   - School identifiers are 15-digit fictitious IDs assigned by the
#     Secretariat of Education; consistent across RA and Aprender.
# ==============================================================================

# ---- Packages ----------------------------------------------------------------
required_pkgs <- c("data.table", "haven", "ggplot2", "here")
to_install <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                    logical(1), quietly = TRUE)]
if (length(to_install)) install.packages(to_install,
                                         repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({
  library(data.table)
  library(haven)
  library(ggplot2)
})

# ---- Project paths (anchored to the project root via here) -------------------
# here::here() locates the project root by its markers (.here / .Rproj), so
# scripts run from any working directory.
DIR_DATA      <- here::here("data")
DIR_PROCESSED <- here::here("data", "processed")
DIR_OUTPUT    <- here::here("output")
DIR_TABLES    <- here::here("output", "tables")
DIR_FIGURES   <- here::here("output", "figures")
for (d in c(DIR_PROCESSED, DIR_TABLES, DIR_FIGURES))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- Global parameters -------------------------------------------------------
RA_YEARS       <- 2011:2025      # RA waves available (always a census)
APRENDER_WAVES <- c(2016, 2018, 2021, 2022, 2023, 2025)  # grade-6 Lang/Math
WAVES          <- sort(APRENDER_WAVES)  # outcome waves sorted; scripts map a year
                                        # to its 1..6 index via match(year, WAVES)
WAVE_TYPE      <- c("2016" = "census", "2018" = "census", "2021" = "census",
                    "2022" = "sample", "2023" = "census", "2025" = "census")
ID_WIDTH       <- 15             # fictitious school ID: 15 digits, character
BITERS         <- 1000           # multiplier-bootstrap replications for att_gt.
                                 # ONE constant for the whole project: this used to be a
                                 # local choice and drifted to four different values across
                                 # ten scripts (200, 300, 500), so the SAME specification
                                 # reported different standard errors depending on which
                                 # script produced it. `biters` controls only the Monte Carlo
                                 # error of the bootstrap, not the estimand, so no value is
                                 # more correct -- but they must all be the SAME, and 1000 is
                                 # the package default, which needs no defending.
NBOOTS         <- as.integer(Sys.getenv("NBOOTS", "1000"))
                                 # bootstrap replications for fect standard errors. Same
                                 # argument as BITERS, and set to the same value so the two
                                 # bootstraps in the project are not silently asymmetric.
                                 # Overridable from the shell (NBOOTS=5 Rscript ...) to
                                 # smoke-test the pipeline cheaply.

# ---- RA column patterns ------------------------------------------------------
# RA column headers change across waves. I resolve columns by ASCII regex
# patterns on make.names()-transformed headers, which is robust to
# both naming conventions and to accented characters.
RA_PATTERNS <- list(
  extended_day          = "jornada.*extendida.*primaria",
  full_day              = "jornada.*completa.*primaria",
  double_shift          = "turno.*doble.*primaria",
  indigenous            = "poblaci.*indigena.*primaria",
  computer              = "tenencia.*computadora.*primaria",
  disability            = "discapacidad.*primaria",     # 3 cols (summed)
  meals                 = "aliment",                    # 5 meal counts (summed)
  lunch                 = "almuerzo",                   # lunch alone: the marker of a
                                                        # day that crosses midday. ASCII
                                                        # pattern, no accent to break on.
  enrolment_grades_1_7  = "^x[._]+[1-7]$",
  enrolment_male_1_7    = "^v[._]+[1-7]$",
  repeaters_g6          = "^r[._]+6$",
  overage_g6            = "^s[._]+6$",
  promoted_g6           = "^promovidos[._]+6$",
  promoted_ex_g6        = "^promovidos_ex[._]+6$",
  not_promoted_g6       = "^nopromo[._]+6$",
  left_transfer_g6      = "^scp[._]+6$",          # left with transfer pass
  dropout_g6            = "^ssp[._]+6$"           # left without pass
)

# RA base file name patterns (names are stable even when base numbers change).
RA_FILES <- list(
  population = "poblaci",                # Base 6 (treatment source)
  enrolment  = "matr",                   # Base 2 (exclude "edad" = Base 7)
  trajectory = "trayect"                 # Base 3 (secondary outcome)
)

# ---- Aprender per-wave item map ----------------------------------------------
# Questionnaire item codes change across waves (verified in the official
# dictionaries). NA = the concept is not measured in that wave.
# `format`: "sav" = labelled SPSS file (aggregate via value labels);
#           "csv_coded" = plain CSV with numeric codes (2025).
APRENDER_MAP <- list(
  "2016" = list(format = "sav", age = "Ap1", sex = "Ap2", preschool = "Ap11",
                repeated = "Ap12", absence = "Ap13", mother_ed = "Ap7",
                father_ed = "Ap8", ses = "isocioa", climate = "iclima",
                context_q = "qvulneraa"),
  "2018" = list(format = "sav", age = "ap1", sex = "ap2", preschool = "ap16",
                repeated = "ap17", absence = NA, mother_ed = "ap9",
                father_ed = "ap10", ses = "isocioa", climate = NA,
                context_q = NA),
  "2021" = list(format = "sav", age = "ap01", sex = "ap03", preschool = "ap24",
                repeated = "ap25", absence = "ap30", mother_ed = "Nivel_Ed_Madre",
                father_ed = "Nivel_Ed_Padre", ses = "NSE_puntaje", climate = NA,
                context_q = NA),
  "2022" = list(format = "sav", age = "ap01", sex = "ap03", preschool = "ap22",
                repeated = "ap23", absence = "ap24", mother_ed = "Nivel_Ed_Madre",
                father_ed = "Nivel_Ed_Padre", ses = "NSE_puntaje", climate = NA,
                context_q = NA),
  "2023" = list(format = "sav", age = "ap01", sex = "ap03", preschool = "ap21",
                repeated = "ap22", absence = "ap23", mother_ed = "ap20a1",
                father_ed = "ap20a2", ses = NA, climate = "clima_escolar",
                context_q = NA),
  # 2025 ships as coded CSV (no value labels). Composition items are computed
  # from standard Aprender codings and their distributions are printed in the
  # build log for verification. mother/father education use the constructed
  # ordinal variables; SES is a quintile ("Q1".."Q5") -> numeric 1-5.
  "2025" = list(format = "csv_coded", age = "ap01", sex = "ap03",
                preschool = "ap06", repeated = "ap07", absence = NA,
                mother_ed = "Nivel_Ed_Madre", father_ed = "Nivel_Ed_Padre",
                ses = "NSE_nivel", climate = NA, context_q = NA)
)

# Stable Aprender outcome columns (all waves)
APRENDER_OUTCOME <- list(score_lang = "lpuntaje", score_math = "mpuntaje",
                         level_lang = "ldesemp",  level_math = "mdesemp")

# ==============================================================================
# Helper functions
# ==============================================================================

# ---- Numeric conversion respecting ',' decimals ------------------------------
to_num <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.numeric(x)) return(x)
  suppressWarnings(as.numeric(gsub(",", ".", trimws(x), fixed = TRUE)))
}

# ---- School ID: robust character conversion + zero-padding -------------------
# IDs are 15-digit; a handful lose a leading zero upstream (~6/19k verified).
# sprintf("%.0f") avoids scientific notation when IDs arrive as doubles.
pad_id <- function(x, width = ID_WIDTH) {
  if (!is.character(x)) {
    out <- rep(NA_character_, length(x)); ok <- !is.na(x)
    out[ok] <- sprintf("%.0f", as.double(x[ok]))
    x <- out
  }
  x <- trimws(x)
  short <- !is.na(x) & nchar(x) < width
  x[short] <- paste0(strrep("0", width - nchar(x[short])), x[short])
  x
}

# The school ID column name is not stable across releases (ID, ID1, cueanexo).
ID_CANDIDATES <- c("ID", "ID1", "id", "cueanexo", "CUEANEXO", "cue", "CUE")

standardize_id <- function(dt, candidates = ID_CANDIDATES, to = "school_id") {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) return(dt)
  if (hit[1] != to) setnames(dt, hit[1], to)
  dt[[to]] <- pad_id(dt[[to]])
  dt
}

# ---- File discovery ------------------------------------------------------------
# RA year folders live in three containers with different layouts.
ra_year_dir <- function(year) {
  candidates <- c(
    file.path(DIR_DATA, "bases_ra_2011-2017", paste("Bases Usuarias", year)),
    file.path(DIR_DATA, "bases_ra_2018-2024", paste("Bases Usuarias", year)),
    file.path(DIR_DATA, sprintf("bases_ra_%d", year))
  )
  hit <- candidates[dir.exists(candidates)]
  if (length(hit)) hit[1] else NA_character_
}

aprender_dir <- function(year) {
  d <- file.path(DIR_DATA, sprintf("aprender_%d", year))
  if (dir.exists(d)) d else NA_character_
}

# Find one file whose basename matches ALL `include` patterns (any order) and
# none of `exclude`. When several match, files not matching `prefer_not` are
# preferred (e.g. prefer census over sample releases when both exist).
find_file <- function(dir, include, exclude = NULL, prefer_not = NULL) {
  if (is.na(dir)) return(NA_character_)
  files <- list.files(dir, full.names = TRUE, recursive = FALSE)
  keep <- files
  for (p in include) keep <- keep[grepl(p, basename(keep), ignore.case = TRUE)]
  for (p in exclude %||% character(0))
    keep <- keep[!grepl(p, basename(keep), ignore.case = TRUE)]
  if (length(keep) > 1 && !is.null(prefer_not)) {
    alt <- keep[!grepl(prefer_not, basename(keep), ignore.case = TRUE)]
    if (length(alt)) keep <- alt
  }
  if (length(keep) == 0) return(NA_character_)
  keep[1]
}
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Readers -----------------------------------------------------------------
# RA CSVs: read everything as character (protects the 15-digit ID and defers
# numeric parsing), normalise headers with make.names(), standardise the ID.
read_ra <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  dt <- fread(path, sep = ";", encoding = "UTF-8", colClasses = "character",
              showProgress = FALSE)
  setnames(dt, make.names(names(dt), unique = TRUE))
  standardize_id(dt)[]
}

# Aprender: .sav via haven (keeps value labels); .csv via fread (';', dec ',').
# For CSVs the ID column is forced to character at read time (a 15-digit ID
# parsed as integer64/double would risk precision and formatting issues).
read_aprender <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NULL)
  ext <- tolower(tools::file_ext(path))
  dt <- if (ext == "sav") {
    as.data.table(haven::read_sav(path))
  } else {
    hdr <- names(fread(path, nrows = 0, sep = ";"))
    idc <- intersect(ID_CANDIDATES, hdr)
    fread(path, sep = ";", dec = ",", encoding = "UTF-8", showProgress = FALSE,
          integer64 = "double",
          colClasses = if (length(idc)) setNames(rep("character", length(idc)), idc))
  }
  standardize_id(dt)[]
}

# Read selected labelled items from a wave's DIRECTOR base, as character labels.
# `items` is a named vector, e.g. c(dir_hours_g6 = "dp47"). Returns a data.table
# with school_id (+ year) and the requested items; NULL if base/ID not found.
read_director_items <- function(year, items) {
  f <- find_file(aprender_dir(year), include = c("primaria", "directivo"),
                 exclude = "secundaria", prefer_not = "muestral")
  if (is.na(f)) return(NULL)
  raw <- haven::read_sav(f)
  idc <- intersect(ID_CANDIDATES, names(raw))
  if (length(idc) == 0) return(NULL)
  out <- data.table(school_id = pad_id(as.character(raw[[idc[1]]])), year = year)
  for (nm in names(items)) {
    col <- items[[nm]]
    if (col %in% names(raw))
      out[, (nm) := as.character(haven::as_factor(raw[[col]]))]
  }
  out
}

# ---- Column resolution -------------------------------------------------------
# Case-insensitive exact name lookup (Aprender item codes vary in case).
get_col <- function(dt, name) {
  if (length(name) != 1 || is.na(name)) return(NULL)
  hit <- names(dt)[tolower(names(dt)) == tolower(name)]
  if (length(hit) == 0) return(NULL)
  dt[[hit[1]]]
}

# All column names matching an ASCII regex (robust to accents/renamings).
grep_cols <- function(dt, pattern)
  grep(pattern, names(dt), ignore.case = TRUE, value = TRUE)

# Row-wise sum over the columns matching `pattern` (NA if none match).
sum_matched <- function(dt, pattern) {
  cols <- grep_cols(dt, pattern)
  if (length(cols) == 0) return(rep(NA_real_, nrow(dt)))
  if (length(cols) == 1) return(to_num(dt[[cols]]))
  m <- vapply(cols, function(cl) to_num(dt[[cl]]), numeric(nrow(dt)))
  rowSums(m, na.rm = TRUE)
}

# First column matching `pattern`, as numeric (NULL if none).
num_matched <- function(dt, pattern) {
  cols <- grep_cols(dt, pattern)
  if (length(cols) == 0) return(NULL)
  to_num(dt[[cols[1]]])
}

# Presence flag from a marker column: RA uses "X"/"" markers for some fields
# (e.g. "Turno Doble"). Returns 1 for non-empty/non-zero, 0 otherwise.
flag_matched <- function(dt, pattern) {
  cols <- grep_cols(dt, pattern)
  if (length(cols) == 0) return(NULL)
  v <- trimws(as.character(dt[[cols[1]]]))
  fifelse(is.na(v), NA_integer_, as.integer(v != "" & v != "0"))
}

# ---- Label helpers (SPSS files) ----------------------------------------------
# Value labels as character; falls back to raw values for unlabelled columns.
labels_of <- function(x) {
  if (is.null(x)) return(NULL)
  tryCatch(as.character(haven::as_factor(x)), error = function(e) as.character(x))
}

mean_safe <- function(x) if (is.null(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# ---- Director-reported grade-6 daily hours -----------------------------------
# Extracts the leading integer from the hours item (dp23/dp18/dp47, waves
# 2021/2022/2023). Labels differ across waves, so the bottom category is mapped
# explicitly. Lives here because 08 and 08b both classify on it.
parse_hours <- function(x) {
  x <- as.character(x)
  n <- suppressWarnings(as.integer(sub(".*?([0-9]+).*", "\\1", x)))
  n[grepl("menos de 3|3 horas o menos", tolower(x))] <- 2L
  n
}

# ---- Covariate sets for the doubly-robust specifications ---------------------
# Defined ONCE here because five scripts (06, 07, 08, 09, 11) condition on the
# same vector; keeping a private copy in each is how the treatment definition
# silently diverged before (see 02_build_estimation_sample.R).
#
# XF_STRUCT : RA structural covariates only. Near-complete, so it fits on the
#             full sample and isolates selection into having an Aprender 2016
#             record.
# XF_COMP   : MAIN. Adds Aprender-2016 composition and the school-level SES
#             proxy. Mother's secondary education is in the main set because
#             socioeconomic disadvantage is not merely a predictor of the
#             outcome but an explicit school-SELECTION criterion of the
#             provincial extended-day programmes (Veleda 2013 lists "social
#             vulnerability" as the selection rule in CABA, Buenos Aires,
#             Cordoba, Mendoza and Rio Negro), which makes it a confounder
#             rather than a precision covariate.
# XF_COMP_NOSES : the previous main set, retained so the SES covariate's
#             contribution is reportable rather than assumed.
#
# These enter CS only. The primary CDH estimator differences within school, so
# a time-invariant covariate has delta-X = 0 and is absorbed by the design.
# The SES covariate is PREDETERMINED: mother's education averaged over the waves
# that precede each school's OWN first adoption (`bl_mom_pre`, built in 02), not
# the 2016 value. The 2016 value is pre-treatment for most of the sample but not
# for all of it, and dating the covariate by each school's own adoption makes that
# contamination visible (schools with no pre-treatment wave get NA) instead of
# absorbing it. XF_COMP_2016 keeps the previous construction so the choice is
# reportable rather than assumed; 06h estimates both.
XF_STRUCT     <- ~ sector + area + log_enroll
XF_COMP_NOSES <- ~ sector + area + log_enroll + bl_pct_female + bl_pct_repeater + bl_pct_preschool
XF_COMP       <- update(XF_COMP_NOSES, ~ . + bl_mom_pre)
XF_COMP_2016  <- update(XF_COMP_NOSES, ~ . + bl_pct_mother_sec)

# ---- Baseline covariates (school-level, time-invariant) ----------------------
# Doubly-robust DiD (Callaway-Sant'Anna) conditions on PRE-TREATMENT covariates:
# the parallel-trends assumption is required to hold only conditional on X, and X
# must be measured before treatment so it is not itself an outcome. We therefore
# build ONE time-invariant covariate vector per school, fixed at `base_year`
# (default 2016, the first outcome wave), and merge it onto every school-wave.
#
# Structural covariates (sector, area, enrolment) come from the RA census, which
# covers ~all schools every year; when a school is absent from RA `base_year` we
# fall back to its earliest available RA record. Composition covariates come from
# Aprender `base_year`, so schools first tested after `base_year` have them
# missing (they still carry the structural covariates). ses_index is kept on the
# single `base_year` scale only (it is not comparable across waves; see build log).
#
# Returns a data.table keyed by school_id with:
#   sector, area          (factors)      log_enroll = log(max(enroll_primary, 1))
#   bl_pct_female, bl_pct_repeater, bl_pct_preschool, bl_ses  (Aprender base_year)
#   bl_struct_complete, bl_comp_complete (completeness flags for the two specs)
make_baseline_covars <- function(panel, ra_panel, base_year = 2016L) {
  ids <- unique(panel$school_id)

  # Structural: RA base_year, else earliest RA year with valid enrolment.
  ra <- ra_panel[!is.na(enroll_primary),
                 .(school_id, year, sector, area, enroll_primary)]
  base_ra <- ra[year == base_year]
  miss    <- setdiff(ids, base_ra$school_id)
  fb      <- ra[school_id %in% miss][order(school_id, year), .SD[1L], by = school_id]
  struct  <- unique(rbind(base_ra, fb, fill = TRUE)[
               , .(school_id, sector, area, enroll_primary)], by = "school_id")

  # Composition: Aprender base_year only (single, comparable ses scale).
  # bl_score_* is the base-year school mean score (SD units): NOT used in the main
  # spec, only in a "quasi value-added" robustness that conditions on it.
  comp <- panel[year == base_year, .(school_id,
                bl_pct_female    = pct_female,
                bl_pct_repeater  = pct_repeater,
                bl_pct_preschool = pct_preschool,
                bl_ses           = ses_index,
                bl_pct_mother_sec = pct_mother_secondary,
                bl_score_lang    = score_lang / 100,
                bl_score_math    = score_math / 100)]

  # Grade-6 enrolment, same base_year-with-fallback rule as the structural block.
  # It is the weight for the pupil-average estimand: the outcome is a grade-6
  # mean, so total primary enrolment weights by the wrong population.
  g6      <- ra_panel[!is.na(enroll_g6) & enroll_g6 > 0, .(school_id, year, enroll_g6)]
  base_g6 <- g6[year == base_year]
  miss_g6 <- setdiff(ids, base_g6$school_id)
  fb_g6   <- g6[school_id %in% miss_g6][order(school_id, year), .SD[1L], by = school_id]
  struct_g6 <- unique(rbind(base_g6, fb_g6)[, .(school_id, enroll_g6)], by = "school_id")

  out <- merge(data.table(school_id = ids), struct, by = "school_id", all.x = TRUE)
  out <- merge(out, struct_g6, by = "school_id", all.x = TRUE)
  out <- merge(out, comp, by = "school_id", all.x = TRUE)
  blank_na <- function(x) { x <- trimws(as.character(x)); x[x == ""] <- NA; x }
  out[, `:=`(log_enroll = log(pmax(enroll_primary, 1)),
             sector     = factor(blank_na(sector)),
             area       = factor(blank_na(area)))]
  out[, bl_struct_complete := !is.na(sector) & !is.na(area) & !is.na(enroll_primary)]
  out[, bl_comp_complete   := bl_struct_complete & !is.na(bl_pct_female) &
                              !is.na(bl_pct_repeater) & !is.na(bl_pct_preschool)]
  out[]
}

# ---- The absorbing subsample Callaway-Sant'Anna requires ----------------------
# CS assumes adoption is absorbing, and this design is not, so every CS estimate
# in the project runs on the same restricted frame: schools that never contract,
# excluding those already treated in their first observed wave (no pre-period),
# with the cohort variable `G` the estimator needs. That derivation used to be
# copied inline into TWELVE scripts (06, 06c, 06f, 06h, 07, 07b, 07c, 08, 09, 11,
# 11b, 11d). It is the SAMPLE DEFINITION of every CS number reported, so a change
# in one copy would have left the other eleven estimating on a different
# population with no error raised -- the same failure mode that made 02 the single
# source of the treatment definitions.
#
# `drop_reverters = FALSE` is for callers whose input is already non-reverting
# (11b), where the filter is a no-op but stating it is clearer than relying on it.
# The row order is set here rather than assumed: `reverts` is computed with
# diff(), which reads the physical order of the rows, and an undeclared ordering
# assumption is how six scripts once depended on 01 happening to sort its output.
cs_sample <- function(dat, tv = "main_treatment", drop_reverters = TRUE) {
  x <- copy(dat)
  setorder(x, gid, wave)
  g <- x[, .(reverts = .N >= 2 && any(diff(get(tv)) < 0),
             ever    = any(get(tv) == 1),
             firstG  = if (any(get(tv) == 1)) min(wave[get(tv) == 1]) else 0L),
         by = gid]
  if (drop_reverters) g <- g[reverts == FALSE]
  merge(dat, g[!(ever == TRUE & firstG == 1), .(gid, G = as.numeric(firstG))],
        by = "gid")
}

# ---- How many placebos a design can carry ------------------------------------
# Bounded twice: by the periods left before the switch, and by the package's own
# rule that there be no more placebos than effects. Computed rather than
# discovered, because did_multiplegt_dyn truncates an over-request with a message
# instead of failing, and a truncation nobody reads is how a table ends up
# claiming more pre-treatment evidence than it has.
max_placebo <- function(effects, n_periods)
  max(min(effects, n_periods - effects - 1L), 0L)

# ---- Reproducible randomness -------------------------------------------------
# Bootstrapped standard errors are reproducible only if the RNG state at the
# moment of each call is fixed. Seeding once at the top of a script does NOT
# achieve that: every result then depends on how much randomness the preceding
# calls happened to consume, so adding, removing or reordering a specification
# silently shifts the standard errors of all the ones that follow. (Observed in
# practice: refactoring moved bootstrap SEs by up to 0.003 SD while leaving every
# point estimate identical.)
#
# set_call_seed() is therefore called immediately BEFORE each stochastic
# estimation. Each labelled call gets its own deterministic stream derived from
# the project seed, so a result depends only on its own label and never on its
# position in the script.
SEED <- 20260720L

set_call_seed <- function(label = NULL) {
  s <- if (is.null(label) || !nzchar(label)) SEED
       else SEED + sum(utf8ToInt(label))
  set.seed(s)
  invisible(s)
}
