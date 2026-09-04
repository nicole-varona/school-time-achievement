# ==============================================================================
# 07b_pretrend_decomposition.R — Where does the Maths k=-2 placebo come from?
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# The Maths event-study has ONE significant pre-treatment placebo at k=-2
# (~+0.10 SD); Language is clean. Rather than trying to remove it, we CHARACTERISE
# it (per methodological advice): map the k=-2 coefficient back to the adoption
# cohorts / calendar comparisons and the provinces that identify it.
#   (1) Cohort/calendar decomposition: Callaway-Sant'Anna group-time ATT(g,t) at
#       the k=-2 pre-period, by cohort g (each maps to a specific calendar gap
#       because outcome waves are irregularly spaced).
#   (2) By region: the CDH k=-2 placebo estimated within each region, to see
#       whether a single region drives the pooled violation.
# NO province trends / cohort restriction (those would tune the specification).
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(did); library(DIDmultiplegtDYN); library(ggplot2); library(polars)})
# Randomness is seeded PER CALL, immediately before each bootstrap, via
# set_call_seed() (00_setup.R) -- not once here. A script-level seed makes every
# result depend on how much randomness the preceding calls happened to consume,
# so the mathematics fit below would inherit whatever the language fit left
# behind, and reordering the two would move both.

panel <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
# Treatment comes from 02 (the single source of truth) instead of being rebuilt
# here. The persistence rule used to live in nine copies, and shift() silently
# depends on the panel being sorted -- an invariant no script declared.
panel <- merge(panel,
               readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))[
                 , .(school_id, year, main_treatment)],
               by = c("school_id", "year"), all.x = TRUE)

d <- panel[year %in% WAVES & !is.na(score_lang),
           .(school_id, year, y_lang = score_lang / 100, y_math = score_math / 100, main_treatment)]
d <- d[!is.na(main_treatment)]
d[, `:=`(wave = match(year, WAVES), gid = as.integer(factor(school_id)))]

# ---- Province -> region (fix mojibake duplicates, then map) ------------------
# Collapse both encodings (proper UTF-8 and the mojibake duplicates) to a common
# ASCII skeleton by dropping every non-[A-Za-z ] BYTE. Accented vowels vanish, so
# region patterns below match the skeleton (e.g. "Cordoba"/"CArdoba" -> "crdoba").
normalize_prov <- function(x) {
  x <- gsub("[^A-Za-z ]", "", as.character(x), useBytes = TRUE)  # ASCII skeleton
  x <- tolower(trimws(x))
  gsub(" +", " ", x)
}
prov <- panel[!is.na(province), .N, by = .(school_id, province)][order(school_id, -N)]
prov <- prov[, .SD[1L], by = school_id][, .(school_id, pn = normalize_prov(province))]
prov[, region := fcase(
  grepl("jujuy|salta|tucum|catamarca|santiago|rioja", pn),      "NOA",
  grepl("formosa|chaco|corrientes|misiones", pn),               "NEA",
  grepl("mendoza|san juan|san luis", pn),                       "Cuyo",
  grepl("neuqu|negro|chubut|santa cruz|fuego|tierra", pn),      "Patagonia",
  grepl("ciudad|caba|capital|auton", pn),                       "CABA",
  grepl("rdoba|santa fe|entre r|pampa|buenos aires", pn),       "Centro/Pampeana",
  default = "Otra/NA")]
d <- merge(d, prov[, .(school_id, region)], by = "school_id", all.x = TRUE)

# ---- Absorbing subsample + cohort G ------------------------------------------
d_cs <- cs_sample(d)

# ============================================================================
# (1) Cohort / calendar decomposition — Callaway-Sant'Anna group-time ATT(g,t)
# ============================================================================
# With a universal base period, ATT(g, t<g) are the pre-period placebos. The
# pooled k=-2 placebo is a weighted average of ATT(g, g-2) across cohorts g; each
# g maps to a specific calendar comparison because waves are irregularly spaced.
gt_decomp <- function(yname, subj) {
  set_call_seed(paste("gt decomp", subj))
  fit <- tryCatch(att_gt(yname = yname, tname = "wave", idname = "gid", gname = "G",
                         data = d_cs[!is.na(get(yname))], control_group = "notyettreated",
                         allow_unbalanced_panel = TRUE, bstrap = TRUE, biters = BITERS, base_period = "universal", est_method = "dr"),
                  error = function(e) { NULL })
  if (is.null(fit)) return(NULL)
  gt <- data.table(g = fit$group, t = fit$t, att = fit$att, se = fit$se)
  gt[, `:=`(k = t - g, g_year = WAVES[g], t_year = WAVES[t], subject = subj)]
  gt[]
}
gt <- rbindlist(list(gt_decomp("y_lang", "language"), gt_decomp("y_math", "mathematics")), fill = TRUE)
k2 <- gt[k == -2][, .(subject, cohort_year = g_year, placebo_year = t_year,
                      calendar_gap = paste0(t_year, "->", g_year),
                      att = round(att, 4), se = round(se, 4),
                      sig = ifelse(abs(att) > 1.96 * se, "*", ""))]
fwrite(gt, file.path(DIR_TABLES, "07b_cohort_decomposition.csv"))

# ============================================================================
# (2) By-region k=-2 placebo (does ONE region drive it?) — CDH per region
# ============================================================================
cdh_k2 <- function(dat, yname) {
  r <- tryCatch(did_multiplegt_dyn(df = as.data.frame(dat[!is.na(get(yname))]), outcome = yname,
                  group = "gid", time = "wave", treatment = "main_treatment", effects = 2, placebo = 2,
                  cluster = "gid", graph_off = TRUE),
                error = function(e) NULL)
  if (is.null(r) || is.null(r$results$Placebos)) return(c(NA_real_, NA_real_, NA_real_))
  p <- as.data.table(r$results$Placebos)
  p[, k := -as.integer(sub("Placebo_", "", rownames(r$results$Placebos)))]
  row <- p[k == -2]
  if (!nrow(row)) return(c(NA_real_, NA_real_, NA_real_))
  c(row$Estimate, row$Estimate - 1.96 * row$SE, row$Estimate + 1.96 * row$SE)
}
regs <- d[!is.na(region) & region != "Otra/NA", unique(region)]
by_reg <- rbindlist(lapply(regs, function(rg) {
  v <- cdh_k2(d[region == rg], "y_math")
  data.table(region = rg, schools = uniqueN(d[region == rg]$gid),
             k2_math = v[1], lo = v[2], hi = v[3])
}), fill = TRUE)
by_reg_all <- cdh_k2(d, "y_math")   # pooled, for reference
by_reg <- rbind(by_reg, data.table(region = "ALL (pooled)", schools = uniqueN(d$gid),
                                    k2_math = by_reg_all[1], lo = by_reg_all[2], hi = by_reg_all[3]))
by_reg[, `:=`(k2_math = round(k2_math, 4), lo = round(lo, 4), hi = round(hi, 4),
              sig = ifelse(!is.na(lo) & (lo > 0 | hi < 0), "*", ""))]
fwrite(by_reg, file.path(DIR_TABLES, "07b_kminus2_by_region.csv"))

