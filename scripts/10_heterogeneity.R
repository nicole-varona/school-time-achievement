# ==============================================================================
# 10_heterogeneity.R — Heterogeneity as a PRE-REGISTERED POLICY EXTENSION
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Pre-registered design (see §sec-heterogeneidad). PRIMARY (confirmatory) =
# socioeconomic disadvantage, measured by the OFFICIAL 2016 context-vulnerability
# quartile (a baseline classification not derived from the assessed cohort).
# Split-sample: the ATT is estimated WITHIN each quartile with the primary
# estimator (did_multiplegt_dyn, binary >=50% in the year), plus ONE supplementary
# interaction-style test (independent-sample z on the most- vs least-disadvantaged
# quartiles). Directional hypothesis: larger effects for more disadvantaged
# schools (compensatory). Aprender-2016 SES terciles = robustness. SECONDARY
# (descriptive, NOT hypothesis tests): region, sector, urban/rural.
#
# TWO CONSTRAINTS CARRIED FROM THE ESTIMATOR FINDING (11/11b/11c), and they are
# why this is presented as DESCRIPTIVE / EXPLORATORY rather than causal-by-cell:
#   (1) Subgroup ATTs inherit the identification problem of the overall estimate,
#       and worse, because switchers collapse further once the sample is split.
#       The number of Effect_1 switchers is therefore reported for EVERY cell, so
#       "no effect" can be told apart from "no identifying variation".
#   (2) The 2023-cohort DIAGNOSTIC EXCLUSION is repeated WITHIN each subgroup. A
#       subgroup "effect" that vanishes once the 2023 cohort is removed is the
#       same artefact documented in 11c, not genuine heterogeneity. Every cell is
#       reported both on the full sample and excluding the 2023 cohort.
#
# Disciplined-language rule (pre-registered): estimated points are "larger/smaller
# in ...", never "the programme works better in ...".
# Region x identification caveat: the Maths pre-trend lives in the Northeast, so a
# "Northeast effect" for Maths is exactly where identification is weakest; it is
# flagged, reported with the 2023-excluded column, never dropped or suppressed.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
for (p in c("DIDmultiplegtDYN", "ggplot2"))
  if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({ library(polars); library(ggplot2) })
# Randomness: per-call via set_call_seed() (00_setup.R), not a global seed.

# ---- Load the estimation sample ---------------------------------------------
# Panel, treatment definitions and baseline covariates are built ONCE in
# 02_build_estimation_sample.R. Rebuilding them here (as this script used to do)
# is what allowed the same quantity to be measured differently across scripts.
#   d                : school x outcome wave, ready to estimate
#   schools          : one row per school (reversal, adoption cohort, subsamples)
#   treatment_annual : treatment on the full annual census grid
d                <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
schools          <- readRDS(file.path(DIR_PROCESSED, "school_attributes.rds"))
treatment_annual <- readRDS(file.path(DIR_PROCESSED, "treatment_annual.rds"))

# ---- 2023 adoption cohort (for the within-subgroup diagnostic exclusion) -----
ann <- treatment_annual[!is.na(main_treatment)][order(school_id, year)]
coh23_ids  <- schools[cohort_2023 == TRUE, school_id]        # 2023 adopters

# ---- Subgroup variables: one per school, fixed at 2016 baseline --------------
# ASCII-skeleton normalisation of province (the RA province strings carry mojibake
# accents that break literal matching; matching on the accent-stripped skeleton is
# robust). Region = the five standard Argentine regions.
skel <- function(x) gsub("[^A-Za-z ]", "", iconv(as.character(x), to = "ASCII//TRANSLIT"), useBytes = TRUE)
region_of <- function(prov) {
  s <- skel(prov)
  fcase(
    grepl("Jujuy|Salta|Tucum|Catamarca|Rioja|Santiago", s), "NOA",
    grepl("Formosa|Chaco|Corrientes|Misiones", s),          "NEA",
    grepl("Mendoza|San Juan|San Luis", s),                  "Cuyo",
    grepl("Neuqu|Negro|Chubut|Santa Cruz|Fuego", s),        "Patagonia",
    grepl("Buenos Aires|Ciudad|rdoba|Santa Fe|Entre|Pampa", s), "Centro/Pampeana",
    default = NA_character_)
}
# The subgroup variables are NOT part of the estimation sample: the vulnerability
# quartile, the Aprender SES index and the province live in the 2016 wave of the
# unified panel. Reading them here is not a re-derivation of the treatment (which
# comes from 02); it is simply pulling extra school characteristics.
panel <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
base <- panel[year == 2016, .(school_id, context_quartile, ses_index,
                              province, sector, area)]
base[, `:=`(
  region     = region_of(province),
  sector     = factor(trimws(sector)),
  area       = factor(trimws(area)),
  ses_tercile= as.integer(cut(ses_index, breaks = quantile(ses_index, c(0, 1/3, 2/3, 1), na.rm = TRUE),
                              include.lowest = TRUE, labels = FALSE)))]
# `d` already carries sector/area from the baseline covariates (RA-sourced, with a
# fallback year). The subgroup cuts use the 2016 Aprender-frame versions instead,
# so drop the duplicates before merging or the join would produce .x/.y columns
# and `run_dimension("sector")` would fail to find its variable.
d[, c("sector", "area") := NULL]
d <- merge(d, base[, .(school_id, context_quartile, ses_tercile, region, sector, area)],
           by = "school_id", all.x = TRUE)

# ---- Direction check: which context quartile is the most disadvantaged? ------
dir_chk <- panel[year == 2016 & !is.na(context_quartile),
                 .(n = .N, mean_ses = round(mean(ses_index, na.rm = TRUE), 3),
                   mean_lang = round(mean(score_lang, na.rm = TRUE), 1),
                   mean_math = round(mean(score_math, na.rm = TRUE), 1)),
                 by = context_quartile][order(context_quartile)]
fwrite(dir_chk, file.path(DIR_TABLES, "10_quartile_direction.csv"))

# ============================================================================
# Estimation engine: CDH within a subsample, returning ATT + SE + switchers
# ============================================================================
# effects = 2 / placebo = 1 as in 06. The Effect_1 switcher count is the key
# identifying-variation measure per cell (see constraint 1 above). Analytical
# clustered SEs (CDH default for a binary treatment); no bootstrap, so fast.
run_cdh_cell <- function(dat, yname) {
  set_call_seed(paste(yname, nrow(dat)))
  if (uniqueN(dat$gid) < 50 || dat[main_treatment == 1, uniqueN(gid)] < 20)
    return(data.table(att = NA_real_, se = NA_real_, sw1 = NA_integer_, n_sch = uniqueN(dat$gid)))
  out <- tryCatch(DIDmultiplegtDYN::did_multiplegt_dyn(
    df = as.data.frame(dat), outcome = yname, group = "gid", time = "wave",
    treatment = "main_treatment", effects = 2, placebo = 1, cluster = "gid",
    graph_off = TRUE), error = function(e) NULL)
  if (is.null(out)) return(data.table(att = NA_real_, se = NA_real_,
                                      sw1 = NA_integer_, n_sch = uniqueN(dat$gid)))
  eff <- as.data.frame(out$results$Effects)
  data.table(att = out$results$ATE[1], se = out$results$ATE[2],
             sw1 = as.integer(eff[["Switchers"]][1]), n_sch = uniqueN(dat$gid))
}

# One dimension: loop its groups x {full, excl-2023} x {language, maths}.
run_dimension <- function(dimvar, dim_label, exclude_2023 = TRUE) {
  grps <- sort(unique(na.omit(d[[dimvar]])))
  samples <- if (exclude_2023) c("full", "excl. 2023 cohort") else "full"
  out <- rbindlist(lapply(grps, function(g) {
    rbindlist(lapply(samples, function(sm) {
      dd <- d[get(dimvar) == g]
      if (sm == "excl. 2023 cohort") dd <- dd[!school_id %in% coh23_ids]
      rbindlist(lapply(c("y_lang", "y_math"), function(yn) {
        r <- run_cdh_cell(dd, yn)
        cbind(data.table(dimension = dim_label, group = as.character(g), sample = sm,
                         subject = if (yn == "y_lang") "language" else "mathematics"), r)
      }))
    }))
  }))
  out
}

res_quart  <- run_dimension("context_quartile", "context quartile")
res_ses    <- run_dimension("ses_tercile", "SES tercile")
res_region <- run_dimension("region", "region")
res_sector <- run_dimension("sector", "sector", exclude_2023 = FALSE)
res_area   <- run_dimension("area", "area", exclude_2023 = FALSE)

het <- rbindlist(list(res_quart, res_ses, res_region, res_sector, res_area), fill = TRUE)
het[, `:=`(att = round(att, 4), se = round(se, 4))]
fwrite(het, file.path(DIR_TABLES, "10_heterogeneity_results.csv"))

# ============================================================================
# One supplementary interaction test: most- vs least-disadvantaged quartile
# ============================================================================
# Independent-sample z on the split-sample ATTs (the two subsamples are disjoint,
# so their estimates are independent). Reported on the full sample AND excluding
# the 2023 cohort, since constraint (2) applies to the contrast too.
z_test <- function(tab, g_lo, g_hi, sm) {
  a <- tab[group == g_lo & sample == sm]; b <- tab[group == g_hi & sample == sm]
  rbindlist(lapply(c("language", "mathematics"), function(su) {
    x <- a[subject == su]; y <- b[subject == su]
    if (nrow(x) == 0 || nrow(y) == 0 || is.na(x$att) || is.na(y$att))
      return(data.table(subject = su, sample = sm, diff = NA_real_, z = NA_real_, p = NA_real_))
    diff <- x$att - y$att; se <- sqrt(x$se^2 + y$se^2); z <- diff / se
    data.table(subject = su, sample = sm, diff = round(diff, 4),
               z = round(z, 3), p = round(2 * pnorm(-abs(z)), 4))
  }))
}
q_lo <- min(res_quart$group); q_hi <- max(res_quart$group)   # extreme quartiles
inter <- rbindlist(lapply(c("full", "excl. 2023 cohort"),
                          function(sm) z_test(res_quart, q_lo, q_hi, sm)), fill = TRUE)
inter[, `:=`(comparison = sprintf("Q%s vs Q%s", q_lo, q_hi))]

# The same contrast on the ROBUSTNESS measure (Aprender SES terciles). Reported
# because the two disadvantage measures behave differently under the diagnostic
# exclusion -- the official quartile contrast weakens, the tercile one does not --
# and that comparison cannot be made if only the quartile test is written out.
t_lo <- min(res_ses$group); t_hi <- max(res_ses$group)
inter_ses <- rbindlist(lapply(c("full", "excl. 2023 cohort"),
                              function(sm) z_test(res_ses, t_lo, t_hi, sm)), fill = TRUE)
inter_ses[, `:=`(comparison = sprintf("T%s vs T%s (SES tercile)", t_lo, t_hi))]

fwrite(rbind(inter, inter_ses), file.path(DIR_TABLES, "10_interaction_test.csv"))

# Figure 10_heterogeneity.png is rendered from 10_heterogeneity_results.csv by
# 91_figures.R in the shared style (points + 95% CI, faceted by dimension).

# ---- Switcher-count summary: where is there identifying variation at all? -----
sw_summary <- het[sample == "full" & subject == "language",
                  .(dimension, group, n_schools = n_sch, switchers_eff1 = sw1)]
fwrite(sw_summary, file.path(DIR_TABLES, "10_switchers_per_cell.csv"))

# Fitted model objects are NOT saved. They were, and nothing ever read them:
# 855 MB of caches sat in data/processed alongside the three real tables and
# made the data flow hard to follow. Every number this script produces is in
# its CSV outputs, and the registry (92) indexes them.
