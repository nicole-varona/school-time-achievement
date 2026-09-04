# ==============================================================================
# 06h_composition_check.R — Does the treatment change who sits the test?
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/estimation_sample.rds, panel.rds
# Outputs : output/tables/06h_composition_by_group.csv   raw and by wave pair
#           output/tables/06h_composition_test.csv       the contrast, with FE
#           output/tables/06h_ses_covariate.csv          predetermined vs 2016 SES
#
# The design needs the treatment not to change the socioeconomic composition of
# the school. If a longer day pulls different families in or out, then any
# contemporaneous socioeconomic measure sits DOWNSTREAM of the treatment, and the
# repeated-cross-section outcome is partly measuring a different population
# rather than a different level of learning. Maternal education is observed in all
# six waves, so the channel is directly testable.
#
# The test compares the wave-to-wave CHANGE in the share of mothers with at least
# secondary education across schools grouped by what their treatment just did.
# Two things make the raw comparison misleading and both are handled:
#   - maternal education rises steeply over time for everyone, so only the
#     DIFFERENTIAL change is informative;
#   - adoption is concentrated in particular wave pairs (above all 2021 to 2023),
#     so a raw contrast confounds "just adopted" with "measured over a wave pair
#     in which everyone rose". The contrast is therefore run WITH wave-pair fixed
#     effects, which is what identifies it, and the by-pair table is reported so
#     the reader can see the same thing without the regression.
#
# Section 3 additionally compares the two ways of building the SES covariate: the
# predetermined measure now in use, built in 02 from the waves preceding each
# school's own adoption, and the 2016 value it replaced. Same estimator, same
# subsample rule, so the difference is the covariate and nothing else.
# ============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))
suppressPackageStartupMessages({library(did); library(fixest)})

d     <- readRDS(file.path(DIR_PROCESSED, "estimation_sample.rds"))
panel <- readRDS(file.path(DIR_PROCESSED, "panel.rds"))
if (!"bl_mom_pre" %in% names(d))
  stop("bl_mom_pre missing from the estimation sample; run 02 first")

# ---- Build the wave-to-wave change and the switching groups ------------------
d <- merge(d, panel[, .(school_id, year, mom = pct_mother_secondary)],
           by = c("school_id", "year"), all.x = TRUE)
setorder(d, gid, wave)
d[, `:=`(dD    = main_treatment - shift(main_treatment),
         dmom  = mom - shift(mom),
         lagT  = shift(main_treatment),
         wpair = paste0(shift(year), "-", year)), by = gid]
# Re-adoption: an increase that follows an earlier decrease in the same school.
d[, prior_exit := cumsum(fifelse(!is.na(dD) & dD < 0, 1L, 0L)) -
                  fifelse(!is.na(dD) & dD < 0, 1L, 0L), by = gid]
d[, grp := fifelse(is.na(dD), NA_character_,
            fifelse(dD > 0 & prior_exit > 0, "re-adopts",
            fifelse(dD > 0, "just adopts",
            fifelse(dD < 0, "just exits",
            fifelse(lagT == 0, "still untreated", "still treated")))))]

x <- d[!is.na(dmom) & !is.na(grp)]

# ---- 1. The change by group, pooled and by wave pair ------------------------
pooled <- x[, .(wave_pair = "all (pooled)", n = .N,
                mean_change = mean(dmom), se = sd(dmom) / sqrt(.N)), by = grp]
bypair <- x[, .(n = .N, mean_change = mean(dmom), se = sd(dmom) / sqrt(.N)),
            by = .(grp, wave_pair = wpair)]
comp <- rbind(pooled, bypair)[, `:=`(mean_change = round(mean_change, 4),
                                     se = round(se, 4))]
setcolorder(comp, c("wave_pair", "grp", "n", "mean_change", "se"))
fwrite(comp[order(wave_pair, grp)], file.path(DIR_TABLES, "06h_composition_by_group.csv"))

# ---- 2. The contrast that identifies it -------------------------------------
# Each switching group against schools that did not change treatment and were
# untreated, with and without wave-pair fixed effects. Without them the estimate
# is the raw difference in the table above; with them it is the differential
# change within the same wave pair, which is the quantity of interest.
x[, grp := relevel(factor(grp), ref = "still untreated")]
fits <- list(`no wave-pair FE` = feols(dmom ~ grp, data = x, cluster = ~gid),
             `wave-pair FE`    = feols(dmom ~ grp | wpair, data = x, cluster = ~gid))
test <- rbindlist(lapply(names(fits), function(nm) {
  m <- fits[[nm]]; co <- coef(m); se <- se(m)
  k <- grep("^grp", names(co), value = TRUE)
  data.table(spec = nm, group = sub("^grp", "", k),
             estimate = round(co[k], 4), se = round(se[k], 4),
             t = round(co[k] / se[k], 2))
}))
# The scale the estimate has to be read against: dispersion BETWEEN schools.
test[, sd_between_schools := round(sd(d$mom, na.rm = TRUE), 4)]
fwrite(test, file.path(DIR_TABLES, "06h_composition_test.csv"))

# ---- 3. The two SES covariates, same estimator and same subsample -----------
# The absorbing subsample, built exactly as in 06_estimation.R.
d_cs <- cs_sample(d)

run_cs <- function(yname, xformla, label) {
  set_call_seed(label)
  cv <- all.vars(xformla)
  dd <- d_cs[!is.na(get(yname))][complete.cases(d_cs[!is.na(get(yname)), ..cv])]
  fit <- tryCatch(att_gt(yname = yname, tname = "wave", idname = "gid", gname = "G",
                    xformla = xformla, data = dd, control_group = "notyettreated",
                    allow_unbalanced_panel = TRUE, bstrap = TRUE, cband = FALSE,
                    biters = BITERS, est_method = "dr"), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  a <- suppressWarnings(aggte(fit, type = "simple", na.rm = TRUE))
  data.table(att = round(a$overall.att, 4), se = round(a$overall.se, 4),
             n_schools = uniqueN(dd$gid))
}

SETS <- list(`predetermined (pre-adoption, in use)` = XF_COMP,
             `2016 wave (previous construction)` = XF_COMP_2016)
ses <- rbindlist(lapply(names(SETS), function(nm) {
  rbindlist(lapply(c("y_lang", "y_math"), function(y) {
    r <- run_cs(y, SETS[[nm]], paste(nm, y))
    if (is.null(r)) return(NULL)
    r[, `:=`(ses_covariate = nm,
             subject = fifelse(y == "y_lang", "language", "mathematics"))]
  }), fill = TRUE)
}), fill = TRUE)

# Coverage is part of the comparison: the predetermined measure is missing for
# schools with no pre-treatment wave, and that cost belongs next to the estimate.
cov_tab <- data.table(
  quantity = c("schools with the 2016 SES value", "schools with the predetermined value",
               "schools with no pre-treatment wave", "median pre-treatment waves used"),
  value = c(d[!is.na(bl_pct_mother_sec), uniqueN(gid)],
            d[!is.na(bl_mom_pre), uniqueN(gid)],
            d[is.na(bl_mom_pre), uniqueN(gid)],
            as.numeric(median(unique(d[, .(gid, n_pre_waves)])$n_pre_waves, na.rm = TRUE))))
fwrite(rbind(ses, cov_tab, fill = TRUE), file.path(DIR_TABLES, "06h_ses_covariate.csv"))
