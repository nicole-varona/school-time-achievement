# ==============================================================================
# 08b_recording_rates.R — How well the census records each length of school day
#
# Extended School Time and Student Achievement in Argentine Primary Education
# Dissertation — MSc in Social Research Methods (LSE, Department of Methodology)
# Author  : Candidate 62628
# ------------------------------------------------------------------------------
# Inputs  : data/processed/panel.rds
# Outputs : output/tables/08b_recording_by_hours.csv
#           output/tables/08b_recording_by_province.csv
#           output/tables/08b_five_hour_margin.csv
#
# 08 measures how much 1h+ contaminates the treated group. This measures the
# instrument itself: taking the directors' reported grade-6 hours as the length
# of the day, how often does the census record a school of that length as
# extended day? The recording rate rises with reported hours, so the census is a
# good measure of genuine extension and a poor one of the five-hour day -- which
# makes the mismeasurement non-classical and concentrated at one margin.
#
# A school counts as recorded when its expanded share reaches the treatment
# threshold, i.e. the same rule the treatment indicator uses (02), so the rates
# describe the variable the estimates are built on rather than a looser one.
#
# No estimation: reads the panel, tabulates, writes three tables.
# ==============================================================================

if (!requireNamespace("here", quietly = TRUE))
  install.packages("here", repos = "https://cloud.r-project.org")
source(here::here("scripts", "00_setup.R"))

HOURS_WAVES <- c(2021L, 2022L, 2023L)   # waves carrying the directors' hours item
THR_MAIN    <- 0.50                     # headline treatment threshold, as in 02
MIN_CASES   <- 30L                      # province cells below this are not rated

panel <- as.data.table(readRDS(file.path(DIR_PROCESSED, "panel.rds")))

d <- panel[year %in% HOURS_WAVES & !is.na(dir_hours_g6),
           .(school_id, year, province, share_expanded,
             hours = parse_hours(dir_hours_g6))]
d <- d[!is.na(hours)]
d[, recorded := as.integer(!is.na(share_expanded) & share_expanded >= THR_MAIN)]

# Coarse bins are the only ones comparable across waves: 2021 and 2022 top out at
# "6 hours or more", and only 2023 separates 6 from 7+.
d[, hours_bin := fifelse(hours <= 4, "4h or less",
                 fifelse(hours == 5, "5h", "6h or more"))]
d[, hours_detail := fifelse(hours <= 4, "4h or less", paste0(hours, "h"))]
d[hours >= 7, hours_detail := "7h or more"]

rate <- function(x, by, scale)
  x[, .(scale = scale, schools = .N, recorded = sum(recorded),
        rate_pct = round(100 * mean(recorded), 1)), by = by]

by_hours <- rbind(
  rate(d,             c("year", "hours" = "hours_bin"),    "comparable across waves"),
  rate(d[year == 2023], c("year", "hours" = "hours_detail"), "2023 detail"))
setorder(by_hours, scale, year, hours)
fwrite(by_hours, file.path(DIR_TABLES, "08b_recording_by_hours.csv"))

# Two provincial cuts, and they answer different questions. The 2021 six-hour cut
# comes BEFORE the federal programme and asks whether schools record an
# unambiguous object consistently wherever they are; the 2023 five-hour cut is
# the recording rate of the programme itself, and is the column the appendix
# table of provincial programmes reports.
by_prov <- rbind(
  d[year == 2021 & hours_bin == "6h or more"][, cut := "2021, 6h or more"],
  d[year == 2023 & hours_bin == "5h"][,        cut := "2023, 5h"]
)[, .(schools = .N, recorded = sum(recorded),
      rate_pct = round(100 * mean(recorded), 1)), by = .(cut, province)]
by_prov[, rated := schools >= MIN_CASES]
setorder(by_prov, cut, -rate_pct)
fwrite(by_prov, file.path(DIR_TABLES, "08b_recording_by_province.csv"))

# The five-hour margin as scalars, because the manuscript quotes them in prose:
# how many five-hour schools there are, how many the census records, how many
# therefore sit in the comparison group, and how much of what the census calls
# extended day in 2023 is a five-hour school.
f23  <- d[year == 2023 & hours_bin == "5h"]
ext23 <- d[year == 2023 & recorded == 1]
margin <- data.table(
  quantity = c("five-hour schools (directors, 2023)",
               "of which recorded as extended day",
               "of which left in the comparison group",
               "% of census extended day that is a five-hour school",
               "% of census extended day that is six hours or more"),
  value = c(nrow(f23), sum(f23$recorded), sum(f23$recorded == 0),
            round(100 * mean(ext23$hours_bin == "5h"), 1),
            round(100 * mean(ext23$hours_bin == "6h or more"), 1)))
fwrite(margin, file.path(DIR_TABLES, "08b_five_hour_margin.csv"))

# ---- What the labelling rule costs, and in which direction ------------------
# 08 labels each school time-invariantly by the PEAK grade-6 hours it ever
# reports. That rule is not neutral: a school reporting five hours in one wave
# and six in another is labelled genuine, never five-hour. The direction is
# therefore known by construction -- peak can only move a school UP -- and what
# is not known without counting is how many schools it moves, which is what
# decides whether the five-hour margin above is a tight bound or a loose one.
#
# Only schools with more than one observation can be affected, so the count is
# reported against that base rather than against all schools. The alternative
# rules are the LAST reported value and the MODAL one; neither is proposed as
# better, they exist to size the choice.
cls <- function(h) fifelse(is.na(h), NA_character_,
                   fifelse(h >= 6, "genuine", fifelse(h == 5, "five-hour", "simple")))
mode1 <- function(x) { u <- unique(x); u[which.max(tabulate(match(x, u)))] }

lab <- d[order(school_id, year), .(n_obs = .N,
                                   peak = max(hours),
                                   last = hours[.N],
                                   modal = mode1(hours)), by = school_id]
lab[, `:=`(c_peak = cls(peak), c_last = cls(last), c_modal = cls(modal))]
multi <- lab[n_obs > 1]

rule <- data.table(
  quantity = c("schools with directors' hours",
               "of which observed more than once",
               "five-hour under the peak rule",
               "five-hour under the last-reported rule",
               "five-hour under the modal rule",
               "multi-observation schools the peak rule moves out of five-hour",
               "multi-observation schools the peak rule moves into five-hour"),
  value = c(nrow(lab), nrow(multi),
            lab[c_peak  == "five-hour", .N],
            lab[c_last  == "five-hour", .N],
            lab[c_modal == "five-hour", .N],
            multi[c_peak != "five-hour" & (c_last == "five-hour" | c_modal == "five-hour"), .N],
            multi[c_peak == "five-hour" & c_last != "five-hour" & c_modal != "five-hour", .N]))
fwrite(rule, file.path(DIR_TABLES, "08b_labelling_rule.csv"))

# The claim the section rests on is the monotone gradient, so it is asserted
# rather than left to the reader: recording must rise with reported hours in
# every wave, and the five-hour rate must stay below the six-hour one.
chk <- dcast(by_hours[scale == "comparable across waves"],
             year ~ hours, value.var = "rate_pct")
if (nrow(chk) != length(HOURS_WAVES))
  stop("08b: expected one row per hours wave")
if (any(chk[["4h or less"]] >= chk[["5h"]]) ||
    any(chk[["5h"]] >= chk[["6h or more"]]))
  stop("08b: recording rate is not monotone in reported hours")

# ---- Reported hours by management sector, before and after the programme -----
# The two sectors sit under different rules. State schools run the schedule their
# jurisdiction sets; privately managed ones are authorised and supervised by the
# same jurisdiction but formulate their own study plans, so the provincial rule
# is a floor for them rather than a timetable (LEN 26.206, arts. 62-63).
# The prediction is that the private sector should sit ABOVE the floor where the
# floor is short -- and that the federal programme, which funds state management
# only, should move the state sector and not the private one. Both are testable
# on the directors' hours, and the ordering reverses between 2021 and 2023.
sec <- panel[year %in% HOURS_WAVES & !is.na(dir_hours_g6) & !is.na(sector),
             .(school_id, year, sector, hours = parse_hours(dir_hours_g6))]
sec <- sec[!is.na(hours), .(schools = .N,
                            pct_4h_or_less = round(100 * mean(hours <= 4), 1),
                            pct_5h         = round(100 * mean(hours == 5), 1),
                            pct_6h_or_more = round(100 * mean(hours >= 6), 1),
                            mean_hours     = round(mean(hours), 2)),
           by = .(year, sector)]
setorder(sec, year, sector)
fwrite(sec, file.path(DIR_TABLES, "08b_hours_by_sector.csv"))
