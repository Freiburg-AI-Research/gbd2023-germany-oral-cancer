# 07_population_sensitivity.R
# Sensitivity analysis: replace static-population assumption in Nordpred
# projection with UN World Population Prospects 2024 medium-variant projection
# for Germany (5-year age x sex x year, all comparable to GBD/Nordpred bins).
#
# Output: eFigure 4 (or eTable) comparing static vs WPP absolute count projections.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr)
  library(scales); library(nordpred); library(wpp2024)
})

DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause",
                      "IHME-GBD_2023_DATA-684738bc-1.csv")
OUT_DIR  <- file.path(DATA_DIR, "output")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TBL_DIR  <- file.path(OUT_DIR, "tables")
DATA_OUT <- file.path(OUT_DIR, "data")

W_CM <- 18; PNG_DPI <- 250

AGE_5YR_RAW <- c("<5 years","5-9 years","10-14 years","15-19 years","20-24 years",
                 "25-29 years","30-34 years","35-39 years","40-44 years","45-49 years",
                 "50-54 years","55-59 years","60-64 years","65-69 years","70-74 years",
                 "75-79 years","80-84 years","85-89 years","90-94 years","95+ years")
AGE_5YR <- c(AGE_5YR_RAW[1:17], "85+ years")

HIST_PERIODS <- list(c(1994,1998),c(1999,2003),c(2004,2008),c(2009,2013),c(2014,2018),c(2019,2023))
HIST_PERIOD_LBL <- c("1994-1998","1999-2003","2004-2008","2009-2013","2014-2018","2019-2023")
HIST_PERIOD_MID <- c(1996, 2001, 2006, 2011, 2016, 2021)

FCST_PERIOD_LBL <- c("2024-2028","2029-2033","2034-2038","2039-2043")
FCST_PERIOD_MID <- c(2026, 2031, 2036, 2041)
NPROJ <- length(FCST_PERIOD_LBL)
START_EST_AGE <- 6; START_USE_AGE <- 6

# ---- WPP 2024 Germany medium-variant population projection -----------------
data(popMprojMed); data(popFprojMed)
# Each column is 1 July of the given year, in THOUSANDS

# Map our 18 Nordpred age bins to WPP age groups
wpp_age_map <- c("<5 years" = "0-4", "5-9 years" = "5-9", "10-14 years" = "10-14",
                 "15-19 years" = "15-19", "20-24 years" = "20-24", "25-29 years" = "25-29",
                 "30-34 years" = "30-34", "35-39 years" = "35-39", "40-44 years" = "40-44",
                 "45-49 years" = "45-49", "50-54 years" = "50-54", "55-59 years" = "55-59",
                 "60-64 years" = "60-64", "65-69 years" = "65-69", "70-74 years" = "70-74",
                 "75-79 years" = "75-79", "80-84 years" = "80-84")
# For 85+, sum WPP 85-89, 90-94, 95-99, 100+
wpp_85plus <- c("85-89","90-94","95-99","100+")

# Pull Germany rows
ger_M <- as.data.table(popMprojMed[popMprojMed$name == "Germany",])
ger_F <- as.data.table(popFprojMed[popFprojMed$name == "Germany",])

# Period midpoint years we need: 2026, 2031, 2036, 2041
# WPP has 2025, 2030, 2035, 2040 — close enough. Linear interpolation between
# adjacent 5-year columns gives the midpoint.
interp_pop <- function(ger_dt, age_wpp, year_target) {
  row <- ger_dt[age == age_wpp]
  if (nrow(row) == 0) return(NA_real_)
  # Get adjacent 5-year columns
  yrs <- as.integer(grep("^[0-9]+$", names(ger_dt), value = TRUE))
  lower_y <- max(yrs[yrs <= year_target])
  upper_y <- min(yrs[yrs >= year_target])
  if (lower_y == upper_y) return(as.numeric(row[[as.character(lower_y)]]))
  v_lo <- as.numeric(row[[as.character(lower_y)]])
  v_hi <- as.numeric(row[[as.character(upper_y)]])
  v_lo + (v_hi - v_lo) * (year_target - lower_y) / (upper_y - lower_y)
}

build_wpp_pyr <- function(sex) {
  ger_dt <- if (sex == "Male") ger_M else ger_F
  # For each of our 18 age bins and 4 future period-midpoint years
  res <- matrix(NA_real_, nrow = 18, ncol = NPROJ,
                dimnames = list(AGE_5YR, FCST_PERIOD_LBL))
  for (a in AGE_5YR) {
    for (i in seq_along(FCST_PERIOD_MID)) {
      yr <- FCST_PERIOD_MID[i]
      if (a == "85+ years") {
        # Sum WPP 85-89, 90-94, 95-99, 100+
        vals <- sapply(wpp_85plus, function(w) interp_pop(ger_dt, w, yr))
        res[a, i] <- sum(vals, na.rm = TRUE)
      } else {
        res[a, i] <- interp_pop(ger_dt, wpp_age_map[[a]], yr)
      }
    }
  }
  # Convert thousands -> persons, then to 5-year person-years (multiply by 5)
  res * 1000 * 5
}

wpp_pyr_M <- build_wpp_pyr("Male")
wpp_pyr_F <- build_wpp_pyr("Female")
wpp_pyr_B <- wpp_pyr_M + wpp_pyr_F

cat("WPP 2024 Germany medium-variant population projection by 5-year period:\n")
cat("\nBoth sexes (thousands of person-years per 5-yr period):\n")
print(round(wpp_pyr_B / 1000))

# ---- Re-run Nordpred for Both sexes with WPP populations -------------------
d <- fread(CAUSE_CSV)
dg <- d[location_name == "Germany" & age_name %in% AGE_5YR_RAW]
dg_n <- dg[metric_name == "Number",
           .(measure_name, sex_name, year, age_name, n = val)]
dg_r <- dg[metric_name == "Rate",
           .(measure_name, sex_name, year, age_name, r = val)]
dg_p <- merge(dg_n, dg_r, by = c("measure_name","sex_name","year","age_name"))
dg_p[, pop := ifelse(r > 0, n / r * 1e5, NA_real_)]
dg_p[, age_bin := fifelse(age_name %in% c("85-89 years","90-94 years","95+ years"),
                          "85+ years", age_name)]
dg <- dg_p[, .(n = sum(n, na.rm = TRUE), pop = sum(pop, na.rm = TRUE)),
           by = .(measure_name, sex_name, year, age_bin)]
setnames(dg, "age_bin", "age_name")
dg[, age_name := factor(age_name, levels = AGE_5YR)]

build_mats_with_wpp <- function(sex, measure, wpp_future_pyr) {
  df <- dg[sex_name == sex & measure_name == measure,
           .(year, age_name, cases = n, pyr_hist = pop)]
  df[, period := NA_character_]
  for (i in seq_along(HIST_PERIODS)) {
    rng <- HIST_PERIODS[[i]]
    df[year >= rng[1] & year <= rng[2], period := HIST_PERIOD_LBL[i]]
  }
  agg <- df[!is.na(period), .(cases = sum(cases, na.rm = TRUE),
                              pyr = sum(pyr_hist, na.rm = TRUE)),
            by = .(age_name, period)]
  cm <- dcast(agg, age_name ~ period, value.var = "cases", drop = FALSE)
  cm <- as.matrix(cm[, HIST_PERIOD_LBL, with = FALSE]); rownames(cm) <- AGE_5YR
  pm <- dcast(agg, age_name ~ period, value.var = "pyr", drop = FALSE)
  pm <- as.matrix(pm[, HIST_PERIOD_LBL, with = FALSE]); rownames(pm) <- AGE_5YR

  cm[is.na(cm)] <- 0; cm <- round(cm); storage.mode(cm) <- "integer"
  pm[is.na(pm) | pm <= 0] <- 1
  # Stack historical + WPP future
  pm_full <- cbind(pm, wpp_future_pyr)
  list(cases = as.data.frame(cm), pyr = as.data.frame(pm_full))
}

run_one <- function(sex, measure, wpp_future_pyr) {
  mats <- build_mats_with_wpp(sex, measure, wpp_future_pyr)
  fit <- nordpred(cases = mats$cases, pyr = mats$pyr,
                  startestage = START_EST_AGE, startuseage = START_USE_AGE,
                  noperiods = NPROJ,
                  cuttrend = c(0, 0.25, 0.5, 0.75),
                  linkfunc = "power5")
  counts <- nordpred.getpred(fit, incidence = FALSE, standpop = NULL)
  # counts is age x periods (historical + future)
  total_counts <- colSums(counts, na.rm = TRUE)
  data.table(period = c(HIST_PERIOD_LBL, FCST_PERIOD_LBL),
             period_mid = c(HIST_PERIOD_MID, FCST_PERIOD_MID),
             total_count = as.numeric(total_counts),
             sex = sex, measure = measure)
}

# Run for Both sexes (use combined WPP pyr) for Deaths and DALYs
res_wpp <- rbind(
  run_one("Both","Deaths", wpp_pyr_B),
  run_one("Both","DALYs (Disability-Adjusted Life Years)", wpp_pyr_B),
  run_one("Both","Incidence", wpp_pyr_B)
)
res_wpp[, scenario := "WPP 2024 medium variant"]

# Also load static-population baseline from earlier run
res_static <- fread(file.path(DATA_OUT, "nordpred_projections.csv"))
res_static <- res_static[sex == "Both" & measure %in% c("Deaths","DALYs (Disability-Adjusted Life Years)","Incidence"),
                         .(period, period_mid, total_count, sex, measure)]
res_static[, scenario := "Static population (baseline)"]

cmp <- rbind(res_wpp, res_static)
fwrite(cmp, file.path(TBL_DIR, "population_sensitivity.csv"))

# Compare in projection horizon
cat("\nProjection horizon comparison (Both sexes):\n")
proj_only <- cmp[period_mid >= 2024]
print(dcast(proj_only, measure + period ~ scenario, value.var = "total_count")[, lapply(.SD, function(x) ifelse(is.numeric(x), round(x), x))])

# ---- Plot static vs WPP for absolute counts (Both sexes, Deaths & DALYs) ----
cmp[, measure_lbl := dplyr::recode(measure,
    "Deaths"="Deaths","DALYs (Disability-Adjusted Life Years)"="DALYs",
    "Incidence"="Incidence")]
cmp[, measure_lbl := factor(measure_lbl, levels = c("Incidence","Deaths","DALYs"))]
cmp[, annual_count := total_count / 5]  # convert 5-yr total -> annual avg

p4 <- ggplot(cmp, aes(x = period_mid, y = annual_count, colour = scenario,
                      linetype = scenario)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.3) +
  facet_wrap(~ measure_lbl, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = seq(2000, 2040, 10)) +
  scale_y_continuous(labels = comma_format(accuracy = 1)) +
  scale_colour_manual(values = c("Static population (baseline)" = "grey50",
                                 "WPP 2024 medium variant" = "#1f77b4")) +
  scale_linetype_manual(values = c("Static population (baseline)" = "dashed",
                                   "WPP 2024 medium variant" = "solid")) +
  labs(x = "Year (period midpoint)", y = "Annual count",
       colour = NULL, linetype = NULL,
       title = "Population-projection sensitivity for absolute-count forecasts",
       subtitle = "Both sexes, Germany; static-population baseline vs UN WPP 2024 medium-variant population projection") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom",
        plot.title.position = "plot",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "figureS4_population_sensitivity.pdf"),
       p4, width = W_CM, height = 9, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figureS4_population_sensitivity.png"),
       p4, width = W_CM, height = 9, units = "cm", dpi = PNG_DPI)
cat("\neFigure 4 (population sensitivity) written.\n")
