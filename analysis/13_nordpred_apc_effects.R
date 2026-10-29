# 13_nordpred_apc_effects.R
# Extract and report what the Nordpred age-period-cohort model actually
# estimates, rather than using it only as a black-box projector:
#   - goodness-of-fit p-value and the "recent trend" test (Nordpred's rule for
#     basing the projection on the last two periods)
#   - the linear period drift
#   - birth-cohort effects as relative rates vs a reference cohort
# Fits: sex (Male, Female) x measure (Incidence, Deaths, DALYs), same
# specification as 01_nordpred_germany.R.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(nordpred) })

DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause", "IHME-GBD_2023_DATA-684738bc-1.csv")
OUT_DIR   <- file.path(DATA_DIR, "output")
FIG_DIR   <- file.path(OUT_DIR, "figures"); TBL_DIR <- file.path(OUT_DIR, "tables")
W_CM <- 18; PNG_DPI <- 250

AGE_RAW <- c("<5 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years",
             "30-34 years","35-39 years","40-44 years","45-49 years","50-54 years","55-59 years",
             "60-64 years","65-69 years","70-74 years","75-79 years","80-84 years",
             "85-89 years","90-94 years","95+ years")
AGE18 <- c(AGE_RAW[1:17], "85+ years")
AGE_MID <- c(seq(2.5, 82.5, by = 5), 90)            # midpoints of the 18 groups
HIST <- list(c(1994,1998),c(1999,2003),c(2004,2008),c(2009,2013),c(2014,2018),c(2019,2023))
HIST_LBL <- c("1994-1998","1999-2003","2004-2008","2009-2013","2014-2018","2019-2023")
PERIOD_MID <- c(1996, 2001, 2006, 2011, 2016, 2021)
FCST_LBL <- c("2024-2028","2029-2033","2034-2038","2039-2043")
NPROJ <- 4; START <- 6

d <- fread(CAUSE_CSV)
dg <- d[location_name == "Germany" & age_name %in% AGE_RAW]
n <- dg[metric_name == "Number", .(measure = measure_name, sex = sex_name, year, age_name, n = val)]
r <- dg[metric_name == "Rate",   .(measure = measure_name, sex = sex_name, year, age_name, r = val)]
p <- merge(n, r, by = c("measure","sex","year","age_name")); p[, pop := ifelse(r > 0, n / r * 1e5, NA_real_)]
p[, age := fifelse(age_name %in% c("85-89 years","90-94 years","95+ years"), "85+ years", age_name)]
a <- p[, .(n = sum(n), pop = sum(pop, na.rm = TRUE)), by = .(measure, sex, year, age)]
a[, period := NA_character_]
for (i in seq_along(HIST)) a[year >= HIST[[i]][1] & year <= HIST[[i]][2], period := HIST_LBL[i]]

fit_one <- function(sex_i, meas_i) {
  agg <- a[!is.na(period) & sex == sex_i & measure == meas_i, .(cases = sum(n), pyr = sum(pop)), by = .(age, period)]
  cm <- dcast(agg, age ~ period, value.var = "cases"); cm <- as.data.frame(cm[match(AGE18, age), HIST_LBL, with = FALSE]); rownames(cm) <- AGE18
  pm <- dcast(agg, age ~ period, value.var = "pyr");   pm <- as.data.frame(pm[match(AGE18, age), HIST_LBL, with = FALSE]); rownames(pm) <- AGE18
  cm[is.na(cm)] <- 0; cm <- round(cm); pm[is.na(pm) | pm <= 0] <- 1
  fut <- as.data.frame(matrix(rep(pm[, 6], NPROJ), nrow = 18)); colnames(fut) <- FCST_LBL; pm <- cbind(pm, fut)
  nordpred(cases = cm, pyr = pm, startestage = START, startuseage = START, noperiods = NPROJ,
           cuttrend = c(0, 0.25, 0.5, 0.75), linkfunc = "power5")
}

combos <- CJ(sex = c("Male","Female"), measure = c("Incidence","Deaths","DALYs (Disability-Adjusted Life Years)"))
diag_rows <- list(); cohort_rows <- list()
for (i in seq_len(nrow(combos))) {
  sx <- combos$sex[i]; ms <- combos$measure[i]
  fit <- fit_one(sx, ms); cf <- coef(fit$glm); se <- sqrt(diag(vcov(fit$glm)))
  mf <- fit$glm$model                                     # model frame with Age, Period, Cohort indices
  # map cohort index -> mean birth year (period midpoint - age midpoint)
  map <- data.table(Age = as.integer(as.character(mf[["as.factor(Age)"]])),
                    Period = as.integer(mf[["Period"]]),
                    Cohort = as.integer(as.character(mf[["as.factor(Cohort)"]])))
  map[, birth := PERIOD_MID[Period] - AGE_MID[Age]]
  cmap <- map[, .(birth_mid = mean(birth), n_cells = .N), by = Cohort][order(Cohort)]
  # cohort effects (power-5 scale); cohorts without a coefficient are the reference (0)
  cmap[, `:=`(delta = 0, delta_se = 0)]
  for (k in cmap$Cohort) { nm <- paste0("as.factor(Cohort)", k)
    if (nm %in% names(cf)) cmap[Cohort == k, `:=`(delta = unname(cf[nm]), delta_se = unname(se[nm]))] }
  # linear predictor at reference age 60-64 (index 13) and period 4 (2009-2013)
  eta_base <- unname(cf["as.factor(Age)13"]) + unname(cf["Period"]) * 4 + unname(cf["as.factor(Period)4"])
  ref <- cmap[which.min(abs(birth_mid - 1950)), Cohort]             # reference: cohort born ~1950
  d_ref <- cmap[Cohort == ref, delta]
  eta0 <- eta_base + d_ref
  cmap[, rr    := ((eta0 - d_ref + delta) / eta0)^5]
  cmap[, rr_lo := ((eta0 - d_ref + delta - 1.96 * delta_se) / eta0)^5]
  cmap[, rr_hi := ((eta0 - d_ref + delta + 1.96 * delta_se) / eta0)^5]
  cmap[, `:=`(sex = sx, measure = ms, ref_cohort = ref)]
  cohort_rows[[i]] <- cmap
  # drift: per-period change on eta scale -> approximate annual % change in rate at reference
  drift <- unname(cf["Period"]); drift_pct_yr <- (((eta0 + drift / 5) / eta0)^5 - 1) * 100
  diag_rows[[i]] <- data.table(sex = sx, measure = ms,
    gof_p = fit$gofpvalue, recent_trend_used = fit$recent, recent_test_p = fit$pvaluerecent,
    drift_eta_per_period = drift, drift_pct_per_year_at_ref = drift_pct_yr,
    n_cells_fitted = nrow(mf), residual_df = fit$glm$df.residual, deviance = fit$glm$deviance)
}
diag_tbl <- rbindlist(diag_rows); cohort_tbl <- rbindlist(cohort_rows)
diag_tbl[, measure := fifelse(measure == "DALYs (Disability-Adjusted Life Years)", "DALYs", measure)]
cohort_tbl[, measure := fifelse(measure == "DALYs (Disability-Adjusted Life Years)", "DALYs", measure)]
fwrite(diag_tbl, file.path(TBL_DIR, "nordpred_diagnostics.csv"))
fwrite(cohort_tbl, file.path(TBL_DIR, "nordpred_cohort_effects.csv"))

cat("=== Nordpred diagnostics (GoF, recent-trend rule, drift) ===\n")
print(diag_tbl[, .(sex, measure, gof_p = signif(gof_p, 2), recent_trend_used, recent_test_p = signif(recent_test_p, 2),
                   drift_pct_per_year_at_ref = round(drift_pct_per_year_at_ref, 2), deviance = round(deviance, 1), residual_df)])
cat("\n=== Cohort relative rates (vs cohort born ~1950), selected cohorts ===\n")
print(cohort_tbl[order(sex, measure, birth_mid),
                 .(sex, measure, cohort_born = round(birth_mid), rr = round(rr, 2), rr_lo = round(rr_lo, 2), rr_hi = round(rr_hi, 2))])

# ---- eFigure: cohort relative rates -----------------------------------------
cohort_tbl[, measure := factor(measure, levels = c("Incidence","Deaths","DALYs"))]
cohort_tbl[, sex := factor(sex, levels = c("Male","Female"))]
pc <- ggplot(cohort_tbl, aes(birth_mid, rr, colour = sex, fill = sex)) +
  geom_ribbon(aes(ymin = rr_lo, ymax = rr_hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ measure, ncol = 3) +
  scale_colour_manual(values = c(Male = "#1f77b4", Female = "#d62728")) +
  scale_fill_manual(values = c(Male = "#1f77b4", Female = "#d62728")) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .06))) +
  labs(x = "Birth cohort (central year)", y = "Relative rate vs cohort born ~1950", colour = NULL, fill = NULL,
       title = "Birth-cohort effects from the Nordpred age-period-cohort model, Germany",
       subtitle = "Relative rate at age 60-64 in 2009-2013; approximate 95% CIs ignore parameter covariance") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(), plot.title.position = "plot")
ggsave(file.path(FIG_DIR, "figureS6_nordpred_cohort_effects.pdf"), pc, width = W_CM, height = 8, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figureS6_nordpred_cohort_effects.png"), pc, width = W_CM, height = 8, units = "cm", dpi = PNG_DPI)
ggsave(file.path(FIG_DIR, "figureS6_nordpred_cohort_effects.tif"), pc, width = W_CM, height = 8, units = "cm", dpi = 600, device = "tiff", compression = "lzw", type = "cairo")
cat("\neFigure 6 (cohort effects) written.\n")
