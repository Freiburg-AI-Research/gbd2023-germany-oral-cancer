# 11_crossvalidation_rebuild.R
# Rebuild the GBD-vs-ZfKD cross-validation to address reviewer critiques:
#   (a) compare on the SAME standard population: re-standardise GBD age-specific
#       rates to Segi (world) and Europe-1976 standards, which ZfKD provides;
#   (b) replace Pearson correlation (co-movement, not agreement) with agreement
#       statistics: annual GBD/ZfKD ratio, mean ratio with 95% CI, Bland-Altman
#       limits on the log scale, and a direct comparison of log-linear slopes
#       (EAPC) with a test for difference;
#   (c) report ALL four series (incidence/deaths x male/female), not a subset.
#
# Inputs : GBD 2023 cause CSV (age-specific rates with UIs), zfkd_long.csv
# Outputs: eTable crossval_agreement.csv, eFigure 3 (trajectories + ratios),
#          crossval_restandardised_series.csv

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork); library(scales)
})

DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause", "IHME-GBD_2023_DATA-684738bc-1.csv")
ZFKD_LONG <- file.path(DATA_DIR, "output", "tables", "zfkd_long.csv")
OUT_DIR   <- file.path(DATA_DIR, "output")
FIG_DIR   <- file.path(OUT_DIR, "figures")
TBL_DIR   <- file.path(OUT_DIR, "tables")

YR_LO <- 1999L; YR_STAT_HI <- 2022L; YR_HI <- 2023L   # 2023 ZfKD provisional -> plotted, excluded from stats
W_CM <- 18; PNG_DPI <- 250                             # 18 cm * 250/2.54 = 1772 px < 2000

AGE_RAW <- c("<5 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years",
             "30-34 years","35-39 years","40-44 years","45-49 years","50-54 years","55-59 years",
             "60-64 years","65-69 years","70-74 years","75-79 years","80-84 years",
             "85-89 years","90-94 years","95+ years")
AGE18 <- c(AGE_RAW[1:17], "85+ years")

# Standard populations (per 100,000), 18 groups with 85+ pooled
SEGI <- c(12000,10000,9000,9000,8000,8000,6000,6000,6000,6000,5000,4000,4000,3000,2000,1000,500,500)
EURO1976 <- c(8000,7000,7000,7000,7000,7000,7000,7000,7000,7000,7000,6000,5000,4000,3000,2000,1000,1000)
names(SEGI) <- names(EURO1976) <- AGE18
stopifnot(sum(SEGI) == 1e5, sum(EURO1976) == 1e5)

# ---- GBD age-specific rates -> pooled 85+ -> re-standardised ASRs -----------
d <- fread(CAUSE_CSV)
gall <- d[location_name == "Germany" & sex_name %in% c("Male","Female") &
          age_name %in% AGE_RAW & year >= YR_LO & year <= YR_HI]
# Population per (sex, year, age) is measure-independent: back-derive it from
# whichever measure has a non-zero rate (young age groups have zero oral-cancer
# cases, so deriving it from Incidence alone yields 0/0).
pop_n <- gall[metric_name == "Number", .(measure_name, sex = sex_name, year, age_name, n = val)]
pop_r <- gall[metric_name == "Rate",   .(measure_name, sex = sex_name, year, age_name, r = val)]
pop_m <- merge(pop_n, pop_r, by = c("measure_name","sex","year","age_name"))[r > 0, pop_i := n / r * 1e5]
pop <- pop_m[!is.na(pop_i), .(pop = median(pop_i)), by = .(sex, year, age_name)]

g  <- gall[measure_name %in% c("Incidence","Deaths")]
gn <- g[metric_name == "Number", .(measure = measure_name, sex = sex_name, year, age_name,
                                   n = val, n_lo = lower, n_hi = upper)]
gp <- merge(gn, pop, by = c("sex","year","age_name"), all.x = TRUE)
# Cells with no derivable population have zero cases in every GBD measure; their
# contribution to the ASR is zero regardless of denominator. Verify, then zero them.
no_pop <- gp[is.na(pop)]
cat(sprintf("Age-cells without derivable population: %d (all with zero cases: %s)\n",
            nrow(no_pop), all(no_pop$n == 0)))
stopifnot(all(no_pop$n == 0))
gp[, age := fifelse(age_name %in% c("85-89 years","90-94 years","95+ years"), "85+ years", age_name)]
# Pool 85+: rate = sum(n)/sum(pop); UI bounds pooled the same way (conservative)
ga <- gp[, .(n = sum(n), n_lo = sum(n_lo), n_hi = sum(n_hi), pop = sum(pop, na.rm = TRUE)),
         by = .(measure, sex, year, age)]
ga[pop > 0,  `:=`(rate = n / pop * 1e5, rate_lo = n_lo / pop * 1e5, rate_hi = n_hi / pop * 1e5)]
ga[pop == 0, `:=`(rate = 0, rate_lo = 0, rate_hi = 0)]
fwrite(pop, file.path(TBL_DIR, "gbd_germany_population_by_age.csv"))

restandardise <- function(dt, w, label) {
  dt[, .(asr    = sum(rate    * w[age]) / 1e5,
         asr_lo = sum(rate_lo * w[age]) / 1e5,
         asr_hi = sum(rate_hi * w[age]) / 1e5,
         standard_short = label),
     by = .(measure, sex, year)]
}
gbd_std <- rbind(restandardise(ga, SEGI, "Segi"), restandardise(ga, EURO1976, "Europe1976"))
gbd_std[, source := "GBD 2023"]

# Sanity: GBD's own world-standard ASR vs our Segi re-standardisation (should differ,
# since GBD uses its own 2023 world standard; printed for transparency)
gbd_own <- d[location_name == "Germany" & sex_name %in% c("Male","Female") &
             measure_name %in% c("Incidence","Deaths") & age_name == "Age-standardized" &
             metric_name == "Rate" & year == 2022, .(measure = measure_name, sex = sex_name, gbd_world_std = val)]
chk <- merge(gbd_std[year == 2022 & standard_short == "Segi", .(measure, sex, segi = asr)], gbd_own, by = c("measure","sex"))
cat("\n2022 check: GBD own world-standard ASR vs re-standardised to Segi\n"); print(chk[, .(measure, sex, gbd_world_std = round(gbd_world_std, 2), segi = round(segi, 2))])

# ---- ZfKD ------------------------------------------------------------------
z <- fread(ZFKD_LONG)
z <- z[year >= YR_LO & year <= YR_HI, .(measure, sex, year, asr = asr_total, asr_lo = NA_real_, asr_hi = NA_real_,
                                        standard_short, source = "ZfKD")]

series <- rbind(gbd_std[, .(measure, sex, year, asr, asr_lo, asr_hi, standard_short, source)], z)
series[, provisional := source == "ZfKD" & year == 2023L]
fwrite(series, file.path(TBL_DIR, "crossval_restandardised_series.csv"))

# ---- Agreement statistics (1999-2022, both standards) -----------------------
wide <- dcast(series[year <= YR_STAT_HI], measure + sex + year + standard_short ~ source, value.var = "asr")
wide <- wide[!is.na(`GBD 2023`) & !is.na(ZfKD)]
wide[, `:=`(ratio = `GBD 2023` / ZfKD, logratio = log(`GBD 2023` / ZfKD))]

eapc <- function(y, x) {   # log-linear slope -> EAPC with 95% CI (unname: c() would concatenate names)
  f <- lm(log(y) ~ x); b <- unname(coef(f)[2]); se <- unname(sqrt(vcov(f)[2, 2]))
  c(eapc = (exp(b) - 1) * 100, lo = (exp(b - 1.96 * se) - 1) * 100, hi = (exp(b + 1.96 * se) - 1) * 100)
}
slope_diff_p <- function(dt) {  # interaction test: do the log-linear slopes differ between sources?
  long <- melt(dt, id.vars = "year", measure.vars = c("GBD 2023","ZfKD"), variable.name = "src", value.name = "asr")
  f <- lm(log(asr) ~ year * src, data = long)
  b <- unname(coef(f)["year:srcZfKD"]); se <- unname(sqrt(vcov(f)["year:srcZfKD","year:srcZfKD"]))
  2 * pt(-abs(b / se), df = f$df.residual)
}

agree <- wide[, {
  n <- .N; mlr <- mean(logratio); sdlr <- sd(logratio); se <- sdlr / sqrt(n)
  eg <- eapc(`GBD 2023`, year); ez <- eapc(ZfKD, year)
  list(n_years = n,
       mean_ratio = exp(mlr), ratio_lo = exp(mlr - qt(.975, n - 1) * se), ratio_hi = exp(mlr + qt(.975, n - 1) * se),
       ba_lower = exp(mlr - 1.96 * sdlr), ba_upper = exp(mlr + 1.96 * sdlr),
       ratio_min = min(ratio), ratio_max = max(ratio),
       eapc_gbd = eg[["eapc"]], eapc_gbd_lo = eg[["lo"]], eapc_gbd_hi = eg[["hi"]],
       eapc_zfkd = ez[["eapc"]], eapc_zfkd_lo = ez[["lo"]], eapc_zfkd_hi = ez[["hi"]],
       slope_diff_p = slope_diff_p(.SD),
       pearson_r = cor(`GBD 2023`, ZfKD))
}, by = .(standard_short, measure, sex)]
setorder(agree, standard_short, measure, sex)
fwrite(agree, file.path(TBL_DIR, "crossval_agreement.csv"))

cat("\n=== Agreement GBD 2023 vs ZfKD, 1999-2022, same standard population ===\n")
pr <- copy(agree)[, `:=`(
  `mean ratio (95% CI)` = sprintf("%.2f (%.2f-%.2f)", mean_ratio, ratio_lo, ratio_hi),
  `BA limits`           = sprintf("%.2f-%.2f", ba_lower, ba_upper),
  `EAPC GBD`            = sprintf("%+.2f (%+.2f to %+.2f)", eapc_gbd, eapc_gbd_lo, eapc_gbd_hi),
  `EAPC ZfKD`           = sprintf("%+.2f (%+.2f to %+.2f)", eapc_zfkd, eapc_zfkd_lo, eapc_zfkd_hi),
  `slope diff p`        = sprintf("%.3f", slope_diff_p),
  r                     = sprintf("%.2f", pearson_r))]
print(pr[, .(standard_short, measure, sex, `mean ratio (95% CI)`, `BA limits`, `EAPC GBD`, `EAPC ZfKD`, `slope diff p`, r)])

# ---- Sensitivity: complete-registration era only (2009-2022) ----------------
# German population-based registration reached nationwide coverage around 2009;
# earlier ZfKD incidence is partly estimated and rising completeness can mimic a
# rising trend. If slopes still diverge over 2009-2022, the disagreement is not a
# completeness artefact.
wide09 <- wide[year >= 2009L & standard_short == "Segi"]
agree09 <- wide09[, {
  n <- .N; mlr <- mean(logratio); sdlr <- sd(logratio); se <- sdlr / sqrt(n)
  eg <- eapc(`GBD 2023`, year); ez <- eapc(ZfKD, year)
  list(n_years = n, mean_ratio = exp(mlr), ratio_lo = exp(mlr - qt(.975, n - 1) * se), ratio_hi = exp(mlr + qt(.975, n - 1) * se),
       eapc_gbd = eg[["eapc"]], eapc_gbd_lo = eg[["lo"]], eapc_gbd_hi = eg[["hi"]],
       eapc_zfkd = ez[["eapc"]], eapc_zfkd_lo = ez[["lo"]], eapc_zfkd_hi = ez[["hi"]],
       slope_diff_p = slope_diff_p(.SD))
}, by = .(measure, sex)]
setorder(agree09, measure, sex)
fwrite(agree09, file.path(TBL_DIR, "crossval_agreement_2009_2022.csv"))
cat("\n=== Sensitivity: 2009-2022 (nationwide registration era), Segi ===\n")
print(agree09[, .(measure, sex,
  `mean ratio (95% CI)` = sprintf("%.2f (%.2f-%.2f)", mean_ratio, ratio_lo, ratio_hi),
  `EAPC GBD`  = sprintf("%+.2f (%+.2f to %+.2f)", eapc_gbd, eapc_gbd_lo, eapc_gbd_hi),
  `EAPC ZfKD` = sprintf("%+.2f (%+.2f to %+.2f)", eapc_zfkd, eapc_zfkd_lo, eapc_zfkd_hi),
  `slope diff p` = sprintf("%.3f", slope_diff_p))])

# ---- eFigure 3: trajectories (Segi) + ratio panels --------------------------
plot_dt <- series[standard_short == "Segi"]
plot_dt[, measure := factor(measure, levels = c("Incidence","Deaths"))]
plot_dt[, sex := factor(sex, levels = c("Male","Female"))]
cols <- c("GBD 2023" = "#1f77b4", "ZfKD" = "#d62728")

p_traj <- ggplot(plot_dt, aes(year, asr, colour = source)) +
  geom_ribbon(data = plot_dt[source == "GBD 2023"], aes(ymin = asr_lo, ymax = asr_hi, fill = source),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(shape = provisional), size = 1.3) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  scale_colour_manual(values = cols) + scale_fill_manual(values = cols) +
  facet_grid(measure ~ sex, scales = "free_y") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .05))) +
  scale_x_continuous(breaks = seq(2000, 2020, 10)) +
  labs(x = NULL, y = "ASR per 100,000 (Segi world standard)", colour = NULL, fill = NULL,
       title = "A  Age-standardised rates on a common standard population") +
  guides(fill = "none") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(), plot.title.position = "plot")

ratio_dt <- dcast(plot_dt[source %in% c("GBD 2023","ZfKD")], measure + sex + year + provisional ~ source, value.var = "asr")
ratio_dt <- ratio_dt[!is.na(`GBD 2023`) & !is.na(ZfKD)][, ratio := `GBD 2023` / ZfKD]
band <- agree[standard_short == "Segi", .(measure = factor(measure, levels = c("Incidence","Deaths")),
                                           sex = factor(sex, levels = c("Male","Female")), mean_ratio, ba_lower, ba_upper)]
p_ratio <- ggplot(ratio_dt, aes(year, ratio)) +
  geom_rect(data = band, aes(xmin = -Inf, xmax = Inf, ymin = ba_lower, ymax = ba_upper), inherit.aes = FALSE,
            fill = "grey85", alpha = 0.6) +
  geom_hline(data = band, aes(yintercept = mean_ratio), linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = 1, colour = "black", linewidth = 0.3) +
  geom_line(colour = "#6a3d9a", linewidth = 0.6) +
  geom_point(aes(shape = provisional), colour = "#6a3d9a", size = 1.3) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  facet_grid(measure ~ sex) +
  scale_y_continuous(limits = c(0.5, 2.0), breaks = c(0.5, 1, 1.5, 2)) +
  scale_x_continuous(breaks = seq(2000, 2020, 10)) +
  labs(x = "Year", y = "GBD / ZfKD ratio",
       title = "B  Ratio of GBD 2023 to ZfKD (dashed: mean ratio; grey: Bland-Altman 95% limits)") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(), plot.title.position = "plot")

fig <- p_traj / p_ratio + plot_layout(heights = c(1.15, 1))
ggsave(file.path(FIG_DIR, "figureS3_gbd_vs_zfkd.pdf"), fig, width = W_CM, height = 17, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figureS3_gbd_vs_zfkd.png"), fig, width = W_CM, height = 17, units = "cm", dpi = PNG_DPI)
cat("\neFigure 3 (rebuilt) and crossval_agreement.csv written.\n")
