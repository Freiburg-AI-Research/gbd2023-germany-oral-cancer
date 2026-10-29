# 15_age_specific_trends.R
# Age-specific divergence of lip and oral cavity cancer trends in Germany:
#   - rates by age band (20-49, 50-69, 70+) and 10-year bands, by sex, 1990-2023
#   - EAPC per band (log-linear), incidence and deaths
#   - mean age at diagnosis (incidence-weighted) and at death, by sex and year
# The aggregate trends mask an age shift: incidence collapses under 50 and
# rises at 70+, consistent with the birth-cohort pattern from Nordpred.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork); library(scales) })
DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause", "IHME-GBD_2023_DATA-684738bc-1.csv")
OUT_DIR   <- file.path(DATA_DIR, "output"); FIG_DIR <- file.path(OUT_DIR, "figures"); TBL_DIR <- file.path(OUT_DIR, "tables")
W_CM <- 18; PNG_DPI <- 250

AGE_RAW <- c("<5 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years",
             "40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years",
             "75-79 years","80-84 years","85-89 years","90-94 years","95+ years")
AGE_MID <- c(seq(2.5, 92.5, by = 5), 97.5); names(AGE_MID) <- AGE_RAW
BAND3  <- c(rep("<20", 4), rep("20-49", 6), rep("50-69", 4), rep("70+", 6)); names(BAND3) <- AGE_RAW
BAND10 <- c(rep("<30", 6), rep("30-39", 2), rep("40-49", 2), rep("50-59", 2), rep("60-69", 2), rep("70-79", 2), rep("80+", 4)); names(BAND10) <- AGE_RAW

d <- fread(CAUSE_CSV)
g <- d[location_name == "Germany" & age_name %in% AGE_RAW & sex_name %in% c("Male","Female") &
       measure_name %in% c("Incidence","Deaths")]
n <- g[metric_name == "Number", .(measure = measure_name, sex = sex_name, year, age_name, n = val)]
# population per (sex, year, age) from any measure with a non-zero rate
pn <- d[location_name == "Germany" & age_name %in% AGE_RAW & sex_name %in% c("Male","Female") & metric_name == "Number", .(measure_name, sex = sex_name, year, age_name, n = val)]
pr <- d[location_name == "Germany" & age_name %in% AGE_RAW & sex_name %in% c("Male","Female") & metric_name == "Rate",   .(measure_name, sex = sex_name, year, age_name, r = val)]
pop <- merge(pn, pr, by = c("measure_name","sex","year","age_name"))[r > 0][, .(pop = median(n / r * 1e5)), by = .(sex, year, age_name)]
x <- merge(n, pop, by = c("sex","year","age_name"), all.x = TRUE)
x[, `:=`(band3 = BAND3[age_name], band10 = BAND10[age_name], age_mid = AGE_MID[age_name])]

eapc <- function(y, t) { ok <- y > 0; f <- lm(log(y[ok]) ~ t[ok]); b <- unname(coef(f)[2]); s <- unname(sqrt(vcov(f)[2, 2]))
  list(eapc = (exp(b) - 1) * 100, lo = (exp(b - 1.96 * s) - 1) * 100, hi = (exp(b + 1.96 * s) - 1) * 100) }

band_rates <- function(bandcol) {
  br <- x[!is.na(pop), .(rate = sum(n) / sum(pop) * 1e5, n = sum(n)), by = c("measure","sex","year", bandcol)]
  setnames(br, bandcol, "band"); br[band != "<20" & band != "<30"] }
b3 <- band_rates("band3"); b10 <- band_rates("band10")
t3  <- b3[,  c(eapc(rate, year), rate_1990 = rate[year == 1990], rate_2023 = rate[year == 2023]), by = .(measure, sex, band)]
t10 <- b10[, c(eapc(rate, year), rate_1990 = rate[year == 1990], rate_2023 = rate[year == 2023]), by = .(measure, sex, band)]
t3[,  rel_change_pct := (rate_2023 / rate_1990 - 1) * 100]; t10[, rel_change_pct := (rate_2023 / rate_1990 - 1) * 100]
fwrite(t3,  file.path(TBL_DIR, "age_band3_eapc.csv")); fwrite(t10, file.path(TBL_DIR, "age_band10_eapc.csv"))
fwrite(b3,  file.path(TBL_DIR, "age_band3_rates_1990_2023.csv"))
cat("=== EAPC 1990-2023 by 10-year age band (crude within-band rate) ===\n")
print(t10[order(measure, sex, band), .(measure, sex, band, rate_1990 = round(rate_1990, 2), rate_2023 = round(rate_2023, 2),
        change = sprintf("%+.0f%%", rel_change_pct), EAPC = sprintf("%+.2f (%+.2f to %+.2f)", eapc, lo, hi))])

# Mean age at diagnosis / death (count-weighted age-group midpoint)
ma <- x[, .(mean_age = sum(n * age_mid) / sum(n)), by = .(measure, sex, year)]
fwrite(ma, file.path(TBL_DIR, "mean_age_diagnosis_death.csv"))
cat("\n=== Mean age at diagnosis (Incidence) and at death (Deaths), years ===\n")
print(dcast(ma[year %in% c(1990, 2000, 2010, 2023)], measure + sex ~ year, value.var = "mean_age")[, lapply(.SD, function(v) if (is.numeric(v)) round(v, 1) else v)])

# Figure: (A) incidence rate by age band over time (log scale), by sex; (B) mean age at diagnosis and death
b10[, band := factor(band, levels = c("30-39","40-49","50-59","60-69","70-79","80+"))]
b10[, sex := factor(sex, levels = c("Male","Female"))]
pal <- c("30-39" = "#fee08b", "40-49" = "#fdae61", "50-59" = "#f46d43", "60-69" = "#d73027", "70-79" = "#a50026", "80+" = "#4d004b")
pA <- ggplot(b10[measure == "Incidence"], aes(year, rate, colour = band)) +
  geom_line(linewidth = 0.7) + facet_wrap(~ sex) +
  scale_colour_manual(values = pal, name = "Age (years)") +
  scale_y_log10(breaks = c(0.5, 1, 2, 5, 10, 20, 50)) +
  scale_x_continuous(breaks = seq(1990, 2020, 10)) +
  labs(x = NULL, y = "Incidence per 100,000 (log scale)", title = "A  Age-specific incidence, 1990–2023") +
  theme_minimal(base_size = 9) + theme(strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(), plot.title.position = "plot", legend.position = "right")
ma[, measure := factor(fifelse(measure == "Incidence", "Age at diagnosis", "Age at death"), levels = c("Age at diagnosis","Age at death"))]
ma[, sex := factor(sex, levels = c("Male","Female"))]
pB <- ggplot(ma, aes(year, mean_age, colour = sex, linetype = measure)) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(Male = "#1f77b4", Female = "#d62728"), name = NULL) +
  scale_linetype_manual(values = c(`Age at diagnosis` = "solid", `Age at death` = "dashed"), name = NULL) +
  scale_x_continuous(breaks = seq(1990, 2020, 10)) +
  labs(x = "Year", y = "Mean age (years)", title = "B  Mean age at diagnosis and at death") +
  theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot", legend.position = "right")
fig <- pA / pB + plot_layout(heights = c(1.3, 1))
ggsave(file.path(FIG_DIR, "figure_age_specific_trends.pdf"), fig, width = W_CM, height = 14, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figure_age_specific_trends.png"), fig, width = W_CM, height = 14, units = "cm", dpi = PNG_DPI)
ggsave(file.path(FIG_DIR, "figure_age_specific_trends.tif"), fig, width = W_CM, height = 14, units = "cm", dpi = 600, device = "tiff", compression = "lzw", type = "cairo")
cat("\nfigure_age_specific_trends written.\n")
