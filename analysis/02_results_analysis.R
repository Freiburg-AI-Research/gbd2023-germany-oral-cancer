# 02_results_analysis.R
# Comprehensive analyses for the Results section:
#   1. EAPC (estimated annual percent change) for 1990-2023 trends
#   2. Figure 2: Germany 2023 age pattern by sex
#   3. Figure 3: Risk-factor decomposition (tobacco + alcohol over time)
#   4. Table 1: Master burden table (1994 / 2023 / 2040 + EAPC), by sex
#   5. Comparison table: Germany vs Western Europe vs Global, 2023
#   6. Supplementary Figure: Sex-stratified Figure 1

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr)
  library(tidyr); library(scales); library(patchwork)
})

DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause",
                      "IHME-GBD_2023_DATA-684738bc-1.csv")
RISK_CSV  <- file.path(DATA_DIR, "data/gbd_risk",
                      "IHME-GBD_2023_DATA-6a8e77e9-1.csv")
OUT_DIR  <- file.path(DATA_DIR, "output")
FIG_DIR  <- file.path(OUT_DIR, "figures")
TBL_DIR  <- file.path(OUT_DIR, "tables")
DATA_OUT <- file.path(OUT_DIR, "data")

W_CM <- 18; PNG_DPI <- 250

# 5-year age bins (canonical order)
AGE_5YR_RAW <- c("<5 years","5-9 years","10-14 years","15-19 years","20-24 years",
                 "25-29 years","30-34 years","35-39 years","40-44 years","45-49 years",
                 "50-54 years","55-59 years","60-64 years","65-69 years","70-74 years",
                 "75-79 years","80-84 years","85-89 years","90-94 years","95+ years")

# ---- Load data ------------------------------------------------------------
cat("Loading data...\n")
d  <- fread(CAUSE_CSV)
dr <- fread(RISK_CSV)
proj <- fread(file.path(DATA_OUT, "nordpred_projections.csv"))

# ============================================================================
# 1. EAPC analysis (annual ASR; log-linear regression)
# ============================================================================
cat("\n[1/6] EAPC analysis (log-linear regression on annual ASR)...\n")

eapc_input <- d[location_name == "Germany" &
                age_name == "Age-standardized" & metric_name == "Rate" &
                measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)",
                                    "Incidence","Prevalence","YLLs (Years of Life Lost)",
                                    "YLDs (Years Lived with Disability)") &
                sex_name %in% c("Both","Male","Female"),
                .(measure = measure_name, sex = sex_name, year, asr = val)]

compute_eapc <- function(df) {
  df <- df[asr > 0]   # log requires positive
  if (nrow(df) < 5) return(data.table(eapc = NA, eapc_lo = NA, eapc_hi = NA))
  fit <- lm(log(asr) ~ year, data = df)
  b   <- coef(fit)["year"]
  se  <- summary(fit)$coefficients["year","Std. Error"]
  data.table(eapc     = (exp(b)        - 1) * 100,
             eapc_lo  = (exp(b - 1.96*se) - 1) * 100,
             eapc_hi  = (exp(b + 1.96*se) - 1) * 100)
}

eapc_all <- eapc_input[, compute_eapc(.SD), by = .(measure, sex)]
eapc_all[, eapc_fmt := sprintf("%+.2f (%+.2f to %+.2f)", eapc, eapc_lo, eapc_hi)]
fwrite(eapc_all, file.path(TBL_DIR, "table_eapc_1990_2023.csv"))
cat("  EAPC table written.\n")
print(eapc_all[, .(measure, sex, eapc_fmt)])

# ============================================================================
# 2. Figure 2: Age pattern in Germany 2023, by sex
# ============================================================================
cat("\n[2/6] Figure 2: age pattern...\n")

age_2023 <- d[location_name == "Germany" & year == 2023 &
              age_name %in% AGE_5YR_RAW & metric_name == "Rate" &
              sex_name %in% c("Male","Female") &
              measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)",
                                  "Incidence")]
age_2023[, age_name := factor(age_name, levels = AGE_5YR_RAW)]
age_2023[, age_mid  := as.numeric(age_name) * 5 - 2.5]   # 0-4 -> 2.5, 5-9 -> 7.5, etc.

age_2023[, measure_lbl := dplyr::recode(measure_name,
            "Deaths" = "Deaths",
            "DALYs (Disability-Adjusted Life Years)" = "DALYs",
            "Incidence" = "Incidence")]
age_2023[, measure_lbl := factor(measure_lbl, levels = c("Incidence","Deaths","DALYs"))]

p2 <- ggplot(age_2023, aes(x = age_mid, y = val, colour = sex_name, group = sex_name)) +
  geom_line(linewidth = 0.7) + geom_point(size = 1.3) +
  facet_wrap(~ measure_lbl, scales = "free_y", ncol = 3) +
  # X-axis: place breaks at actual age-group midpoints so the highest label is
  # "95+" rather than the artefactual "100-104"
  scale_x_continuous(breaks = c(2.5, 22.5, 42.5, 62.5, 82.5, 97.5),
                     labels = c("0–4","20–24","40–44","60–64","80–84","95+"),
                     limits = c(0, 100),
                     expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(labels = comma_format(accuracy = 1), limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.06))) +
  scale_colour_manual(values = c(Male = "#1f77b4", Female = "#d62728")) +
  labs(x = "Age group (years)", y = "Rate per 100,000",
       colour = NULL,
       title = "Age pattern of lip and oral cavity cancer in Germany, 2023",
       subtitle = "Sex-specific rates per 100,000 population by 5-year age group; y-axes differ between panels") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        plot.title.position = "plot",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0))

ggsave(file.path(FIG_DIR, "figure2_germany_age_pattern_2023.pdf"),
       p2, width = W_CM, height = 9, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figure2_germany_age_pattern_2023.png"),
       p2, width = W_CM, height = 9, units = "cm", dpi = PNG_DPI)
cat("  Figure 2 written.\n")

# ============================================================================
# 3. Figure 3: Risk-factor decomposition over time
# ============================================================================
cat("\n[3/6] Figure 3: risk-factor decomposition...\n")

# Use Tobacco (parent) + High alcohol use (alcohol's only level) for
# decomposition; this avoids double-counting Smoking + Chewing.
risk_d <- dr[location_name == "Germany" &
             rei_name %in% c("Tobacco","High alcohol use") &
             measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)") &
             metric_name == "Number" &
             age_name == "All ages" &
             sex_name %in% c("Male","Female")]

risk_d[, rei_name := factor(rei_name, levels = c("Tobacco","High alcohol use"))]
risk_d[, measure_lbl := ifelse(measure_name == "Deaths", "Attributable deaths",
                              "Attributable DALYs")]

# Also pull total burden to compute residual (non-attributable)
total_b <- d[location_name == "Germany" &
             measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)") &
             metric_name == "Number" &
             age_name == "All ages" &
             sex_name %in% c("Male","Female"),
             .(year, sex = sex_name, measure = measure_name, total = val)]

# Sum attributable across the two risks (note: risks for the same cause are NOT
# strictly additive due to overlapping mechanisms — but the GBD comparative
# risk assessment framework treats them as joint PAFs. For visualization we
# stack them as separate categories.)
risk_long <- risk_d[, .(year, sex = sex_name, measure = measure_name,
                       risk = rei_name, attributable = val)]
risk_pct  <- merge(risk_long, total_b, by = c("year","sex","measure"))
risk_pct[, paf := attributable / total * 100]
fwrite(risk_pct, file.path(DATA_OUT, "paf_germany_1990_2023.csv"))

# Plot: attributable deaths over time, stacked by risk, by sex (panel rows)
p3 <- ggplot(risk_long[measure == "Deaths"],
            aes(x = year, y = attributable, fill = risk)) +
  geom_area(alpha = 0.85, position = "stack") +
  facet_wrap(~ sex, ncol = 2) +
  scale_fill_manual(values = c(Tobacco = "#5b3a29", `High alcohol use` = "#c1462a")) +
  scale_x_continuous(breaks = seq(1990, 2020, 10)) +
  scale_y_continuous(labels = comma_format(accuracy = 1)) +
  labs(x = "Year", y = "Attributable deaths per year",
       fill = NULL,
       title = "Risk-factor attributable lip and oral cavity cancer deaths in Germany, 1990–2023",
       subtitle = "Annual deaths attributable to tobacco (smoking + chewing) and high alcohol use, by sex") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        plot.title.position = "plot",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "figure3_germany_risk_attributable.pdf"),
       p3, width = W_CM, height = 10, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figure3_germany_risk_attributable.png"),
       p3, width = W_CM, height = 10, units = "cm", dpi = PNG_DPI)
cat("  Figure 3 written.\n")

# Also report PAFs in 1990 vs 2023 for the Results text
paf_summary <- risk_pct[year %in% c(1990, 2023),
                        .(year, sex, measure, risk, attributable, paf = round(paf, 1))]
fwrite(paf_summary, file.path(TBL_DIR, "paf_1990_vs_2023.csv"))
cat("\nPAF summary (1990 vs 2023):\n")
print(paf_summary[measure == "Deaths"][order(sex, year, risk)])

# ============================================================================
# 4. Table 1: Master burden table
# ============================================================================
cat("\n[4/6] Table 1: master burden table...\n")

# Period totals: 1994-98, 2019-23 (observed); 2039-43 (projected)
key_periods <- c("1994-1998","2019-2023","2039-2043")

# From proj (Nordpred output already aggregated to periods)
proj_key <- proj[period %in% key_periods,
                 .(measure, sex, period, asr_per_100k, total_count)]

# Format
proj_key[, total_count := round(total_count)]
proj_key[, asr_per_100k := round(asr_per_100k, 2)]
proj_key[, measure_lbl := dplyr::recode(measure,
            "Deaths" = "Deaths",
            "DALYs (Disability-Adjusted Life Years)" = "DALYs",
            "Incidence" = "Incidence")]
proj_key[, measure_lbl := factor(measure_lbl,
                                 levels = c("Incidence","Deaths","DALYs"))]

# Wide format
t1_count <- dcast(proj_key, measure_lbl + sex ~ period, value.var = "total_count")
t1_asr   <- dcast(proj_key, measure_lbl + sex ~ period, value.var = "asr_per_100k")

# Merge EAPC
eapc_t1 <- eapc_all[measure %in% c("Deaths","DALYs (Disability-Adjusted Life Years)","Incidence"),
                    .(measure_lbl = dplyr::recode(measure,
                                                  "Deaths"="Deaths",
                                                  "DALYs (Disability-Adjusted Life Years)"="DALYs",
                                                  "Incidence"="Incidence"),
                      sex, eapc_fmt)]
t1_asr <- merge(t1_asr, eapc_t1, by = c("measure_lbl","sex"))

t1_count[, type := "Count (5-yr period total)"]
t1_asr[, type := "ASR per 100,000"]
t1 <- rbind(t1_count[, .(measure_lbl, sex, type,
                        `1994-1998` = `1994-1998`,
                        `2019-2023` = `2019-2023`,
                        `2039-2043` = `2039-2043`,
                        eapc_fmt = NA_character_)],
            t1_asr[, .(measure_lbl, sex, type,
                       `1994-1998` = `1994-1998`,
                       `2019-2023` = `2019-2023`,
                       `2039-2043` = `2039-2043`,
                       eapc_fmt)])

setorder(t1, measure_lbl, sex, type)
setnames(t1, "eapc_fmt", "EAPC 1990-2023 % (95% CI)")
fwrite(t1, file.path(TBL_DIR, "table1_master_burden.csv"))
cat("  Table 1 written. First 10 rows:\n")
print(head(t1, 12))

# ============================================================================
# 5. Comparison: Germany vs Western Europe vs Global, 2023
# ============================================================================
cat("\n[5/6] Germany vs Western Europe vs Global comparison (2023)...\n")

comp <- d[year == 2023 & age_name == "Age-standardized" & metric_name == "Rate" &
          sex_name == "Both" &
          measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)",
                              "Incidence","Prevalence")]
comp_w <- dcast(comp, measure_name ~ location_name, value.var = "val")
comp_w[, measure_name := dplyr::recode(measure_name,
        "Deaths"="Deaths","DALYs (Disability-Adjusted Life Years)"="DALYs",
        "Incidence"="Incidence","Prevalence"="Prevalence")]
comp_w[, Germany   := round(Germany, 2)]
comp_w[, `Western Europe` := round(`Western Europe`, 2)]
comp_w[, Global    := round(Global, 2)]
fwrite(comp_w, file.path(TBL_DIR, "table_germany_vs_we_vs_global_2023.csv"))
cat("  Comparison table:\n")
print(comp_w)

# ============================================================================
# 6. Figure 1: sex-stratified observed + projected ASR (Incidence, Deaths, DALYs)
#    Replaces the former both-sexes Figure 1 (reviewer: pooled-sex results are
#    uninformative given the sex divergence). Zero-based y-axes.
# ============================================================================
cat("\n[6/6] Figure 1: sex-stratified observed and projected ASR...\n")

MEAS3 <- c("Incidence","Deaths","DALYs (Disability-Adjusted Life Years)")
lbl3  <- function(x) factor(ifelse(x == "Incidence", "Incidence",
                            ifelse(x == "Deaths", "Deaths", "DALYs")),
                            levels = c("Incidence","Deaths","DALYs"))
obs_sx <- d[location_name == "Germany" &
           age_name == "Age-standardized" & metric_name == "Rate" &
           sex_name %in% c("Male","Female") & measure_name %in% MEAS3,
           .(sex = sex_name, measure = measure_name, year, asr = val,
             asr_lo = lower, asr_hi = upper)]
proj_sx <- proj[sex %in% c("Male","Female") & measure %in% MEAS3 &
                period_type == "forecast"]
obs_sx[,  measure_lbl := lbl3(measure)]
proj_sx[, measure_lbl := lbl3(measure)]

p1 <- ggplot() +
  geom_ribbon(data = obs_sx,
              aes(x = year, ymin = asr_lo, ymax = asr_hi, fill = sex),
              alpha = 0.2, colour = NA) +
  geom_line(data = obs_sx,
            aes(x = year, y = asr, colour = sex), linewidth = 0.7) +
  geom_line(data = proj_sx,
            aes(x = period_mid, y = asr_per_100k, colour = sex),
            linewidth = 0.7, linetype = "dashed") +
  geom_point(data = proj_sx,
             aes(x = period_mid, y = asr_per_100k, colour = sex), size = 1.5) +
  facet_wrap(~ measure_lbl, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c(Male = "#1f77b4", Female = "#d62728")) +
  scale_fill_manual(values = c(Male = "#1f77b4", Female = "#d62728")) +
  scale_x_continuous(breaks = seq(1990, 2040, 10)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06)),
                     labels = comma_format(accuracy = 1)) +
  labs(x = "Year", y = "Age-standardised rate per 100,000",
       colour = NULL, fill = NULL,
       title = "Lip and oral cavity cancer in Germany by sex: observed 1990–2023 and projected to 2040",
       subtitle = "GBD 2023 age-standardised rates with 95% UI (shaded); Nordpred projection at period midpoints (dashed).\ny-axes differ between panels.") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom",
        plot.title.position = "plot",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "figure1_germany_sex_stratified.pdf"),
       p1, width = W_CM, height = 8, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figure1_germany_sex_stratified.png"),
       p1, width = W_CM, height = 8, units = "cm", dpi = PNG_DPI)
ggsave(file.path(FIG_DIR, "figure1_germany_sex_stratified.tif"),   # submission-grade raster (cairo, 600 dpi)
       p1, width = W_CM, height = 8, units = "cm", dpi = 600, device = "tiff", compression = "lzw", type = "cairo")
cat("  Figure 1 (sex-stratified) written.\n")

cat("\nAll analyses complete.\n")
