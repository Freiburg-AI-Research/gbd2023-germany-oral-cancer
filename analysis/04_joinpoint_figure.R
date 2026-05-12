# 04_joinpoint_figure.R
# Supplementary Figure: Joinpoint regression with segments, by sex × measure

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(dplyr); library(tidyr)
  library(scales); library(segmented)
})

DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause",
                      "IHME-GBD_2023_DATA-684738bc-1.csv")
TBL_DIR  <- file.path(DATA_DIR, "analysis", "output", "tables")
FIG_DIR  <- file.path(DATA_DIR, "analysis", "output", "figures")

W_CM <- 18; PNG_DPI <- 250

d <- fread(CAUSE_CSV)
jp_segs <- fread(file.path(TBL_DIR, "joinpoint_segments.csv"))

asr <- d[location_name == "Germany" &
         age_name == "Age-standardized" & metric_name == "Rate" &
         measure_name %in% c("Deaths","DALYs (Disability-Adjusted Life Years)",
                             "Incidence") &
         sex_name %in% c("Male","Female"),
         .(measure = measure_name, sex = sex_name, year, asr = val,
           asr_lo = lower, asr_hi = upper)]
asr[, sex := factor(sex, levels = c("Male","Female"))]
jp_segs <- jp_segs[sex %in% c("Male","Female")]
jp_segs[, sex := factor(sex, levels = c("Male","Female"))]
asr[, measure_lbl := dplyr::recode(measure,
    "Deaths" = "Deaths",
    "DALYs (Disability-Adjusted Life Years)" = "DALYs",
    "Incidence" = "Incidence")]
asr[, measure_lbl := factor(measure_lbl, levels = c("Incidence","Deaths","DALYs"))]

# For each (measure × sex) segment, generate a smooth fitted line
fit_lines <- jp_segs[, {
  # Predict fitted log-linear within each segment
  # APC = (exp(beta)-1)*100  ->  beta = log(1 + APC/100)
  beta <- log(1 + apc/100)
  # We need the intercept; refit log-linear on the segment
  seg_data <- asr[measure == .BY$measure & sex == .BY$sex &
                  year >= seg_start & year <= seg_end]
  if (nrow(seg_data) >= 2 && all(seg_data$asr > 0)) {
    lm_fit <- lm(log(asr) ~ year, data = seg_data)
    yrs <- seq(seg_start, seg_end)
    fitted <- exp(predict(lm_fit, newdata = data.table(year = yrs)))
    data.table(year = yrs, fitted = fitted, seg_id = paste(seg_start, seg_end, sep = "-"))
  } else {
    NULL
  }
}, by = .(measure, sex, seg_start, seg_end, apc, apc_fmt)]

fit_lines[, measure_lbl := dplyr::recode(measure,
    "Deaths" = "Deaths",
    "DALYs (Disability-Adjusted Life Years)" = "DALYs",
    "Incidence" = "Incidence")]
fit_lines[, measure_lbl := factor(measure_lbl, levels = c("Incidence","Deaths","DALYs"))]

# Breakpoint markers
bp <- unique(jp_segs[, .(measure, sex, breakpoint = seg_end)])
bp <- bp[breakpoint < 2023]  # exclude end of data
bp[, measure_lbl := dplyr::recode(measure,
    "Deaths" = "Deaths",
    "DALYs (Disability-Adjusted Life Years)" = "DALYs",
    "Incidence" = "Incidence")]
bp[, measure_lbl := factor(measure_lbl, levels = c("Incidence","Deaths","DALYs"))]

# Per-sex per-measure plot
pJP <- ggplot() +
  geom_point(data = asr,
             aes(x = year, y = asr, colour = sex), size = 0.9, alpha = 0.6) +
  geom_line(data = fit_lines,
            aes(x = year, y = fitted, colour = sex, group = interaction(sex, seg_id)),
            linewidth = 0.8) +
  geom_vline(data = bp,
             aes(xintercept = breakpoint),
             colour = "grey50", linetype = "dotted", linewidth = 0.4) +
  facet_grid(measure_lbl ~ sex, scales = "free_y") +
  scale_colour_manual(values = c(Male = "#1f77b4", Female = "#d62728")) +
  scale_x_continuous(breaks = seq(1990, 2020, 10)) +
  labs(x = "Year", y = "Age-standardized rate per 100,000",
       colour = NULL,
       title = "Joinpoint regression of lip and oral cavity cancer ASRs in Germany, 1990–2023",
       subtitle = "Solid lines: piecewise log-linear fit; dotted lines: detected change-points") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom",
        plot.title.position = "plot",
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(FIG_DIR, "figureS2_joinpoint.pdf"),
       pJP, width = W_CM, height = 14, units = "cm", device = "pdf")
ggsave(file.path(FIG_DIR, "figureS2_joinpoint.png"),
       pJP, width = W_CM, height = 14, units = "cm", dpi = PNG_DPI)

cat("Supplementary Figure 2 (joinpoint) written.\n")
