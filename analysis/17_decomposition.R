# 17_decomposition.R
# Decompose the 1990 -> 2023 change in annual incident cases and deaths, by
# sex, into population growth, population ageing, and change in age-specific
# rates (Das Gupta-type three-factor decomposition, averaged over all factor
# orderings). Also split the change in deaths into tobacco-/alcohol-attributable
# and non-attributable components (GBD 2023 comparative risk assessment).

suppressPackageStartupMessages(library(data.table))
DATA_DIR <- "."  # repository root
CAUSE_CSV <- file.path(DATA_DIR, "data/gbd_cause", "IHME-GBD_2023_DATA-684738bc-1.csv")
RISK_CSV  <- file.path(DATA_DIR, "data/gbd_risk", "IHME-GBD_2023_DATA-6a8e77e9-1.csv")
TBL_DIR   <- file.path(DATA_DIR, "output", "tables")
AGE_RAW <- c("<5 years","5-9 years","10-14 years","15-19 years","20-24 years","25-29 years","30-34 years","35-39 years",
             "40-44 years","45-49 years","50-54 years","55-59 years","60-64 years","65-69 years","70-74 years",
             "75-79 years","80-84 years","85-89 years","90-94 years","95+ years")
d <- fread(CAUSE_CSV)
pn <- d[location_name == "Germany" & age_name %in% AGE_RAW & sex_name %in% c("Male","Female","Both") & metric_name == "Number", .(measure_name, sex = sex_name, year, age_name, n = val)]
pr <- d[location_name == "Germany" & age_name %in% AGE_RAW & sex_name %in% c("Male","Female","Both") & metric_name == "Rate",   .(measure_name, sex = sex_name, year, age_name, r = val)]
pop <- merge(pn, pr, by = c("measure_name","sex","year","age_name"))[r > 0][, .(pop = median(n / r * 1e5)), by = .(sex, year, age_name)]
cnt <- pn[measure_name %in% c("Incidence","Deaths"), .(measure = measure_name, sex, year, age_name, n)]
x <- merge(cnt, pop, by = c("sex","year","age_name"), all.x = TRUE)[year %in% c(1990, 2023)]
x[is.na(pop), pop := 0]

decomp <- function(dt) {   # dt: one measure x sex, years 1990 and 2023, columns age_name, n, pop
  y1 <- dt[year == 1990][order(match(age_name, AGE_RAW))]; y2 <- dt[year == 2023][order(match(age_name, AGE_RAW))]
  P1 <- sum(y1$pop); P2 <- sum(y2$pop); S1 <- y1$pop / P1; S2 <- y2$pop / P2
  R1 <- ifelse(y1$pop > 0, y1$n / y1$pop, 0); R2 <- ifelse(y2$pop > 0, y2$n / y2$pop, 0)
  N <- function(P, S, R) P * sum(S * R)
  f <- list(P = list(P1, P2), S = list(S1, S2), R = list(R1, R2)); eff <- c(P = 0, S = 0, R = 0)
  for (ord in list(c("P","S","R"), c("P","R","S"), c("S","P","R"), c("S","R","P"), c("R","P","S"), c("R","S","P"))) {
    cur <- list(P = P1, S = S1, R = R1); prev <- N(cur$P, cur$S, cur$R)
    for (k in ord) { cur[[k]] <- f[[k]][[2]]; nxt <- N(cur$P, cur$S, cur$R); eff[k] <- eff[k] + (nxt - prev) / 6; prev <- nxt } }
  data.table(n_1990 = sum(y1$n), n_2023 = sum(y2$n), total_change = sum(y2$n) - sum(y1$n),
             population_growth = eff[["P"]], population_ageing = eff[["S"]], rate_change = eff[["R"]]) }
res <- x[, decomp(.SD), by = .(measure, sex)]
res[, `:=`(pct_growth = population_growth / n_1990 * 100, pct_ageing = population_ageing / n_1990 * 100, pct_rate = rate_change / n_1990 * 100)]
fwrite(res, file.path(TBL_DIR, "decomposition_1990_2023.csv"))
cat("=== Decomposition of change in annual counts, 1990 -> 2023 (absolute; % of 1990 count) ===\n")
print(res[order(measure, sex), .(measure, sex, n_1990 = round(n_1990), n_2023 = round(n_2023),
      growth = sprintf("%+.0f (%+.0f%%)", population_growth, pct_growth), ageing = sprintf("%+.0f (%+.0f%%)", population_ageing, pct_ageing),
      rate = sprintf("%+.0f (%+.0f%%)", rate_change, pct_rate))])

# Attributable vs non-attributable change in deaths
dr <- fread(RISK_CSV)
att <- dr[location_name == "Germany" & age_name == "All ages" & metric_name == "Number" & measure_name == "Deaths" &
          rei_name %in% c("Tobacco","High alcohol use") & sex_name %in% c("Male","Female","Both") & year %in% c(1990, 2023),
          .(sex = sex_name, year, risk = rei_name, attributable = val)]
tot <- d[location_name == "Germany" & age_name == "All ages" & metric_name == "Number" & measure_name == "Deaths" & sex_name %in% c("Male","Female","Both") & year %in% c(1990, 2023), .(sex = sex_name, year, total = val)]
aw <- dcast(att, sex + year ~ risk, value.var = "attributable"); aw <- merge(aw, tot, by = c("sex","year"))
aw[, non_tobacco := total - Tobacco]
ch <- aw[, .(total_change = total[year == 2023] - total[year == 1990], tobacco_change = Tobacco[year == 2023] - Tobacco[year == 1990],
             alcohol_change = `High alcohol use`[year == 2023] - `High alcohol use`[year == 1990],
             non_tobacco_change = non_tobacco[year == 2023] - non_tobacco[year == 1990]), by = sex]
fwrite(ch, file.path(TBL_DIR, "attributable_change_deaths_1990_2023.csv"))
cat("\n=== Change in annual deaths 1990 -> 2023: tobacco-attributable vs non-tobacco (alcohol overlaps tobacco; shown separately) ===\n")
print(ch[, lapply(.SD, function(v) if (is.numeric(v)) round(v) else v)])
