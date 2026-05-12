# Data inputs

This repository contains analysis code and derived outputs for a national
GBD 2023 analysis of lip and oral cavity cancer in Germany. The raw input
data are not redistributed here. The instructions below let you reproduce
the analysis end-to-end after a short manual download step.

## 1. GBD 2023 cause-of-death and risk-factor data

Source: [Global Burden of Disease Study 2023 (GBD 2023) Results](https://vizhub.healthdata.org/gbd-results/)
(Institute for Health Metrics and Evaluation).

**Cause query** — Save the resulting CSV as
`data/gbd_cause/IHME-GBD_2023_DATA-684738bc-1.csv`:

- Estimate: Cause of death or injury
- Cause: Lip and oral cavity cancer
- Location: Germany, Western Europe, Global
- Year: 1990–2023 (all)
- Age: all 5-year groups (Under 5 through 95+) plus All ages and Age-standardized
- Sex: Both, Male, Female
- Measure: Deaths, DALYs, YLLs, YLDs, Incidence, Prevalence
- Metric: Number, Rate, Percent

**Risk query** — Save as `data/gbd_risk/IHME-GBD_2023_DATA-6a8e77e9-1.csv`:

- Estimate: Risk factor
- Cause: Lip and oral cavity cancer
- Risks: Tobacco, Smoking, Chewing tobacco, High alcohol use (the four
  risks GBD 2023 attributes to lip and oral cavity cancer)
- Location: Germany, Western Europe, Global
- Same year/age/sex/metric dimensions as above
- Measure: Deaths, DALYs, YLLs, YLDs (incidence/prevalence not modelled for
  risk attribution)

For very wide queries you may need an account with a higher download
quota; the public limit is 100,000 rows per request, which can be reached
when including all risks across all years and ages. In that case run the
queries one location at a time or per measure.

## 2. ZfKD / Robert Koch-Institute cross-validation data (optional)

Source: [Centre for Cancer Registry Data (ZfKD)](https://www.krebsdaten.de/Krebs/EN/Database/databasequery_step1_node.html),
Robert Koch-Institute, Berlin.

Two queries (Inzidenz and Mortalität) with these settings:

- Kennzahl: Altersstandardisierte Rate
- Diagnose: C00–C08 (Lippe, Zungengrund, sonstige Teile der Zunge,
  Zahnfleisch, Mundboden, Gaumen, sonstige Teile des Mundes, Parotis,
  sonstige große Speicheldrüsen)
- Geschlecht: männlich + weiblich
- Jahre: 1999–2023
- Bevölkerungsstandard: Welt-Standard (Segi) plus Europa-Standard (alt, 1976)
- Layout: rows = Standard × Sex × Year, columns = Diagnose

Save the two CSVs as `data/zfkd/zfkd_incidence_c00_c08.csv` and
`data/zfkd/zfkd_mortality_c00_c08.csv`. Used only by `analysis/06_zfkd_comparison.R`.

## 3. UN World Population Prospects 2019

Source: the [`wpp2019`](https://cran.r-project.org/package=wpp2019) R package
on CRAN. Installed automatically by `analysis/00_setup.R`. Used by
`analysis/07_population_sensitivity.R`.

## Citations

When using the GBD 2023 inputs cite:

> Global Burden of Disease Collaborative Network. Global Burden of Disease
> Study 2023 (GBD 2023) Results. Seattle, United States: Institute for
> Health Metrics and Evaluation (IHME), 2024. Available from
> https://vizhub.healthdata.org/gbd-results/.

When using the ZfKD inputs cite:

> Zentrum für Krebsregisterdaten im Robert Koch-Institut: Datenbankabfrage
> mit Schätzung der Inzidenz, Prävalenz und des Überlebens von Krebs in
> Deutschland auf Basis der epidemiologischen Landeskrebsregisterdaten.
> Mortalitätsdaten bereitgestellt vom Statistischen Bundesamt.
> www.krebsdaten.de/abfrage.

When using the WPP 2019 inputs cite:

> United Nations, Department of Economic and Social Affairs, Population
> Division. World Population Prospects 2019. New York: UN, 2019.
