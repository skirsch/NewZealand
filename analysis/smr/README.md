# New Zealand SMR analysis

This directory contains an independent all-cause mortality SMR analysis using the Barry Young-derived New Zealand vaccination record file in this repository.

## Inputs

- Record-level vaccine data:
  - `data/nz-record-level-data-4M-records.csv.gz`
- Stats NZ age-specific death-rate baseline:
  - `analysis/smr/data/statsnz_2021_age_specific_death_rates_total_population.csv`
  - Source: Stats NZ Infoshare, `Death Rates - DMM`, `Age-specific death rates by sex, December years (total population) (Annual-Dec)`, table reference displayed as `DMM001AA`.
  - Selected dimensions: `Total Population`, 5-year age groups, year `2021`.
  - Units: deaths per 1,000 mean estimated population in each age group.

## Scripts

Run from the repository root:

```powershell
python analysis\smr\code\nz_smr_2021_age_baseline.py
python analysis\smr\code\nz_smr_by_last_dose_number.py
python analysis\smr\code\nz_time_since_last_dose_oe.py
python analysis\smr\code\nz_time_since_last_dose_calendar_adjusted.py
```

## Method summary

The main SMR script:

- collapses vaccination records to one row per `mrn`;
- defines cohort entry as each person's last vaccination date;
- censors follow-up at `2023-10-31`;
- counts deaths on or after last vaccination and on or before censoring;
- splits person-time across 5-year age bands as people age;
- computes expected deaths from Stats NZ 2021 total-population age-specific death rates;
- reports observed deaths, expected deaths, SMR, and a large-count Poisson log confidence interval.

The last-dose-number script uses the same cohort entry, censoring, age-band
splitting, and 2021 baseline, but groups people by the dose number attached to
their final recorded vaccination.

The time-since-last-dose script splits the same follow-up into windows:

- `0-7`
- `8-21`
- `22-42`
- `43-90`
- `91-180`
- `181-365`
- `366+`

The calendar-quarter script further splits person-time by calendar quarter as a diagnostic. It still uses the 2021 age-specific rates inside each quarter cell, so it is not a final true calendar-specific mortality baseline.

## Headline result

Using the 2021 age-only Stats NZ baseline:

- included people: `2,215,095`
- observed deaths: `37,038`
- expected deaths: `24,765.44`
- SMR: `1.4956`
- 95% CI: `1.4804` to `1.5109`

## Interpretation caution

This is an all-cause mortality signal after last recorded vaccination. It does not by itself distinguish vaccine-caused deaths from COVID mortality, frailty/selection, pandemic-era mortality, or other confounding. The time-since-last-dose analyses show non-flat mortality, but calendar-quarter diagnostics show strong calendar-period structure and small edge cells.

Best current use: a transparent anomaly / follow-up analysis, not a standalone vaccine-attribution estimate.
