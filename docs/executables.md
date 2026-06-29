# Executables

Run commands from the repository root unless noted otherwise.

## `analysis/smr/code/nz_smr_2021_age_baseline.py`

- Purpose: age-only all-cause mortality SMR after final recorded vaccination.
- Typical command: `python analysis\smr\code\nz_smr_2021_age_baseline.py`
- Inputs: `data/nz-record-level-data-4M-records.csv.gz`; `analysis/smr/data/statsnz_2021_age_specific_death_rates_total_population.csv`.
- Outputs: `analysis/smr/outputs/nz_smr_2021_age_baseline_summary.csv`; `analysis/smr/outputs/nz_smr_2021_age_baseline_by_age_band.csv`.

## `analysis/smr/code/nz_smr_by_last_dose_number.py`

- Purpose: age-only all-cause mortality SMR grouped by the final recorded dose number.
- Typical command: `python analysis\smr\code\nz_smr_by_last_dose_number.py`
- Inputs: `data/nz-record-level-data-4M-records.csv.gz`; `analysis/smr/data/statsnz_2021_age_specific_death_rates_total_population.csv`.
- Outputs: `analysis/smr/outputs/nz_smr_by_last_dose_number_2021_age_baseline.csv`; `analysis/smr/figures/nz_smr_by_last_dose_number_2021_age_baseline.png`.

## `analysis/smr/code/nz_time_since_last_dose_oe.py`

- Purpose: observed/expected mortality by days since final recorded dose using the same 2021 age-specific baseline.
- Typical command: `python analysis\smr\code\nz_time_since_last_dose_oe.py`
- Inputs: `data/nz-record-level-data-4M-records.csv.gz`; `analysis/smr/data/statsnz_2021_age_specific_death_rates_total_population.csv`.
- Outputs: `analysis/smr/outputs/nz_time_since_last_dose_oe_2021_age_baseline.csv`; `analysis/smr/figures/nz_time_since_last_dose_oe_2021_age_baseline.png`.

## `analysis/smr/code/nz_time_since_last_dose_calendar_adjusted.py`

- Purpose: diagnostic observed/expected mortality by days since final recorded dose and calendar quarter.
- Typical command: `python analysis\smr\code\nz_time_since_last_dose_calendar_adjusted.py`
- Inputs: `data/nz-record-level-data-4M-records.csv.gz`; `analysis/smr/data/statsnz_2021_age_specific_death_rates_total_population.csv`.
- Outputs: `analysis/smr/outputs/nz_time_since_last_dose_by_calendar_quarter_cells.csv`; `analysis/smr/outputs/nz_time_since_last_dose_by_calendar_quarter_summary.csv`; `analysis/smr/outputs/nz_calendar_quarter_oe_totals.csv`; `analysis/smr/figures/nz_time_since_last_dose_calendar_quarter_heatmap.png`.
