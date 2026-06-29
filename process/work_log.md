# Work Log

## 2026-06-29 - SMR by final recorded dose number

### What we did

Added and ran an age-standardized all-cause mortality SMR analysis grouped by each person's final recorded dose number.

### Command / executable

- `python analysis\smr\code\nz_smr_by_last_dose_number.py`
- Regression check: `python analysis\smr\code\nz_smr_2021_age_baseline.py`

### Outputs

- `analysis/smr/outputs/nz_smr_by_last_dose_number_2021_age_baseline.csv`
- `analysis/smr/figures/nz_smr_by_last_dose_number_2021_age_baseline.png`
- Updated executable guide: `docs/executables.md`

### Results

Using the same 2021 age-only Stats NZ baseline and censor date `2023-10-31`, the well-populated final-dose groups did not show monotonically increasing SMRs. SMRs were approximately: dose 1 `2.8013`, dose 2 `1.9468`, dose 3 `1.7180`, dose 4 `1.6266`, dose 5 `0.6485`, dose 6 `1.3975`. Higher dose-number rows had very small expected death counts and unstable estimates.

The aggregate baseline regression check remained unchanged: included people `2,215,095`, observed deaths `37,038`, expected deaths `24,765.44`, SMR `1.4956`.

### Next steps

Consider a calendar-period-stratified version by final dose number, because final dose number is strongly entangled with rollout timing, age eligibility, and healthy-vaccinee selection.
