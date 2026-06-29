#!/usr/bin/env python3
"""Age-standardized mortality by final recorded dose number.

This uses the same cohort definition and 2021 Stats NZ age-specific baseline as
nz_smr_2021_age_baseline.py, but groups people by the dose number attached to
their final vaccination record.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from datetime import date
from pathlib import Path

import matplotlib.pyplot as plt

from nz_smr_2021_age_baseline import (
    Person,
    age_on,
    band_for_age,
    load_people,
    load_rates,
    next_band_boundary,
    parse_date,
)


def empty_row(last_dose_number: int | None) -> dict:
    return {
        "last_dose_number": last_dose_number if last_dose_number is not None else "",
        "included_people": 0,
        "person_years": 0.0,
        "observed": 0,
        "expected": 0.0,
    }


def add_exposure(
    person: Person,
    start: date,
    end: date,
    bands: list[dict],
    row: dict,
) -> None:
    cursor = start
    while cursor < end:
        age = age_on(person.dob, cursor)
        band = band_for_age(age, bands)
        boundary = next_band_boundary(person.dob, cursor, band)
        segment_end = min(end, boundary) if boundary else end
        if segment_end <= cursor:
            raise RuntimeError(f"Non-advancing segment at {cursor} for {person}")

        py = (segment_end - cursor).days / 365.25
        row["person_years"] += py
        row["expected"] += py * band["rate_per_py"]
        cursor = segment_end


def add_smr_ci(row: dict) -> None:
    observed = row["observed"]
    expected = row["expected"]
    if not expected:
        row["smr"] = ""
        row["ci_low"] = ""
        row["ci_high"] = ""
        return

    smr = observed / expected
    row["smr"] = smr
    if observed == 0:
        row["ci_low"] = 0.0
        row["ci_high"] = 3.688879454 / expected
        return

    se_log = 1 / math.sqrt(observed)
    row["ci_low"] = math.exp(math.log(smr) - 1.96 * se_log)
    row["ci_high"] = math.exp(math.log(smr) + 1.96 * se_log)


def compute_by_last_dose(people: dict[str, Person], bands: list[dict], censor: date):
    rows: dict[int | None, dict] = defaultdict(lambda: None)

    def get_row(last_dose_number: int | None) -> dict:
        if rows[last_dose_number] is None:
            rows[last_dose_number] = empty_row(last_dose_number)
        return rows[last_dose_number]

    included = 0
    for person in people.values():
        start = person.last_vax
        if start > censor:
            continue

        row = get_row(person.last_dose_number)
        row["included_people"] += 1
        included += 1

        end = censor
        died = person.death is not None and start <= person.death <= censor
        if died:
            end = person.death  # type: ignore[assignment]
            row["observed"] += 1

        add_exposure(person, start, end, bands, row)

    materialized = {k: v for k, v in rows.items() if v is not None}
    for row in materialized.values():
        add_smr_ci(row)
    return included, materialized


def dose_sort_key(value: int | None) -> tuple[int, int]:
    return (1, 0) if value is None else (0, value)


def write_csv(path: Path, rows: dict[int | None, dict], censor_date: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "last_dose_number",
        "included_people",
        "person_years",
        "observed_deaths",
        "expected_deaths",
        "smr",
        "ci_low",
        "ci_high",
        "censor_date",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for dose in sorted(rows, key=dose_sort_key):
            row = rows[dose]
            writer.writerow(
                {
                    "last_dose_number": row["last_dose_number"],
                    "included_people": row["included_people"],
                    "person_years": f"{row['person_years']:.6f}",
                    "observed_deaths": row["observed"],
                    "expected_deaths": f"{row['expected']:.6f}",
                    "smr": f"{row['smr']:.6f}" if row["smr"] != "" else "",
                    "ci_low": f"{row['ci_low']:.6f}" if row["ci_low"] != "" else "",
                    "ci_high": f"{row['ci_high']:.6f}" if row["ci_high"] != "" else "",
                    "censor_date": censor_date,
                }
            )


def plot_rows(path: Path, rows: dict[int | None, dict], min_expected: float) -> None:
    dose_numbers = [
        dose
        for dose in sorted(rows, key=dose_sort_key)
        if dose is not None and rows[dose]["expected"] >= min_expected
    ]
    labels = [str(dose) for dose in dose_numbers]
    smr = [rows[dose]["smr"] for dose in dose_numbers]
    lo = [rows[dose]["ci_low"] for dose in dose_numbers]
    hi = [rows[dose]["ci_high"] for dose in dose_numbers]
    x = list(range(len(labels)))
    yerr = [
        [smr[i] - lo[i] for i in range(len(labels))],
        [hi[i] - smr[i] for i in range(len(labels))],
    ]

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.errorbar(x, smr, yerr=yerr, fmt="o-", color="black", ecolor="black", capsize=3)
    ax.axhline(1.0, color="#b22222", linestyle="--", linewidth=1.5)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("Final recorded dose number")
    ax.set_ylabel("Observed / Expected deaths")
    ax.set_title(
        "New Zealand age-standardized mortality by final dose number"
        f"\nRows with expected deaths >= {min_expected:g}"
    )
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--records",
        default=r"data/nz-record-level-data-4M-records.csv.gz",
    )
    parser.add_argument(
        "--rates",
        default="analysis/smr/data/statsnz_2021_age_specific_death_rates_total_population.csv",
    )
    parser.add_argument("--censor-date", default="2023-10-31")
    parser.add_argument(
        "--output",
        default="analysis/smr/outputs/nz_smr_by_last_dose_number_2021_age_baseline.csv",
    )
    parser.add_argument(
        "--plot",
        default="analysis/smr/figures/nz_smr_by_last_dose_number_2021_age_baseline.png",
    )
    parser.add_argument(
        "--min-expected-for-plot",
        type=float,
        default=20.0,
        help="Hide very sparse dose-number rows from the plot only.",
    )
    args = parser.parse_args()

    people = load_people(Path(args.records))
    bands = load_rates(Path(args.rates))
    censor = parse_date(args.censor_date) or date(2023, 10, 31)
    included, rows = compute_by_last_dose(people, bands, censor)

    write_csv(Path(args.output), rows, args.censor_date)
    plot_rows(Path(args.plot), rows, args.min_expected_for_plot)

    print(f"included_people={included}")
    for dose in sorted(rows, key=dose_sort_key):
        row = rows[dose]
        label = row["last_dose_number"] if row["last_dose_number"] != "" else "missing"
        print(
            f"last_dose={label}: people={row['included_people']} "
            f"observed={row['observed']} expected={row['expected']:.2f} "
            f"smr={row['smr']:.4f} ci={row['ci_low']:.4f},{row['ci_high']:.4f} "
            f"py={row['person_years']:.2f}"
        )
    print(f"output={args.output}")
    print(f"plot={args.plot}")


if __name__ == "__main__":
    main()
