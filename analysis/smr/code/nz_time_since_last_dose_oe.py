#!/usr/bin/env python3
"""Time-since-last-dose O/E check for the New Zealand Barry Young dataset.

This is an internal causal-narrowing check. It uses the same 2021 Stats NZ
age-specific total-population death rates as nz_smr_2021_age_baseline.py, then
splits follow-up into time-since-last-dose windows.

This version adjusts for age and person-time, but not calendar month. Calendar
wave confounding remains an important interpretation limit.
"""

from __future__ import annotations

import argparse
import csv
import math
from datetime import date, timedelta
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


WINDOWS = [
    ("0-7", 0, 7),
    ("8-21", 8, 21),
    ("22-42", 22, 42),
    ("43-90", 43, 90),
    ("91-180", 91, 180),
    ("181-365", 181, 365),
    ("366+", 366, None),
]


def window_for_day(day: int) -> tuple[str, int, int | None]:
    for window in WINDOWS:
        _, start, end = window
        if day >= start and (end is None or day <= end):
            return window
    return WINDOWS[-1]


def next_window_boundary(last_vax: date, cursor: date) -> date | None:
    day = (cursor - last_vax).days
    _, _, end = window_for_day(day)
    if end is None:
        return None
    return last_vax + timedelta(days=end + 1)


def add_exposure(
    person: Person,
    start: date,
    end: date,
    bands: list[dict],
    rows: dict[str, dict],
) -> None:
    cursor = start
    while cursor < end:
        age = age_on(person.dob, cursor)
        band = band_for_age(age, bands)
        age_boundary = next_band_boundary(person.dob, cursor, band)
        win_boundary = next_window_boundary(person.last_vax, cursor)
        candidates = [end]
        if age_boundary:
            candidates.append(age_boundary)
        if win_boundary:
            candidates.append(win_boundary)
        segment_end = min(candidates)
        if segment_end <= cursor:
            raise RuntimeError(f"Non-advancing segment at {cursor} for {person}")

        py = (segment_end - cursor).days / 365.25
        label = window_for_day((cursor - person.last_vax).days)[0]
        rows[label]["person_years"] += py
        rows[label]["expected"] += py * band["rate_per_py"]
        cursor = segment_end


def compute_by_window(people: dict[str, Person], bands: list[dict], censor: date):
    rows = {
        label: {
            "window": label,
            "day_start": start,
            "day_end": end if end is not None else "",
            "person_years": 0.0,
            "observed": 0,
            "expected": 0.0,
        }
        for label, start, end in WINDOWS
    }

    included = 0
    for person in people.values():
        start = person.last_vax
        if start > censor:
            continue
        end = censor
        died = person.death is not None and start <= person.death <= censor
        if died:
            end = person.death  # type: ignore[assignment]
            death_day = (person.death - person.last_vax).days  # type: ignore[operator]
            rows[window_for_day(death_day)[0]]["observed"] += 1
        included += 1
        add_exposure(person, start, end, bands, rows)

    return included, rows


def add_smr_ci(row: dict) -> None:
    observed = row["observed"]
    expected = row["expected"]
    if not expected:
        row["oe"] = ""
        row["ci_low"] = ""
        row["ci_high"] = ""
        return
    oe = observed / expected
    row["oe"] = oe
    if observed == 0:
        row["ci_low"] = 0.0
        row["ci_high"] = 3.688879454 / expected
        return
    se_log = 1 / math.sqrt(observed)
    row["ci_low"] = math.exp(math.log(oe) - 1.96 * se_log)
    row["ci_high"] = math.exp(math.log(oe) + 1.96 * se_log)


def write_csv(path: Path, rows: dict[str, dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "window",
        "day_start",
        "day_end",
        "person_years",
        "observed",
        "expected",
        "oe",
        "ci_low",
        "ci_high",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for label, _, _ in WINDOWS:
            row = rows[label]
            writer.writerow(
                {
                    "window": row["window"],
                    "day_start": row["day_start"],
                    "day_end": row["day_end"],
                    "person_years": f"{row['person_years']:.6f}",
                    "observed": row["observed"],
                    "expected": f"{row['expected']:.6f}",
                    "oe": f"{row['oe']:.6f}" if row["oe"] != "" else "",
                    "ci_low": f"{row['ci_low']:.6f}" if row["ci_low"] != "" else "",
                    "ci_high": f"{row['ci_high']:.6f}" if row["ci_high"] != "" else "",
                }
            )


def plot_windows(path: Path, rows: dict[str, dict]) -> None:
    labels = [label for label, _, _ in WINDOWS]
    y = list(range(len(labels)))
    oe = [rows[label]["oe"] for label in labels]
    lo = [rows[label]["ci_low"] for label in labels]
    hi = [rows[label]["ci_high"] for label in labels]
    xerr = [[oe[i] - lo[i] for i in range(len(labels))], [hi[i] - oe[i] for i in range(len(labels))]]

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.errorbar(oe, y, xerr=xerr, fmt="o", color="black", ecolor="black", capsize=3)
    ax.axvline(1.0, color="#b22222", linestyle="--", linewidth=1.5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.set_xlabel("Observed / Expected deaths")
    ax.set_ylabel("Days since last vaccination")
    ax.set_title("New Zealand age-standardized mortality by time since last dose")
    ax.grid(axis="x", alpha=0.25)
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
        default="analysis/smr/outputs/nz_time_since_last_dose_oe_2021_age_baseline.csv",
    )
    parser.add_argument(
        "--plot",
        default="analysis/smr/figures/nz_time_since_last_dose_oe_2021_age_baseline.png",
    )
    args = parser.parse_args()

    people = load_people(Path(args.records))
    bands = load_rates(Path(args.rates))
    censor = parse_date(args.censor_date) or date(2023, 10, 31)
    included, rows = compute_by_window(people, bands, censor)
    for row in rows.values():
        add_smr_ci(row)

    write_csv(Path(args.output), rows)
    plot_windows(Path(args.plot), rows)

    print(f"included_people={included}")
    for label, _, _ in WINDOWS:
        row = rows[label]
        print(
            f"{label}: observed={row['observed']} expected={row['expected']:.2f} "
            f"oe={row['oe']:.4f} ci={row['ci_low']:.4f},{row['ci_high']:.4f} "
            f"py={row['person_years']:.2f}"
        )
    print(f"output={args.output}")
    print(f"plot={args.plot}")


if __name__ == "__main__":
    main()
