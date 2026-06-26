#!/usr/bin/env python3
"""Calendar-quarter adjusted time-since-last-dose O/E check.

This extends nz_time_since_last_dose_oe.py by splitting person-time across both
time-since-last-dose windows and calendar quarters. Expected deaths still use
the 2021 Stats NZ age-specific baseline, so the calendar adjustment here is
stratification/diagnostics, not a true quarter-specific mortality-rate baseline.
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
from nz_time_since_last_dose_oe import (
    WINDOWS,
    add_smr_ci,
    next_window_boundary,
    window_for_day,
)


def quarter_label(d: date) -> str:
    return f"{d.year}Q{((d.month - 1) // 3) + 1}"


def next_quarter_boundary(d: date) -> date:
    q = ((d.month - 1) // 3) + 1
    next_month = q * 3 + 1
    year = d.year
    if next_month == 13:
        next_month = 1
        year += 1
    return date(year, next_month, 1)


def empty_row(window: str, quarter: str) -> dict:
    return {
        "window": window,
        "calendar_quarter": quarter,
        "person_years": 0.0,
        "observed": 0,
        "expected": 0.0,
    }


def add_exposure(
    person: Person,
    start: date,
    end: date,
    bands: list[dict],
    cells: dict[tuple[str, str], dict],
) -> None:
    cursor = start
    while cursor < end:
        age = age_on(person.dob, cursor)
        band = band_for_age(age, bands)
        age_boundary = next_band_boundary(person.dob, cursor, band)
        win_boundary = next_window_boundary(person.last_vax, cursor)
        q_boundary = next_quarter_boundary(cursor)
        candidates = [end, q_boundary]
        if age_boundary:
            candidates.append(age_boundary)
        if win_boundary:
            candidates.append(win_boundary)
        segment_end = min(candidates)
        if segment_end <= cursor:
            raise RuntimeError(f"Non-advancing segment at {cursor} for {person}")

        window = window_for_day((cursor - person.last_vax).days)[0]
        quarter = quarter_label(cursor)
        cell = cells[(window, quarter)]
        py = (segment_end - cursor).days / 365.25
        cell["person_years"] += py
        cell["expected"] += py * band["rate_per_py"]
        cursor = segment_end


def compute_cells(people: dict[str, Person], bands: list[dict], censor: date):
    cells: dict[tuple[str, str], dict] = defaultdict(lambda: None)

    def get_cell(window: str, quarter: str) -> dict:
        key = (window, quarter)
        if cells[key] is None:
            cells[key] = empty_row(window, quarter)
        return cells[key]

    included = 0
    for person in people.values():
        start = person.last_vax
        if start > censor:
            continue
        end = censor
        died = person.death is not None and start <= person.death <= censor
        if died:
            end = person.death  # type: ignore[assignment]
            death_window = window_for_day((person.death - person.last_vax).days)[0]  # type: ignore[operator]
            death_quarter = quarter_label(person.death)  # type: ignore[arg-type]
            get_cell(death_window, death_quarter)["observed"] += 1
        included += 1
        # Materialize exposure cells lazily.
        cursor = start
        while cursor < end:
            window = window_for_day((cursor - person.last_vax).days)[0]
            quarter = quarter_label(cursor)
            get_cell(window, quarter)
            age = age_on(person.dob, cursor)
            band = band_for_age(age, bands)
            age_boundary = next_band_boundary(person.dob, cursor, band)
            win_boundary = next_window_boundary(person.last_vax, cursor)
            q_boundary = next_quarter_boundary(cursor)
            candidates = [end, q_boundary]
            if age_boundary:
                candidates.append(age_boundary)
            if win_boundary:
                candidates.append(win_boundary)
            segment_end = min(candidates)
            py = (segment_end - cursor).days / 365.25
            cells[(window, quarter)]["person_years"] += py
            cells[(window, quarter)]["expected"] += py * band["rate_per_py"]
            cursor = segment_end
    return included, {k: v for k, v in cells.items() if v is not None}


def summarize_by_window(cells: dict[tuple[str, str], dict]) -> dict[str, dict]:
    pooled = {
        label: {
            "window": label,
            "person_years": 0.0,
            "observed": 0,
            "expected": 0.0,
            "n_quarter_cells": 0,
            "min_expected_cell": "",
            "n_cells_expected_lt_5": 0,
            "n_cells_expected_lt_20": 0,
        }
        for label, _, _ in WINDOWS
    }
    for (window, _quarter), cell in cells.items():
        row = pooled[window]
        row["person_years"] += cell["person_years"]
        row["observed"] += cell["observed"]
        row["expected"] += cell["expected"]
        row["n_quarter_cells"] += 1
        exp = cell["expected"]
        row["min_expected_cell"] = (
            exp if row["min_expected_cell"] == "" else min(row["min_expected_cell"], exp)
        )
        row["n_cells_expected_lt_5"] += int(exp < 5)
        row["n_cells_expected_lt_20"] += int(exp < 20)
    for row in pooled.values():
        add_smr_ci(row)
    return pooled


def write_cells(path: Path, cells: dict[tuple[str, str], dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["window", "calendar_quarter", "person_years", "observed", "expected", "oe", "ci_low", "ci_high"]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for key in sorted(cells):
            row = cells[key]
            add_smr_ci(row)
            writer.writerow(
                {
                    "window": row["window"],
                    "calendar_quarter": row["calendar_quarter"],
                    "person_years": f"{row['person_years']:.6f}",
                    "observed": row["observed"],
                    "expected": f"{row['expected']:.6f}",
                    "oe": f"{row['oe']:.6f}" if row["oe"] != "" else "",
                    "ci_low": f"{row['ci_low']:.6f}" if row["ci_low"] != "" else "",
                    "ci_high": f"{row['ci_high']:.6f}" if row["ci_high"] != "" else "",
                }
            )


def write_summary(path: Path, pooled: dict[str, dict]) -> None:
    fields = [
        "window",
        "person_years",
        "observed",
        "expected",
        "oe",
        "ci_low",
        "ci_high",
        "n_quarter_cells",
        "min_expected_cell",
        "n_cells_expected_lt_5",
        "n_cells_expected_lt_20",
    ]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for label, _, _ in WINDOWS:
            row = pooled[label]
            writer.writerow(
                {
                    "window": row["window"],
                    "person_years": f"{row['person_years']:.6f}",
                    "observed": row["observed"],
                    "expected": f"{row['expected']:.6f}",
                    "oe": f"{row['oe']:.6f}" if row["oe"] != "" else "",
                    "ci_low": f"{row['ci_low']:.6f}" if row["ci_low"] != "" else "",
                    "ci_high": f"{row['ci_high']:.6f}" if row["ci_high"] != "" else "",
                    "n_quarter_cells": row["n_quarter_cells"],
                    "min_expected_cell": f"{row['min_expected_cell']:.6f}" if row["min_expected_cell"] != "" else "",
                    "n_cells_expected_lt_5": row["n_cells_expected_lt_5"],
                    "n_cells_expected_lt_20": row["n_cells_expected_lt_20"],
                }
            )


def write_quarter_totals(path: Path, cells: dict[tuple[str, str], dict]) -> None:
    totals: dict[str, dict] = defaultdict(lambda: {"person_years": 0.0, "observed": 0, "expected": 0.0})
    for (_window, quarter), cell in cells.items():
        totals[quarter]["person_years"] += cell["person_years"]
        totals[quarter]["observed"] += cell["observed"]
        totals[quarter]["expected"] += cell["expected"]
    with path.open("w", newline="", encoding="utf-8") as f:
        fields = ["calendar_quarter", "person_years", "observed", "expected", "oe"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for quarter in sorted(totals):
            row = totals[quarter]
            writer.writerow(
                {
                    "calendar_quarter": quarter,
                    "person_years": f"{row['person_years']:.6f}",
                    "observed": row["observed"],
                    "expected": f"{row['expected']:.6f}",
                    "oe": f"{(row['observed'] / row['expected']):.6f}" if row["expected"] else "",
                }
            )


def plot_heatmap(path: Path, cells: dict[tuple[str, str], dict]) -> None:
    windows = [label for label, _, _ in WINDOWS]
    quarters = sorted({q for (_w, q) in cells})
    matrix = []
    for window in windows:
        row = []
        for quarter in quarters:
            cell = cells.get((window, quarter))
            row.append((cell["observed"] / cell["expected"]) if cell and cell["expected"] else float("nan"))
        matrix.append(row)

    fig, ax = plt.subplots(figsize=(max(8, len(quarters) * 0.55), 4.6))
    im = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=0.5, vmax=2.5)
    ax.set_xticks(range(len(quarters)))
    ax.set_xticklabels(quarters, rotation=45, ha="right")
    ax.set_yticks(range(len(windows)))
    ax.set_yticklabels(windows)
    ax.set_xlabel("Calendar quarter")
    ax.set_ylabel("Days since last dose")
    ax.set_title("New Zealand O/E by time since last dose and calendar quarter")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Observed / Expected")
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
        "--cells-output",
        default="analysis/smr/outputs/nz_time_since_last_dose_by_calendar_quarter_cells.csv",
    )
    parser.add_argument(
        "--summary-output",
        default="analysis/smr/outputs/nz_time_since_last_dose_by_calendar_quarter_summary.csv",
    )
    parser.add_argument(
        "--quarter-output",
        default="analysis/smr/outputs/nz_calendar_quarter_oe_totals.csv",
    )
    parser.add_argument(
        "--heatmap",
        default="analysis/smr/figures/nz_time_since_last_dose_calendar_quarter_heatmap.png",
    )
    args = parser.parse_args()

    people = load_people(Path(args.records))
    bands = load_rates(Path(args.rates))
    censor = parse_date(args.censor_date) or date(2023, 10, 31)
    included, cells = compute_cells(people, bands, censor)
    pooled = summarize_by_window(cells)

    write_cells(Path(args.cells_output), cells)
    write_summary(Path(args.summary_output), pooled)
    write_quarter_totals(Path(args.quarter_output), cells)
    plot_heatmap(Path(args.heatmap), cells)

    print(f"included_people={included}")
    print(f"calendar_quarter_cells={len(cells)}")
    for label, _, _ in WINDOWS:
        row = pooled[label]
        print(
            f"{label}: observed={row['observed']} expected={row['expected']:.2f} "
            f"oe={row['oe']:.4f} cells={row['n_quarter_cells']} "
            f"min_expected={row['min_expected_cell']:.2f} "
            f"cells_exp_lt_20={row['n_cells_expected_lt_20']}"
        )
    print(f"cells_output={args.cells_output}")
    print(f"summary_output={args.summary_output}")
    print(f"quarter_output={args.quarter_output}")
    print(f"heatmap={args.heatmap}")


if __name__ == "__main__":
    main()
