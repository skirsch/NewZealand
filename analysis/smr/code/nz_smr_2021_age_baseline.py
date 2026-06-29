#!/usr/bin/env python3
"""Independent age-only SMR check for the New Zealand Barry Young dataset.

Method:
- use one row per person, entering at last vaccination date;
- exit at death date if death occurs on/after last vaccination, otherwise censor date;
- censor at 2023-10-31 by default;
- split follow-up person-time across 5-year age bands using date_of_birth;
- multiply person-years by Stats NZ 2021 total-population age-specific death
  rates from Infoshare table DMM001AA, per 1,000 mean estimated population;
- compute SMR = observed deaths / expected deaths.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import defaultdict
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path


DATE_FORMAT = "%m-%d-%Y"


@dataclass
class Person:
    dob: date
    last_vax: date
    death: date | None
    last_dose_number: int | None = None


def parse_date(value: str) -> date | None:
    if not value:
        return None
    for fmt in (DATE_FORMAT, "%Y-%m-%d"):
        try:
            return datetime.strptime(value, fmt).date()
        except ValueError:
            pass
    raise ValueError(f"Unsupported date format: {value!r}")


def add_years(d: date, years: int) -> date:
    try:
        return d.replace(year=d.year + years)
    except ValueError:
        return d.replace(month=2, day=28, year=d.year + years)


def age_on(dob: date, when: date) -> int:
    return when.year - dob.year - ((when.month, when.day) < (dob.month, dob.day))


def band_for_age(age: int, bands: list[dict]) -> dict:
    for band in bands:
        if band["age_min"] <= age <= band["age_max"]:
            return band
    return bands[-1]


def next_band_boundary(dob: date, when: date, current_band: dict) -> date | None:
    next_age = current_band["age_max"] + 1
    if next_age > 200:
        return None
    boundary = add_years(dob, next_age)
    if boundary <= when:
        return None
    return boundary


def load_rates(path: Path) -> list[dict]:
    bands = []
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            bands.append(
                {
                    "age_band": row["age_band"],
                    "age_min": int(row["age_min"]),
                    "age_max": int(row["age_max"]),
                    "rate_per_py": float(row["death_rate_per_1000"]) / 1000.0,
                }
            )
    return bands


def load_people(path: Path) -> dict[str, Person]:
    people: dict[str, Person] = {}
    with gzip.open(path, "rt", newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f)
        for row in reader:
            mrn = row["mrn"]
            service = parse_date(row["date_time_of_service"])
            death = parse_date(row["date_of_death"])
            dob = parse_date(row["date_of_birth"])
            dose_number = int(row["dose_number"]) if row["dose_number"] else None
            if service is None or dob is None:
                continue
            existing = people.get(mrn)
            if existing is None:
                people[mrn] = Person(
                    dob=dob,
                    last_vax=service,
                    death=death,
                    last_dose_number=dose_number,
                )
                continue
            if service > existing.last_vax:
                existing.last_vax = service
                existing.dob = dob
                existing.last_dose_number = dose_number
            elif (
                service == existing.last_vax
                and dose_number is not None
                and (
                    existing.last_dose_number is None
                    or dose_number > existing.last_dose_number
                )
            ):
                existing.last_dose_number = dose_number
            if death and (existing.death is None or death < existing.death):
                existing.death = death
    return people


def compute_smr(people: dict[str, Person], bands: list[dict], censor: date):
    by_band = {
        band["age_band"]: {"person_years": 0.0, "expected": 0.0, "observed": 0}
        for band in bands
    }
    included = 0
    observed = 0
    expected = 0.0
    person_years = 0.0

    for person in people.values():
        start = person.last_vax
        if start > censor:
            continue
        end = censor
        died_in_followup = person.death is not None and start <= person.death <= censor
        if died_in_followup:
            end = person.death
            observed += 1

        included += 1
        cursor = start
        while cursor < end:
            age = age_on(person.dob, cursor)
            band = band_for_age(age, bands)
            boundary = next_band_boundary(person.dob, cursor, band)
            segment_end = min(end, boundary) if boundary else end
            py = (segment_end - cursor).days / 365.25
            exp = py * band["rate_per_py"]
            by_band[band["age_band"]]["person_years"] += py
            by_band[band["age_band"]]["expected"] += exp
            expected += exp
            person_years += py
            cursor = segment_end

        if died_in_followup:
            death_age = age_on(person.dob, person.death)  # type: ignore[arg-type]
            death_band = band_for_age(death_age, bands)["age_band"]
            by_band[death_band]["observed"] += 1

    smr = observed / expected
    # Large-count Poisson log CI for observed/expected.
    se_log = 1 / math.sqrt(observed)
    ci_low = math.exp(math.log(smr) - 1.96 * se_log)
    ci_high = math.exp(math.log(smr) + 1.96 * se_log)
    return {
        "included": included,
        "observed": observed,
        "expected": expected,
        "smr": smr,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "person_years": person_years,
        "by_band": by_band,
    }


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
        "--output-summary",
        default="analysis/smr/outputs/nz_smr_2021_age_baseline_summary.csv",
    )
    parser.add_argument(
        "--output-by-band",
        default="analysis/smr/outputs/nz_smr_2021_age_baseline_by_age_band.csv",
    )
    args = parser.parse_args()

    people = load_people(Path(args.records))
    bands = load_rates(Path(args.rates))
    result = compute_smr(people, bands, parse_date(args.censor_date) or date(2023, 10, 31))

    summary_path = Path(args.output_summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "included_people",
                "observed_deaths",
                "expected_deaths",
                "smr",
                "ci_low",
                "ci_high",
                "person_years",
                "censor_date",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "included_people": result["included"],
                "observed_deaths": result["observed"],
                "expected_deaths": f"{result['expected']:.6f}",
                "smr": f"{result['smr']:.6f}",
                "ci_low": f"{result['ci_low']:.6f}",
                "ci_high": f"{result['ci_high']:.6f}",
                "person_years": f"{result['person_years']:.6f}",
                "censor_date": args.censor_date,
            }
        )

    by_band_path = Path(args.output_by_band)
    with by_band_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["age_band", "person_years", "observed", "expected", "smr"]
        )
        writer.writeheader()
        for band in bands:
            row = result["by_band"][band["age_band"]]
            band_smr = row["observed"] / row["expected"] if row["expected"] else ""
            writer.writerow(
                {
                    "age_band": band["age_band"],
                    "person_years": f"{row['person_years']:.6f}",
                    "observed": row["observed"],
                    "expected": f"{row['expected']:.6f}",
                    "smr": f"{band_smr:.6f}" if band_smr != "" else "",
                }
            )

    print(f"included_people={result['included']}")
    print(f"observed_deaths={result['observed']}")
    print(f"expected_deaths={result['expected']:.2f}")
    print(f"smr={result['smr']:.4f}")
    print(f"95ci={result['ci_low']:.4f},{result['ci_high']:.4f}")
    print(f"person_years={result['person_years']:.2f}")
    print(f"summary={summary_path}")
    print(f"by_band={by_band_path}")


if __name__ == "__main__":
    main()
