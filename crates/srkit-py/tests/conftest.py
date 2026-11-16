"""Shared fixtures: the physics oracle table, and a switch that hides NumPy."""

from __future__ import annotations

import csv
import math
import os
import sys
from pathlib import Path

import pytest

# tests/ -> srkit-py -> crates -> srkit -> repository root
REPO_ROOT = Path(__file__).resolve().parents[4]
ORACLE_CSV = Path(
    os.environ.get("SRKIT_ORACLE_CSV", REPO_ROOT / "tools" / "physics-checks" / "oracle.csv")
)


def _sig_figs(stated: str) -> int:
    mantissa = stated.lower().lstrip("+-").split("e")[0]
    digits = mantissa.replace(".", "").lstrip("0")
    if "." not in mantissa:
        digits = digits.rstrip("0") or "0"
    return max(len(digits), 1)


class OracleRow:
    """One row of `tools/physics-checks/oracle.csv`."""

    def __init__(self, record: dict[str, str]):
        self.id = record["id"]
        self.quantity = record["quantity"]
        self.stated_text = record["stated"].replace(",", "")
        self.stated = float(self.stated_text)
        self.unit = record["unit"]
        self.tolerance = float(record["tolerance"])

    def agrees(self, computed: float) -> bool:
        """The oracle's rule (as `oracle_check.py`): within the row's tolerance, or equal to the stated value
        once rounded to the significant figures it is stated to."""
        if not math.isfinite(computed):
            return False
        if self.stated == 0.0:
            return computed == 0.0
        if abs((computed - self.stated) / self.stated) <= self.tolerance:
            return True
        rounded = float(f"{computed:.{_sig_figs(self.stated_text) - 1}e}")
        return math.isclose(rounded, self.stated, rel_tol=1e-9)

    def __repr__(self) -> str:
        unit = f" {self.unit}" if self.unit else ""
        return f"{self.id} ({self.quantity}): stated {self.stated_text}{unit}, tolerance {self.tolerance}"


@pytest.fixture(scope="session")
def oracle() -> dict[str, OracleRow]:
    """The oracle table, keyed by row id."""
    if not ORACLE_CSV.is_file():
        pytest.fail(
            f"the physics oracle was not found at {ORACLE_CSV}: run the tests from a checkout of the "
            "repository, or point SRKIT_ORACLE_CSV at tools/physics-checks/oracle.csv"
        )
    with ORACLE_CSV.open(encoding="utf-8", newline="") as handle:
        rows = [OracleRow(record) for record in csv.DictReader(handle)]
    table = {row.id: row for row in rows}
    assert len(table) == len(rows), "duplicate oracle ids"
    return table


@pytest.fixture()
def no_numpy(monkeypatch):
    """Make `import numpy` fail for the length of a test, as on a machine without NumPy."""
    monkeypatch.setitem(sys.modules, "numpy", None)
