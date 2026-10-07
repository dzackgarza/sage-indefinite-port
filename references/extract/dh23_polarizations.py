#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/enriques_87_polarizations.json from Dutour Sikirić--Hulek.

Source: arXiv:2302.01679 ("Moduli of polarised Enriques surfaces -- computational
aspects", J. London Math. Soc. 2024), Tables table_subgroups1 and table_subgroups2,
read from the vendored TeX. Each row gives, for one of the 87 conjugacy classes of
groups Gamma_h: the subset S of the E10 diagram, #S, |bar Gamma_h|, the numbers of
isotropic lines #I_1, planes #I_2 and flags #I_12, the degree of the smallest
realization h_min, and phi(h_min).

Run from the repository root: uv run references/extract/dh23_polarizations.py
"""

import json
import math
import re
from pathlib import Path

SOURCE = Path("references/vendor/arxiv/2302.01679/Enriques_compu_rev.tex")
TARGET = Path("tests/fixtures/enriques_87_polarizations.json")
TABLE_LABELS = ("table_subgroups1", "table_subgroups2")
HEADER = r"Nr & $S$ & $\# S$ &  $|\bar \Gamma_h|$  & $\# I_1$  & $\# I_2$ & $\# I_{12}$ & $\min \deg$ & $\phi(h_{min})$\\"
FACTOR = re.compile(r"^(\d+)(?:\^\{(\d+)\})?$")


def group_order(cell: str) -> int:
    product = cell.strip().strip("$").replace(" ", "")
    factors = product.split(r"\cdot")
    values = []
    for factor in factors:
        match = FACTOR.match(factor)
        assert match, f"unparsed group-order factor {factor!r} in {cell!r}"
        values.append(int(match.group(1)) ** int(match.group(2) or 1))
    return math.prod(values)


def tabular_blocks(lines: list[str]) -> list[tuple[str, int, int]]:
    r"""Return (label, first line, last line) of each table's tabular body, 1-based."""
    blocks = []
    for label in TABLE_LABELS:
        label_line = next(i for i, line in enumerate(lines) if rf"\label{{{label}}}" in line)
        begin = max(i for i in range(label_line) if r"\begin{tabular}" in lines[i])
        end = min(i for i in range(begin, len(lines)) if r"\end{tabular}" in lines[i])
        assert end < label_line, f"{label}: tabular does not precede its label"
        header = [i for i in range(begin, end) if lines[i].strip() == HEADER]
        assert len(header) == 1, f"{label}: expected the column header once, found {len(header)}"
        blocks.append((label, header[0] + 2, end))
    return blocks


def main() -> None:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    records = []
    for label, first, last in tabular_blocks(lines):
        for index in range(first - 1, last):
            line = lines[index].strip()
            if line in ("", r"\hline"):
                continue
            assert line.endswith(r"\\"), f"{SOURCE}:{index + 1}: row does not end with \\\\: {line!r}"
            cells = [cell.strip() for cell in line[:-2].split("&")]
            assert len(cells) == 9, f"{SOURCE}:{index + 1}: expected 9 cells, found {len(cells)}"
            case, subset, orbit_length, order, lines_, planes, flags, degree, phi = cells
            records.append(
                {
                    "case": int(case),
                    "S": subset,
                    "polarization_orbit_length": int(orbit_length),
                    "group_order": group_order(order),
                    "line_orbits": int(lines_),
                    "plane_orbits": int(planes),
                    "flag_orbits": int(flags),
                    "minimal_degree": int(degree),
                    "phi_h_min": int(phi),
                    "source": {
                        "kind": "published_table",
                        "citation": "Dutour Sikirić--Hulek, arXiv:2302.01679",
                        "file": str(SOURCE),
                        "table": label,
                        "line": index + 1,
                    },
                }
            )
    cases = [record["case"] for record in records]
    assert cases == list(range(1, 88)), f"expected cases 1..87 in order, got {cases}"
    # Case 87 (S empty, trivial bar Gamma_h) is the stable group of N; DH section 3
    # states its isotropic lines and planes independently of the table: 528 and 24242.
    assert (records[-1]["line_orbits"], records[-1]["plane_orbits"]) == (528, 24242)
    TARGET.write_text(json.dumps(records, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {len(records)} records to {TARGET}")


if __name__ == "__main__":
    main()
