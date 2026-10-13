#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/enriques_polarization_orbits.json from Dutour Sikirić--Hulek.

Source: arXiv:2302.01679, Table ListNumberPolarizationsModuli (vendored TeX). For each
degree 2d = 2, ..., 72 (genus g = d + 1) it gives #h, the number of orbits of
primitive vectors h with h^2 = 2d in Num(S) = U + E8(-1) (the paper notes these agree
with the appendix of Ciliberto--Dedieu--Galati--Knutsen), and #bar Gamma_h, the
number of conjugacy classes of the groups bar Gamma_h they define, together with the
cumulative counts over all degrees up to 2d.

The script checks the table against itself: 2d = 2g - 2 in every column, each
cumulative #h row is the running sum of #h, and #bar Gamma_h <= #h.

Run from the repository root: uv run references/extract/dh23_norm_orbits.py
"""

import json
from itertools import accumulate
from pathlib import Path

SOURCE = Path("references/vendor/arxiv/2302.01679/Enriques_compu_rev.tex")
TARGET = Path("tests/fixtures/enriques_polarization_orbits.json")
ROWS = {
    "$g$": "g",
    "$2d$": "two_d",
    r"$\# h$": "h_orbits",
    r"$\# \bar \Gamma_h$": "gamma_classes",
    r"$\# h, h^2 \leq 2d$": "h_orbits_up_to_2d",
    r"$\# \bar \Gamma_{h}, h^2\leq 2d$": "gamma_classes_up_to_2d",
}


def main() -> None:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    label = next(i for i, line in enumerate(lines) if r"\label{ListNumberPolarizationsModuli}" in line)
    begin = max(i for i in range(label) if r"\begin{tabular}" in lines[i])
    columns: dict[str, list[int]] = {key: [] for key in ROWS.values()}
    row_lines: dict[str, list[int]] = {key: [] for key in ROWS.values()}
    for index in range(begin + 1, label):
        line = lines[index].strip()
        if not line.endswith(r"\\"):
            continue
        cells = [cell.strip() for cell in line[:-2].split("&")]
        key = ROWS[cells[0]]
        values = [int(cell) for cell in cells[1:]]
        assert len(values) == 12, f"{SOURCE}:{index + 1}: expected 12 columns"
        columns[key] += values
        row_lines[key].append(index + 1)
    assert all(len(values) == 36 for values in columns.values()), {k: len(v) for k, v in columns.items()}
    assert columns["g"] == list(range(2, 38))
    assert columns["two_d"] == [2 * g - 2 for g in columns["g"]]
    assert columns["h_orbits_up_to_2d"] == list(accumulate(columns["h_orbits"])), "cumulative #h is not the running sum"
    assert all(gamma <= h for gamma, h in zip(columns["gamma_classes"], columns["h_orbits"], strict=True))
    records = [
        {
            "lattice": "U + E8(-1)",
            "two_d": columns["two_d"][i],
            "g": columns["g"][i],
            "primitive_vector_orbits": columns["h_orbits"][i],
            "gamma_conjugacy_classes": columns["gamma_classes"][i],
            "primitive_vector_orbits_up_to_2d": columns["h_orbits_up_to_2d"][i],
            "gamma_conjugacy_classes_up_to_2d": columns["gamma_classes_up_to_2d"][i],
            "source": {
                "kind": "published_table",
                "citation": "Dutour Sikirić--Hulek, arXiv:2302.01679",
                "file": str(SOURCE),
                "table": "ListNumberPolarizationsModuli",
                "lines": ",".join(str(row_lines[key][i // 12]) for key in ("two_d", "h_orbits", "gamma_classes")),
            },
        }
        for i in range(36)
    ]
    TARGET.write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {len(records)} records to {TARGET}")


if __name__ == "__main__":
    main()
