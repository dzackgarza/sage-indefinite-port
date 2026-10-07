#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/e10_vector_orbits.json from Brandhorst--Gonzalez-Alonso.

Source: arXiv:2408.00306, Appendix Table table1 (vendored TeX). E10 is the even
unimodular lattice of signature (1,9) (line 289), i.e. U + E8(-1). For h^2 in
{0, 2, ..., 10} the table lists the O(E10)-orbits of primitive vectors h, by their
phi-invariant phi(h) = min{h.f : 0 != f in E10, f^2 = 0}, and the index of the image
of the stabilizer O(E10, h) in O(E10 tensor F2). The data was computed with OSCAR.

Cross-checks against independent sources, asserted here:
- the text's orbit counts 1, 2, 2, 2, 3 for h^2 = 2, ..., 10 match the table, and
  match Dutour Sikirić--Hulek's #h for 2d = 2, ..., 10 (enriques_polarization_orbits.json);
- the index for h^2 = 0 is 17 * 31 = 527, which Dutour Sikirić--Hulek state for the
  isotropic line orbit of length 527 (unpolarized_enriques.json).

Run from the repository root: uv run references/extract/bga_e10_orbits.py
"""

import json
import math
import re
from collections import Counter
from pathlib import Path
from typing import TypedDict

SOURCE = Path("references/vendor/arxiv/2408.00306/main_arxiv1.tex")
DH_ORBITS = Path("tests/fixtures/enriques_polarization_orbits.json")
TARGET = Path("tests/fixtures/e10_vector_orbits.json")
FACTOR = re.compile(r"^(\d+)(?:\^\{?(\d+)\}?)?$")


class Orbit(TypedDict):
    lattice: str
    h_squared: int
    phi: int | None
    stabilizer_image_index_in_O_E10_F2: int
    source: dict[str, str | int]


def product(cell: str) -> int:
    body = cell.strip().removeprefix(r"\(").removesuffix(r"\)").replace(" ", "")
    values = []
    for factor in body.split(r"\cdot"):
        match = FACTOR.match(factor)
        assert match, f"unparsed factor {factor!r} in {cell!r}"
        values.append(int(match.group(1)) ** int(match.group(2) or 1))
    return math.prod(values)


def main() -> None:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    label = next(i for i, line in enumerate(lines) if r"\label{table1}" in line)
    definition = next(i for i, line in enumerate(lines) if "even unimodular lattice of signature $(1,9)$ and therefore $E_{10}\\cong U \\oplus E_8$" in line)
    statement = next(i for i, line in enumerate(lines) if "There are 1,2,2,2,3 $O(E_{10})$-orbits of primitive vectors" in line)
    records: list[Orbit] = []
    for index in range(label + 1, len(lines)):
        line = lines[index].strip()
        if line.startswith(r"\end{tabular}"):
            break
        if not line.endswith(r"\\") or line.startswith(r"\(h^2"):
            continue
        square, phi, index_cell = (cell.strip() for cell in line[:-2].split("&"))
        records.append(
            {
                "lattice": "E10 = U + E8(-1)",
                "h_squared": int(square),
                "phi": None if phi == "-" else int(phi),
                "stabilizer_image_index_in_O_E10_F2": product(index_cell),
                "source": {"kind": "published_table", "file": str(SOURCE), "table": "table1", "line": index + 1},
            }
        )
    counts = Counter(record["h_squared"] for record in records if record["h_squared"] > 0)
    assert [counts[s] for s in (2, 4, 6, 8, 10)] == [1, 2, 2, 2, 3], counts
    dh = {record["two_d"]: record["primitive_vector_orbits"] for record in json.loads(DH_ORBITS.read_text())}
    assert all(dh[s] == counts[s] for s in (2, 4, 6, 8, 10)), "BGA and Dutour Sikirić--Hulek disagree on orbit counts"
    isotropic = [record for record in records if record["h_squared"] == 0]
    assert len(isotropic) == 1 and isotropic[0]["stabilizer_image_index_in_O_E10_F2"] == 527
    data = {
        "orbits": records,
        "source": {
            "definition": {"kind": "published_text", "file": str(SOURCE), "line": definition + 1},
            "orbit_counts": {"kind": "published_text", "file": str(SOURCE), "line": statement + 1},
        },
    }
    TARGET.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {len(records)} orbits to {TARGET}; counts agree with Dutour Sikirić--Hulek, isotropic index 527")


if __name__ == "__main__":
    main()
