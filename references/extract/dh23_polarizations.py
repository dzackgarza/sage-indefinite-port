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

S is the set of walls of the Weyl chamber H (labelled by the simple roots -1..8) that
contain h; a line over the set means its complement (lines 686-688), and the proof of
Theorem teo:87groups shows Gamma_h^+ depends only on S. #S counts the subsets S that
define the same conjugacy class. For each case the script records S as root labels and
the vector h_S = sum of the extreme rays g_i whose one non-containing wall is not in S,
using the published rays and roots (e10_fundamental_domain.json, run
dh23_e10_domain.py first). It asserts that the chamber is simplicial (each ray lies off
exactly one wall), that h_S is orthogonal to exactly the roots in S, that no h_S is
shorter than the table's h_min, and that h_S for S empty is the norm-1240 Weyl vector.

Run from the repository root: uv run references/extract/dh23_polarizations.py
"""

import json
import math
import re
from collections.abc import Callable
from pathlib import Path
from typing import TypedDict

SOURCE = Path("references/vendor/arxiv/2302.01679/Enriques_compu_rev.tex")
TARGET = Path("tests/fixtures/enriques_87_polarizations.json")
TABLE_LABELS = ("table_subgroups1", "table_subgroups2")
HEADER = r"Nr & $S$ & $\# S$ &  $|\bar \Gamma_h|$  & $\# I_1$  & $\# I_2$ & $\# I_{12}$ & $\min \deg$ & $\phi(h_{min})$\\"
FACTOR = re.compile(r"^(\d+)(?:\^\{(\d+)\})?$")
E10_DOMAIN = Path("tests/fixtures/e10_fundamental_domain.json")
ROOT_LABELS = tuple(range(-1, 9))


class Polarization(TypedDict):
    case: int
    S: str
    walls: list[int]
    subsets_in_class: int
    face_polarization: list[int]
    face_polarization_norm: int
    group_order: int
    line_orbits: int
    plane_orbits: int
    flag_orbits: int
    minimal_degree: int
    phi_h_min: int
    source: dict[str, str | int]


def walls(cell: str) -> list[int]:
    """The root labels of S, resolving the complement bar."""
    body = cell.strip().strip("$")
    if body == r"\emptyset":
        return []
    match = re.fullmatch(r"(\\overline\{)?\\\{([-\d,]+)\\\}(\})?", body)
    assert match and bool(match.group(1)) == bool(match.group(3)), f"unparsed S cell {cell!r}"
    listed = sorted(int(label) for label in match.group(2).split(","))
    assert set(listed) <= set(ROOT_LABELS), cell
    return sorted(set(ROOT_LABELS) - set(listed)) if match.group(1) else listed


def face_polarizations() -> tuple[dict[int, list[int]], Callable[[list[int]], tuple[list[int], int]]]:
    """Map each wall label to the extreme ray lying off it, and return h_S as a function of S."""
    domain = json.loads(E10_DOMAIN.read_text(encoding="utf-8"))
    gram, rays = domain["gram"], domain["extreme_rays"]
    roots = {root["label"]: root["vector"] for root in domain["simple_roots"]}

    def pairing(u: list[int], v: list[int]) -> int:
        return sum(u[i] * gram[i][j] * v[j] for i in range(10) for j in range(10))

    off_wall: dict[int, list[int]] = {}
    for ray in rays:
        off = [label for label in ROOT_LABELS if pairing(ray, roots[label]) != 0]
        assert len(off) == 1, f"ray {ray} lies off walls {off}: the chamber is not simplicial"
        off_wall[off[0]] = ray
    assert sorted(off_wall) == list(ROOT_LABELS)

    def polarization(subset: list[int]) -> tuple[list[int], int]:
        h = [sum(off_wall[label][j] for label in ROOT_LABELS if label not in subset) for j in range(10)]
        orthogonal = [label for label in ROOT_LABELS if pairing(h, roots[label]) == 0]
        assert orthogonal == subset, f"h_S for S = {subset} is orthogonal to {orthogonal}"
        return h, pairing(h, h)

    return off_wall, polarization


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
    _off_wall, polarization = face_polarizations()
    records: list[Polarization] = []
    for label, first, last in tabular_blocks(lines):
        for index in range(first - 1, last):
            line = lines[index].strip()
            if line in ("", r"\hline"):
                continue
            assert line.endswith(r"\\"), f"{SOURCE}:{index + 1}: row does not end with \\\\: {line!r}"
            cells = [cell.strip() for cell in line[:-2].split("&")]
            assert len(cells) == 9, f"{SOURCE}:{index + 1}: expected 9 cells, found {len(cells)}"
            case, subset, subsets_in_class, order, lines_, planes, flags, degree, phi = cells
            h, norm = polarization(walls(subset))
            records.append(
                {
                    "case": int(case),
                    "S": subset,
                    "walls": walls(subset),
                    "subsets_in_class": int(subsets_in_class),
                    "face_polarization": h,
                    "face_polarization_norm": norm,
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
    # h_min has minimal norm among the polarizations defining Gamma_h, so no h_S is shorter;
    # for S empty, h_S is the Weyl vector, of norm 1240 (line weyl_vector of the E10 domain).
    assert all(record["face_polarization_norm"] >= record["minimal_degree"] for record in records)
    assert records[-1]["walls"] == [] and records[-1]["face_polarization_norm"] == 1240
    TARGET.write_text(json.dumps(records, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {len(records)} records to {TARGET}")


if __name__ == "__main__":
    main()
