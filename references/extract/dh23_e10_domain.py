#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/e10_fundamental_domain.json from Dutour Sikirić--Hulek section 3.

Source: arXiv:2302.01679, subsections "Some basic facts and roots" and "Enumerating
polarizations of small degree" (vendored TeX). For U + E8(-1) with the Gram matrix G
of equation equ:G the paper gives the ten simple roots r_{-1}, ..., r_8 (r^2 = -2) of
a Weyl chamber, the extreme rays g_1, ..., g_10 of that chamber (a Z-basis) with their
Gram matrix W (equation equ:W), states O(U + E8(-1)) = W(U + E8(-1)) and that the Weyl
vector has norm 1240, and lists representatives of every orbit of primitive vectors of
norm at most 30 (Table ListVectorsNormAtMost30), grouped by the conjugacy class of
bar Gamma_h they define.

Checks asserted here: W == g G g^T; every simple root has norm -2 and distinct simple
roots pair non-negatively (the Coxeter condition in the norm -2 convention); every
listed representative has the norm of its degree; and the number of
representatives in each degree equals #h of Table ListNumberPolarizationsModuli
(enriques_polarization_orbits.json).

Run from the repository root: uv run references/extract/dh23_e10_domain.py
"""

import json
import re
from pathlib import Path
from typing import TypedDict

SOURCE = Path("references/vendor/arxiv/2302.01679/Enriques_compu_rev.tex")
DH_ORBITS = Path("tests/fixtures/enriques_polarization_orbits.json")
TARGET = Path("tests/fixtures/e10_fundamental_domain.json")


class Root(TypedDict):
    label: int
    vector: list[int]


class Representative(TypedDict):
    degree: int
    gamma_class: int
    g_coefficients: list[int]
    vector: list[int]
    source_line: int


def matrix_after(lines: list[str], start: int, size: int) -> tuple[list[list[int]], int]:
    """Read the first size x size array of integer rows at or after line index start."""
    rows: list[list[int]] = []
    index = start
    while len(rows) < size:
        text = lines[index].strip().removesuffix(r"\\").strip()
        index += 1
        if re.fullmatch(r"-?\d+(\s*&\s*-?\d+)+", text):
            rows.append([int(cell) for cell in text.split("&")])
    return rows, index


def pairing(gram: list[list[int]], u: list[int], v: list[int]) -> int:
    return sum(u[i] * gram[i][j] * v[j] for i in range(len(u)) for j in range(len(v)))


def main() -> None:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    g_line = next(i for i, line in enumerate(lines) if r"\label{equ:G}" in line)
    gram, _ = matrix_after(lines, g_line, 10)
    roots_line = next(i for i, line in enumerate(lines) if "are numbered from $-1$ to $8$ and have the coordinates" in line)
    roots: list[Root] = []
    for index in range(roots_line + 1, roots_line + 20):
        match = re.match(r"^\s*(-?\d+) &=& \(([-\d, ]+)\)", lines[index])
        if match:
            roots.append({"label": int(match.group(1)), "vector": [int(x) for x in match.group(2).split(",")]})
    assert [root["label"] for root in roots] == list(range(-1, 9)), roots
    rays_line = next(i for i, line in enumerate(lines) if "The generators $g_i$ of the extreme rays are the following" in line)
    rays, _ = matrix_after(lines, rays_line, 10)
    w_line = next(i for i, line in enumerate(lines) if r"\label{equ:W}" in line)
    w, _ = matrix_after(lines, w_line, 10)
    assert all(pairing(gram, r["vector"], r["vector"]) == -2 for r in roots), "a simple root does not have norm -2"
    assert all(pairing(gram, a["vector"], b["vector"]) >= 0 for a in roots for b in roots if a is not b), (
        "simple roots must pair non-negatively off the diagonal (norm -2 convention)"
    )
    assert [[pairing(gram, a, b) for b in rays] for a in rays] == w, "W != g G g^T"
    assert "the isometry group and the Coxeter group coincide" in "\n".join(lines[w_line : w_line + 20])
    weyl_line = next(i for i, line in enumerate(lines) if "The minimal norm of integer vectors with trivial stabilizer in $U + E_8(-1)$ is $1240$" in line)

    table_label = next(i for i, line in enumerate(lines) if r"\label{ListVectorsNormAtMost30}" in line)
    table_begin = max(i for i in range(table_label) if r"\begin{tabular}" in lines[i])
    term = re.compile(r"(\d*)g_\{(\d+)\}")
    expected = {r["two_d"]: r["primitive_vector_orbits"] for r in json.loads(DH_ORBITS.read_text())}
    rows: list[tuple[int, str]] = []
    for index in range(table_begin, table_label):
        line = lines[index].strip()
        if re.match(r"^[^&\s]+ & ", line):
            rows.append((index, line))
        elif line.startswith("&") and rows:
            rows[-1] = (rows[-1][0], rows[-1][1].removesuffix("\\\\").rstrip() + " " + line.lstrip("&").strip())
    representatives: list[Representative] = []
    for index, line in rows:
        if line.startswith(("all &", "$\\deg$ &")):
            continue
        match = re.match(r"^(\d+) & (.*)\\\\$", line)
        assert match, f"{SOURCE}:{index + 1}: unparsed row {line!r}"
        degree = int(match.group(1))
        groups = re.findall(r"(\d+) : \$\\\{(.*?)\\\}\$", match.group(2).replace("$, $", ", "))
        count = 0
        for gamma_class, body in groups:
            for expression in body.split(","):
                coefficients = [0] * 10
                for coefficient, index_g in term.findall(expression.replace(" ", "")):
                    coefficients[int(index_g) - 1] += int(coefficient or 1)
                vector = [sum(coefficients[k] * rays[k][j] for k in range(10)) for j in range(10)]
                assert pairing(gram, vector, vector) == degree, f"{SOURCE}:{index + 1}: {expression} has norm {pairing(gram, vector, vector)}"
                representatives.append({"degree": degree, "gamma_class": int(gamma_class), "g_coefficients": coefficients, "vector": vector, "source_line": index + 1})
                count += 1
        assert count == expected[degree], f"degree {degree}: {count} representatives, #h = {expected[degree]}"
    degrees = sorted({r["degree"] for r in representatives})
    assert degrees == list(range(2, 31, 2)), degrees

    def at(index: int) -> dict[str, str | int]:
        return {"kind": "published_text", "file": str(SOURCE), "line": index + 1}

    data = {
        "lattice": "U + E8(-1)",
        "gram": gram,
        "simple_roots": roots,
        "extreme_rays": rays,
        "extreme_ray_gram": w,
        "orthogonal_group_equals_weyl_group": True,
        "weyl_vector_norm": 1240,
        "orbit_representatives": representatives,
        "source": {
            "gram": at(g_line),
            "simple_roots": at(roots_line),
            "extreme_rays": at(rays_line),
            "extreme_ray_gram": at(w_line),
            "weyl_vector": at(weyl_line),
            "representatives": {"kind": "published_table", "file": str(SOURCE), "table": "ListVectorsNormAtMost30"},
        },
    }
    TARGET.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}: 10 roots, W == g G g^T, {len(representatives)} representatives for degrees 2..30, all matching #h")


if __name__ == "__main__":
    main()
