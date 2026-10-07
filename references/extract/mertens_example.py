#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/mertens_generators.json from Mertens arXiv:1303.3478, Example exComplete.

Source (vendored TeX references/vendor/arxiv/1303.3478/hyperbolic.tex): for the
hyperbolic Gram matrix A of determinant -155, Mertens's algorithm finds 9 inequivalent
D-perfect points with their neighbour counts, the stabilizers of those points (trivial
for x2, x5, x6, x7, x8, of order 2 otherwise) and 16 connecting elements. By the
paper's theorem (line 176) the stabilizers and connecting elements generate Omega, the
stabilizer of the positive cone in Aut_Z(A) = {g in GL_n(Z) : g A g^T = A} (line 129,
remark at line 268); hence Omega = O^+(L) and O(L) = <Omega, -I>.

The D-perfect points and neighbour counts are Mertens's construction; they are recorded
as such, not as this port's perfect-domain counts. The script asserts that every
stabilizer and connecting element is in GL_3(Z) and satisfies g A g^T = A; each
matrix is recorded with the paper's row-vector convention.

Run from the repository root: uv run references/extract/mertens_example.py
"""

import json
import re
from itertools import pairwise
from pathlib import Path
from typing import TypedDict

SOURCE = Path("references/vendor/arxiv/1303.3478/hyperbolic.tex")
TARGET = Path("tests/fixtures/mertens_generators.json")
PMATRIX = re.compile(r"\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", re.DOTALL)


class Generator(TypedDict):
    name: str
    matrix: list[list[int]]
    line: int


def parse(body: str) -> list[list[int]]:
    rows = [row for row in body.split(r"\\") if row.strip()]
    return [[int(cell) for cell in row.split("&")] for row in rows]


def multiply(a: list[list[int]], b: list[list[int]]) -> list[list[int]]:
    return [[sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a: list[list[int]]) -> list[list[int]]:
    return [list(row) for row in zip(*a, strict=True)]


def determinant3(m: list[list[int]]) -> int:
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])


def main() -> None:
    text = SOURCE.read_text(encoding="utf-8")
    start = text.index(r"\begin{example}\label{exComplete}")
    end = text.index(r"\end{example}", start)
    example = text[start:end]
    first_line = text[:start].count("\n") + 1
    gram_match = PMATRIX.search(example)
    assert gram_match, f"no pmatrix in {SOURCE} example exComplete at line {first_line}"
    gram = parse(gram_match.group(1))
    assert gram == [[-1, -3, -1], [-3, 14, 8], [-1, 8, 11]] and determinant3(gram) == -155
    # Sylvester: the leading principal minors 1, -1, -23, -155 change sign once, so one negative eigenvalue.
    minors = [1, gram[0][0], gram[0][0] * gram[1][1] - gram[0][1] * gram[1][0], determinant3(gram)]
    negative = sum(1 for a, b in pairwise(minors) if a * b < 0)
    assert all(minors) and negative == 1, minors
    assert "the algorithm finds $9$ inequivalent $D$-perfect points" in example
    points: dict[int, list[int]] = {}
    for label, body in re.findall(r"x_(\d)=\s*&\s*\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", example, re.DOTALL):
        points[int(label)] = [int(cell) for cell in body.split("&")]
    assert sorted(points) == list(range(1, 10))
    neighbours = [8, 4, 6, 8, 4, 3, 4, 3, 6]
    assert "$8,\\,\n4,\\,\n6,\\,\n8,\\,\n4,\\,\n3,\\,\n4,\\,\n3$ and $6$ neighbours" in example
    assert "The stabilizers of $x_2$, $x_5$, $x_6$, $x_7$, $x_8$ are trivial" in example
    generators: list[Generator] = []
    for match in re.finditer(r"(\\Stab\(x_(\d)\)|c_\{(\d),(\d)\}('?))=&\s*(?:\\langle\s*)?\\begin\{pmatrix\}(.*?)\\end\{pmatrix\}", example, re.DOTALL):
        matrix = parse(match.group(6))
        assert abs(determinant3(matrix)) == 1, match.group(1)
        assert multiply(multiply(matrix, gram), transpose(matrix)) == gram, f"{match.group(1)} is not an isometry"
        name = f"Stab(x_{match.group(2)})" if match.group(2) else f"c_{{{match.group(3)},{match.group(4)}}}{match.group(5)}"
        generators.append({"name": name, "matrix": matrix, "line": first_line + example[: match.start()].count("\n")})
    stabilizers = [g for g in generators if g["name"].startswith("Stab")]
    connecting = [g for g in generators if g["name"].startswith("c_")]
    assert sorted(g["name"] for g in stabilizers) == ["Stab(x_1)", "Stab(x_3)", "Stab(x_4)", "Stab(x_9)"], [g["name"] for g in stabilizers]
    assert len(connecting) == 16 and "we find $16$ connecting elements" in example, len(connecting)
    data = {
        "id": "mertens_det_minus_155",
        "gram": gram,
        "signature": [2, 1],
        "isometry_convention": "row vectors: g A g^T = A",
        "omega": "stabilizer of the positive cone in Aut_Z(A), i.e. O^+(L); O(L) = <Omega, -I>",
        "omega_generators": [{"name": g["name"], "matrix": g["matrix"]} for g in generators],
        "mertens_d_perfect_points": [{"point": points[k], "neighbours": neighbours[k - 1]} for k in range(1, 10)],
        "source": {"kind": "published_text", "file": str(SOURCE), "lines": f"{first_line}-{first_line + example.count(chr(10))}", "generation_theorem_line": 176},
    }
    TARGET.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}: {len(stabilizers)} stabilizers + {len(connecting)} connecting elements, all isometries of A")


if __name__ == "__main__":
    main()
