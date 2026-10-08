#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/conway_sloane_cases.json from Conway--Sloane, SPLAG 3rd ed.

Source: Chapter 15, Section 11 ("Computational complexity"), the two inequivalent
indefinite ternary forms (51a), (51b) of determinant -128 in the genus of
diag{-1, 64, 2}. The text is the OCR extraction of Zotero item T2WVLTDB, pinned by
its SHA-256; the file is a copyrighted book and is not committed (see
references/vendor/README.md).

The two Gram matrices are parsed from the displayed equation. The script asserts
what the passage states about them: both have determinant -128 and signature
(2, 1); (51a) has the diagonal basis e1, e1 + e2, e3 with Gram diag{-1, 64, 2};
and (51b) is the Gram matrix of 3 e1, e2 / 3, e3 in the rational span of (51a).

Run from the repository root: uv run references/extract/conway_sloane.py
"""

import hashlib
import json
import re
from fractions import Fraction
from pathlib import Path

SOURCE = Path("references/vendor/zotero/T2WVLTDB-conway-sloane-1999-splag/extracted.md")
SOURCE_SHA256 = "75827afb551624d82a01494977c4335549279897c6879890ded080c171012ec7"
TARGET = Path("tests/fixtures/conway_sloane_cases.json")
CITATION = "J. H. Conway and N. J. A. Sloane, Sphere Packings, Lattices and Groups, 3rd ed., Springer, 1999, Ch. 15, Sec. 11"

Matrix = list[list[Fraction]]


def read_source() -> list[str]:
    data = SOURCE.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    assert digest == SOURCE_SHA256, f"{SOURCE} has SHA-256 {digest}, not the pinned {SOURCE_SHA256}"
    return data.decode("utf-8").splitlines()


def cite(lines: list[str], first: int, last: int, *needles: str) -> dict[str, str]:
    text = "\n".join(lines[first - 1 : last])
    for needle in needles:
        assert needle in text, f"{SOURCE}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "citation": CITATION, "file": str(SOURCE), "lines": f"{first}-{last}"}


def displayed_arrays(line: str) -> list[list[list[int]]]:
    """The integer matrices of every LaTeX array on one OCR line, in order."""
    arrays = []
    for body in re.findall(r"\\begin\{array\}\s*\{[ l]*\}(.*?)\\end\{array\}", line):
        rows = []
        for row in body.split(r"\\"):
            cells = [re.sub(r"[{}\s]", "", cell) for cell in row.split("&")]
            rows.append([int(cell) for cell in cells])
        arrays.append(rows)
    return arrays


def product(left: Matrix, right: Matrix) -> Matrix:
    return [[sum((left[i][k] * right[k][j] for k in range(len(right))), Fraction(0)) for j in range(len(right[0]))] for i in range(len(left))]


def transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix, strict=True)]


def congruent(gram: list[list[int]], basis: list[list[Fraction]]) -> Matrix:
    """The Gram matrix of the given basis rows, in the coordinates of gram."""
    exact = [[Fraction(entry) for entry in row] for row in gram]
    return product(product(basis, exact), transpose(basis))


def leading_minors(gram: list[list[int]]) -> list[Fraction]:
    """The leading principal minors, by Gaussian elimination without pivoting."""
    work = [[Fraction(entry) for entry in row] for row in gram]
    minors, determinant = [], Fraction(1)
    for k in range(len(work)):
        pivot = work[k][k]
        assert pivot != 0, f"leading minor {k + 1} of {gram} vanishes"
        determinant *= pivot
        minors.append(determinant)
        for i in range(k + 1, len(work)):
            factor = work[i][k] / pivot
            work[i] = [a - factor * b for a, b in zip(work[i], work[k], strict=True)]
    return minors


def signature(gram: list[list[int]]) -> list[int]:
    """Sylvester: the signs of d_k / d_{k-1} over the leading principal minors."""
    minors = leading_minors(gram)
    ratios = [minors[0]] + [minors[k] / minors[k - 1] for k in range(1, len(minors))]
    return [sum(ratio > 0 for ratio in ratios), sum(ratio < 0 for ratio in ratios)]


def main() -> None:
    lines = read_source()
    genus = cite(
        lines,
        12933,
        12939,
        "two inequivalent indefinite ternary forms of determinant $^ { - 1 2 8 }$ of the same genus",
        r"f = \mathrm { d i a g } \left\{ - 1 , 6 4 , 2 \right\}",
        "contains two spinor genera and hence two classes",
        "The matrix (51b) is therefore a representative for the second class in this genus.",
    )
    form_a, form_b = displayed_arrays(lines[12936 - 1])

    for gram in (form_a, form_b):
        assert leading_minors(gram)[-1] == -128, f"{gram} does not have determinant -128"
        assert signature(gram) == [2, 1], f"{gram} does not have signature (2, 1)"
    one, third = Fraction(1), Fraction(1, 3)
    assert congruent(form_a, [[one, 0, 0], [one, one, 0], [0, 0, one]]) == [[-1, 0, 0], [0, 64, 0], [0, 0, 2]], (
        "(51a) does not have the diagonal basis e1, e1 + e2, e3 of diag{-1, 64, 2}"
    )
    assert congruent(form_a, [[3 * one, 0, 0], [0, third, 0], [0, 0, one]]) == form_b, "(51b) is not the Gram matrix of 3 e1, e2 / 3, e3"

    record = {
        "spinor_genus_pair_determinant_minus_128": {
            "id": "cs99_ch15_sec11_forms_51ab",
            "genus": "I_{2,1}(2 x 64)",
            "determinant": -128,
            "signature": [2, 1],
            "form_a": {"equation": "(51a)", "gram": form_a},
            "form_b": {"equation": "(51b)", "gram": form_b},
            "same_genus": True,
            "spinor_genera_in_genus": 2,
            "classes_in_genus": 2,
            "integrally_equivalent": False,
            "source": genus,
        }
    }
    TARGET.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}")


if __name__ == "__main__":
    main()
