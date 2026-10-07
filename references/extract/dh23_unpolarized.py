#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/unpolarized_enriques.json from Dutour Sikirić--Hulek section 3.

Source: arXiv:2302.01679, subsection "The Tits building" (vendored TeX, lines cited
per record). N = U + U(2) + E8(-2) with (e1, e2) and (e3, e4) the standard bases of
U and U(2). The isotropic line representatives are Z e1 and Z e3; the plane
representatives are Z e1 + Z e3 and Z (2 e1 + 2 e2 + w) + Z e3 with w of norm 4 in
E8, i.e. w^2 = -8 in E8(-2). The discriminant-image indices of the stabilizers are
1 and 527 (lines), 527 and 23715 (planes); the stable group has 528 lines and 24242
planes. The flag count 72199 is case 87 of Table table_subgroups2.

Run from the repository root: uv run references/extract/dh23_unpolarized.py
"""

import json
from fractions import Fraction
from pathlib import Path

SOURCE = "references/vendor/arxiv/2302.01679/Enriques_compu_rev.tex"
TARGET = Path("tests/fixtures/unpolarized_enriques.json")

# E8 Cartan matrix in the node order used throughout the fixtures: a chain
# 0-1-2-3-4-5-6 with node 7 attached to node 2.
E8_EDGES = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (2, 7))


def e8_cartan() -> list[list[int]]:
    cartan = [[2 if i == j else 0 for j in range(8)] for i in range(8)]
    for i, j in E8_EDGES:
        cartan[i][j] = cartan[j][i] = -1
    return cartan


def block_sum(*blocks: list[list[int]]) -> list[list[int]]:
    size = sum(len(block) for block in blocks)
    gram = [[0] * size for _ in range(size)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, entry in enumerate(row):
                gram[offset + i][offset + j] = entry
        offset += len(block)
    return gram


def determinant(matrix: list[list[int]]) -> int:
    rows = [[Fraction(entry) for entry in row] for row in matrix]
    size, sign, result = len(rows), 1, Fraction(1)
    for column in range(size):
        pivot = next(row for row in range(column, size) if rows[row][column] != 0)
        if pivot != column:
            rows[column], rows[pivot] = rows[pivot], rows[column]
            sign = -sign
        result *= rows[column][column]
        for row in range(column + 1, size):
            factor = rows[row][column] / rows[column][column]
            rows[row] = [a - factor * b for a, b in zip(rows[row], rows[column], strict=True)]
    value = sign * result
    assert value.denominator == 1
    return int(value)


def pairing(gram: list[list[int]], left: list[int], right: list[int]) -> int:
    return sum(left[i] * gram[i][j] * right[j] for i in range(len(left)) for j in range(len(right)))


def unit(index: int) -> list[int]:
    return [1 if position == index else 0 for position in range(12)]


def main() -> None:
    gram = block_sum(
        [[0, 1], [1, 0]],
        [[0, 2], [2, 0]],
        [[-2 * entry for entry in row] for row in e8_cartan()],
    )
    e1, e2, e3 = unit(0), unit(1), unit(2)
    # w = alpha_0 + alpha_2: two orthogonal E8 roots, so w has norm 4 in E8 and -8 in E8(-2).
    w = [0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0]
    assert pairing(gram, w, w) == -8
    v = [2 * a + 2 * b + c for a, b, c in zip(e1, e2, w, strict=True)]
    lines = {"I_1_1": [e1], "I_1_2": [e3]}
    planes = {"I_2_1": [e1, e3], "I_2_2": [v, e3]}
    for name, basis in {**lines, **planes}.items():
        assert all(pairing(gram, x, y) == 0 for x in basis for y in basis), f"{name} is not totally isotropic"
    det = determinant(gram)
    assert det == 1024, det

    def source(lines_: str) -> dict[str, str]:
        return {"kind": "published_text", "citation": "Dutour Sikirić--Hulek, arXiv:2302.01679", "file": SOURCE, "lines": lines_}

    record = {
        "id": "enriques_unpolarized_N",
        "lattice": {
            "name": "N",
            "construction": "U + U(2) + E8(-2)",
            "signature": [2, 10],
            "rank": 12,
            "det": det,
            "is_even": True,
            "gram": gram,
        },
        "component_preserving_group": {
            "internal_name": "O_component",
            "construction": "kernel(component_character)",
            "source_notation": "O^+(N)",
            "counts": {"line_orbits": 2, "plane_orbits": 2},
            "line_representatives": [
                {"id": "I_1_1", "basis": lines["I_1_1"], "description": "Z e1 (from U)", "discriminant_stabilizer_index": 1},
                {"id": "I_1_2", "basis": lines["I_1_2"], "description": "Z e3 (from U(2))", "discriminant_stabilizer_index": 527},
            ],
            "plane_representatives": [
                {"id": "I_2_1", "basis": planes["I_2_1"], "description": "Z e1 + Z e3", "discriminant_stabilizer_index": 527},
                {
                    "id": "I_2_2",
                    "basis": planes["I_2_2"],
                    "description": "Z (2e1 + 2e2 + w) + Z e3, w = alpha_0 + alpha_2 of norm 4 in E8",
                    "discriminant_stabilizer_index": 23715,
                },
            ],
            "source": source("1159-1166"),
        },
        "stable_component_preserving_group": {
            "internal_name": "O_plus_component",
            "construction": "intersection(kernel(discriminant_representation), kernel(component_character))",
            "source_notation": "\\widetilde O^+(N)",
            "counts": {"line_orbits": 528, "plane_orbits": 24242, "flag_orbits": 72199},
            "source": source("1168-1170; flag_orbits from table_subgroups2 case 87, line 814"),
        },
    }
    TARGET.write_text(json.dumps(record, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}: det {det}, all four representatives totally isotropic")


if __name__ == "__main__":
    main()
