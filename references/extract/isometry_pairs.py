"""Extract tests/fixtures/indefinite_isometry_pairs.json from asserted upstream tests.

Pairs of indefinite Gram matrices with a *recorded* isometry verdict from
Hecke.jl@e2ab5716 test/QuadForm/Quad/ZLatticeAutIso.jl, testset "isometry
  testing": U(2)+A2 against its transform by the recorded unimodular u, U against
  U(2), and a recorded isometric pair of 3x3 Grams.

Each verdict is explicitly asserted by Hecke's upstream test and checked
here by a witness W with W * gram1 * W^T == gram2, or by differing determinants.

Indefinite.jl's test_gap.jl lists three TestLor matrix pairs, but
CaseEquiv1/2/3 are all set to false (lines 34-52), and the code does not
assert a verdict for them. They must NOT be assigned 'isometric: true'
via a self-authored genus argument and presented as pinned upstream oracles.

Run from the repository root under Sage's Python:
    "$(dirname $(sage -c 'import sys; print(sys.executable)'))/python3" references/extract/isometry_pairs.py
"""

import json
import re
from pathlib import Path
from typing import NotRequired, TypedDict

# Importing from sage.all runs Sage's session startup.
from sage.all import AA, ZZ, matrix
from sage.matrix.matrix_integer_dense import Matrix_integer_dense

HECKE = Path("references/vendor/Hecke.jl@e2ab5716/test/QuadForm/Quad/ZLatticeAutIso.jl")
TARGET = Path("tests/fixtures/indefinite_isometry_pairs.json")


class IsometryPair(TypedDict):
    id: str
    gram1: list[list[int]]
    gram2: list[list[int]]
    signature: list[int]
    isometric: bool
    certificate: str
    source: dict[str, str | list[str]]
    witness: NotRequired[list[list[int]]]
    witness_convention: NotRequired[str]


def julia_matrix(text: str) -> Matrix_integer_dense:
    rows = [[int(entry) for entry in row.split()] for row in text.split(";")]
    return matrix(ZZ, rows)


def signature(gram: Matrix_integer_dense) -> list[int]:
    roots = gram.charpoly().roots(AA, multiplicities=True)
    positive = sum(m for r, m in roots if r > 0)
    negative = sum(m for r, m in roots if r < 0)
    assert positive + negative == gram.nrows(), "degenerate"
    return [int(positive), int(negative)]


def as_rows(m: Matrix_integer_dense) -> list[list[int]]:
    return [[int(entry) for entry in row] for row in m.rows()]


def record(
    case_id: str,
    gram1: Matrix_integer_dense,
    gram2: Matrix_integer_dense,
    isometric: bool,
    certificate: str,
    witness: Matrix_integer_dense | None,
    source: dict[str, str | list[str]],
) -> IsometryPair:
    sig = signature(gram1)
    assert sig == signature(gram2) and sig[0] > 0 and sig[1] > 0, f"{case_id}: not an indefinite pair of one signature"
    entry: IsometryPair = {
        "id": case_id,
        "gram1": as_rows(gram1),
        "gram2": as_rows(gram2),
        "signature": sig,
        "isometric": isometric,
        "certificate": certificate,
        "source": source,
    }
    if witness is not None:
        assert witness.det().abs() == 1 and witness * gram1 * witness.transpose() == gram2, f"{case_id}: witness fails"
        entry["witness"] = as_rows(witness)
        entry["witness_convention"] = "witness * gram1 * witness^T == gram2"
    return entry


def main() -> None:
    cases: list[IsometryPair] = []
    lines = HECKE.read_text(encoding="utf-8").splitlines()
    start = next(i for i, line in enumerate(lines) if '@testset "isometry testing"' in line)
    block = "\n".join(lines[start : start + 33])
    span = f"{start + 1}-{start + 33}"
    source = {"kind": "independent_implementation_test", "file": str(HECKE), "lines": span}
    u_match = re.search(r"u = ZZ\[(.*?)\]", block)
    assert u_match, f"{HECKE}:{span}: no witness u = ZZ[...]"
    u = julia_matrix(u_match.group(1))
    gram_l_match = re.search(r"L = integer_lattice\(gram=ZZ\[(.*?)\]\)", block)
    assert gram_l_match, f"{HECKE}:{span}: no L = integer_lattice(gram=ZZ[...])"
    gram_l = julia_matrix(gram_l_match.group(1))
    assert "@test Hecke.is_isometric(L, M)" in block
    cases.append(record("hecke_U2_A2_by_u", gram_l, u * gram_l * u.transpose(), True, "recorded witness u", u, source))
    assert "@test !is_isometric(L1, L2)" in block and "L2 = hyperbolic_plane_lattice(2)" in block
    hyperbolic, hyperbolic_2 = matrix(ZZ, [[0, 1], [1, 0]]), matrix(ZZ, [[0, 2], [2, 0]])
    assert hyperbolic.det() != hyperbolic_2.det()
    cases.append(record("hecke_U_vs_U2", hyperbolic, hyperbolic_2, False, "determinants -1 and -4 differ", None, source))
    grams = re.findall(r"G = matrix\(QQ, 3, 3 ,\[([^\]]*)\]\);", block)
    assert len(grams) == 2 and "@test is_isometric(L1,L2)" in block
    g1, g2 = (matrix(ZZ, 3, 3, [int(x) for x in g.split(",")]) for g in grams)
    swap = matrix(ZZ, [[1, 0, 0], [0, 0, 1], [0, 1, 0]])
    cases.append(record("hecke_3x3_swapped_basis", g1, g2, True, "the second Gram is the first with basis vectors 2 and 3 swapped", swap, source))

    TARGET.write_text(json.dumps(cases, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {len(cases)} pairs to {TARGET}")


if __name__ == "__main__":
    main()
