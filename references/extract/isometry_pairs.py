"""Extract tests/fixtures/indefinite_isometry_pairs.json from vendored tests.

Pairs of indefinite Gram matrices with a recorded isometry verdict, from two
independent implementations' test suites:

- Indefinite.jl@374a5ebb TestLor/{U_I3,U_E8,U_2U_2I3}_mat{1,2}, which its driver
  test_gap.jl tests for equivalence (no witness recorded);
- Hecke.jl@e2ab5716 test/QuadForm/Quad/ZLatticeAutIso.jl, testset "isometry
  testing": U(2)+A2 against its transform by the recorded unimodular u, U against
  U(2), and a recorded isometric pair of 3x3 Grams.

Each verdict carries a certificate checked here, independently of the port:
a witness W with W * gram1 * W^T == gram2; or, for a pair without a witness,
equal genera with no spinor generators (one improper spinor genus), which for an
indefinite lattice of rank >= 3 is a single isometry class (Eichler); or, for a
non-isometric pair, different determinants.

Run from the repository root under Sage's Python:
    "$(dirname $(sage -c 'import sys; print(sys.executable)'))/python3" references/extract/isometry_pairs.py
"""

import json
import re
from pathlib import Path

# Importing from sage.all runs Sage's session startup.
from sage.all import AA, ZZ, Genus, matrix

INDEF = Path("references/vendor/Indefinite.jl@374a5ebb")
HECKE = Path("references/vendor/Hecke.jl@e2ab5716/test/QuadForm/Quad/ZLatticeAutIso.jl")
TARGET = Path("tests/fixtures/indefinite_isometry_pairs.json")


def read_indefinite_matrix(path: Path):
    tokens = [int(token) for token in path.read_text(encoding="utf-8").split()]
    rows, cols, entries = tokens[0], tokens[1], tokens[2:]
    assert len(entries) == rows * cols, path
    return matrix(ZZ, rows, cols, entries)


def julia_matrix(text: str):
    rows = [[int(entry) for entry in row.split()] for row in text.split(";")]
    return matrix(ZZ, rows)


def signature(gram) -> list[int]:
    roots = gram.charpoly().roots(AA, multiplicities=True)
    positive = sum(m for r, m in roots if r > 0)
    negative = sum(m for r, m in roots if r < 0)
    assert positive + negative == gram.nrows(), "degenerate"
    return [int(positive), int(negative)]


def as_rows(m) -> list[list[int]]:
    return [[int(entry) for entry in row] for row in m.rows()]


def genus_certificate(gram1, gram2) -> str:
    assert gram1.nrows() >= 3, "the Eichler certificate needs rank >= 3"
    assert Genus(gram1) == Genus(gram2), "the pair is not in one genus"
    assert Genus(gram1).spinor_generators(proper=False) == [], "the genus has more than one spinor genus"
    return "equal genera, no spinor generators (one improper spinor genus); indefinite of rank >= 3, so one isometry class (Eichler)"


def record(case_id, gram1, gram2, isometric, certificate, witness, source):
    sig = signature(gram1)
    assert sig == signature(gram2) and sig[0] > 0 and sig[1] > 0, f"{case_id}: not an indefinite pair of one signature"
    entry = {
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
    cases = []
    driver = (INDEF / "test_gap.jl").read_text(encoding="utf-8")
    for name in ("U_I3", "U_E8", "U_2U_2I3"):
        assert f'"TestLor/{name}_mat1"' in driver and f'"TestLor/{name}_mat2"' in driver
        gram1 = read_indefinite_matrix(INDEF / "TestLor" / f"{name}_mat1")
        gram2 = read_indefinite_matrix(INDEF / "TestLor" / f"{name}_mat2")
        source = {
            "kind": "independent_implementation_test",
            "files": [str(INDEF / "TestLor" / f"{name}_mat{i}") for i in (1, 2)],
            "driver": str(INDEF / "test_gap.jl"),
        }
        cases.append(record(f"indefinite_jl_{name}", gram1, gram2, True, genus_certificate(gram1, gram2), None, source))

    lines = HECKE.read_text(encoding="utf-8").splitlines()
    start = next(i for i, line in enumerate(lines) if '@testset "isometry testing"' in line)
    block = "\n".join(lines[start : start + 33])
    span = f"{start + 1}-{start + 33}"
    source = {"kind": "independent_implementation_test", "file": str(HECKE), "lines": span}
    u = julia_matrix(re.search(r"u = ZZ\[(.*?)\]", block).group(1))
    gram_l = julia_matrix(re.search(r"L = integer_lattice\(gram=ZZ\[(.*?)\]\)", block).group(1))
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
