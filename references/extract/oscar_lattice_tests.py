"""Extract tests/fixtures/oscar_lattice_oracles.json from OSCAR's vendored tests.

OSCAR is an independent implementation; its tests record outputs for indefinite
lattices that this port computes. Sources (references/vendor/Oscar.jl@d135b70b/test):

- NumberTheory/vinberg.jl: Vinberg's algorithm on four reflective lattices -- the
  number of simple roots (independent of the chamber), and for diag(1,-1,-3) the
  simple roots of the chamber the test fixes;
- Groups/spinor_norms.jl: for rank-3 lattices (table from_sage, active rows and the
  rows commented out to speed up OSCAR's CI), |O(q_L)| and the orders of the images
  in O(q_L) of O(L) (image_in_Oq) and of its subgroup of real spinor norm
  (image_in_Oq_signed); plus two cases whose image map is recorded bijective;
- Groups/isometry_group.jl: |O(U)| = 4, and a vector stabilizer of order 2.

The script asserts each recorded value on its source line, that every lattice is
indefinite and nondegenerate, and that explicit roots give integral reflections.

Run from the repository root under Sage's Python:
    "$(dirname $(sage -c 'import sys; print(sys.executable)'))/python3" references/extract/oscar_lattice_tests.py
"""

import json
import re
from pathlib import Path

# Importing from sage.all runs Sage's session startup.
from sage.all import AA, ZZ, diagonal_matrix, matrix

TESTS = Path("references/vendor/Oscar.jl@d135b70b/test")
TARGET = Path("tests/fixtures/oscar_lattice_oracles.json")


def lines_of(relative: str) -> list[str]:
    return (TESTS / relative).read_text(encoding="utf-8").splitlines()


def at(relative: str, line: int, needle: str) -> dict[str, object]:
    text = lines_of(relative)[line - 1]
    assert needle in text, f"{TESTS / relative}:{line} does not contain {needle!r}: {text.strip()!r}"
    return {"kind": "independent_implementation_test", "file": str(TESTS / relative), "line": line}


def find(relative: str, needle: str, start: int = 1) -> int:
    return next(i for i, line in enumerate(lines_of(relative), start=1) if i >= start and needle in line)


def indefinite(gram) -> list[int]:
    roots = gram.charpoly().roots(AA, multiplicities=True)
    positive = sum(m for r, m in roots if r > 0)
    negative = sum(m for r, m in roots if r < 0)
    assert positive + negative == gram.nrows() and positive > 0 and negative > 0, f"not indefinite: {gram}"
    return [int(positive), int(negative)]


def rows(m) -> list[list[int]]:
    return [[int(entry) for entry in row] for row in m.rows()]


def integral_reflection(gram, root) -> None:
    r = matrix(ZZ, [root])
    square = (r * gram * r.transpose())[0, 0]
    assert square != 0 and all((2 * entry) % square == 0 for entry in (gram * r.transpose()).list()), f"{root} gives no integral reflection"


def vinberg() -> list[dict[str, object]]:
    f = "NumberTheory/vinberg.jl"
    cases = []
    i1_10 = diagonal_matrix(ZZ, [1] + [-1] * 10)
    line = find(f, "@test length(roots) == 12")
    last = find(f, "@test roots[12] == ZZ[3 1 1 1 1 1 1 1 1 1 1]")
    integral_reflection(i1_10, [3] + [1] * 10)
    cases.append(
        {
            "id": "vinberg_I_1_10",
            "gram": rows(i1_10),
            "signature": indefinite(i1_10),
            "num_simple_roots": 12,
            "a_simple_root_up_to_sign": [3] + [1] * 10,
            "source": [at(f, line, "== 12"), at(f, last, "ZZ[3 1 1 1 1 1 1 1 1 1 1]")],
        }
    )
    grosek = diagonal_matrix(ZZ, [1, -1, -13])
    line = find(f, "@test l == 8")
    cusps = find(f, "@test cusp_counter == 2")
    cases.append(
        {
            "id": "vinberg_diag_1_-1_-13",
            "gram": rows(grosek),
            "signature": indefinite(grosek),
            "num_simple_roots": 8,
            "cusps": 2,
            "source": [at(f, line, "== 8"), at(f, cusps, "== 2")],
        }
    )
    scharlau = matrix(ZZ, [[0, -1, 0, 0], [-1, 0, 0, 0], [0, 0, -2, -1], [0, 0, -1, -2]])
    line = find(f, "@test l == 4")
    edges = find(f, "@test edge_counter == 6")
    at(f, find(f, "Q = matrix(ZZ, 4, 4, [0 -1 0 0; -1 0 0 0; 0 0 -2 -1; 0 0 -1 -2])"), "0 0 -2 -1")
    cases.append(
        {
            "id": "vinberg_scharlau_U_A2",
            "gram": rows(scharlau),
            "signature": indefinite(scharlau),
            "num_simple_roots": 4,
            "coxeter_diagram_edges": 6,
            "source": [at(f, line, "== 4"), at(f, edges, "== 6")],
        }
    )
    small = diagonal_matrix(ZZ, [1, -1, -3])
    roots, sources = [], []
    for k, root in enumerate(([0, -1, 0], [0, 0, -1], [1, 0, 1], [3, 3, 1]), start=1):
        needle = f"@test roots[{k}] == ZZ[{' '.join(map(str, root))}]"
        sources.append(at(f, find(f, needle), needle))
        integral_reflection(small, root)
        roots.append(root)
    cases.append(
        {
            "id": "vinberg_diag_1_-1_-3",
            "gram": rows(small),
            "signature": indefinite(small),
            "num_simple_roots": 4,
            "simple_roots_of_test_chamber": roots,
            "chamber": "v0 = (1,0,0), direction (0,9,6)",
            "source": sources,
        }
    )
    return cases


def discriminant_images() -> list[dict[str, object]]:
    f = "Groups/spinor_norms.jl"
    lines = lines_of(f)
    start = find(f, "from_sage = [")
    stop = find(f, "for (g,ks,k,n) in from_sage", start)
    entry = re.compile(r"\(\[([-\d, ]+)\],\s*(\d+),\s*(\d+),\s*(\d+)\)")
    cases, commented = [], False
    for number in range(start, stop):
        text = lines[number - 1]
        if "#=" in text:
            commented = True
        for match in entry.finditer(text):
            gram = matrix(ZZ, 3, 3, [int(x) for x in match.group(1).split(",")])
            signed, full, oq = (int(match.group(i)) for i in (2, 3, 4))
            cases.append(
                {
                    "id": f"from_sage_{len(cases)}",
                    "gram": rows(gram),
                    "signature": indefinite(gram),
                    "image_in_Oq_order": full,
                    "image_in_Oq_signed_order": signed,
                    "O_qL_order": oq,
                    "run_in_oscar_ci": not commented,
                    "source": {"kind": "independent_implementation_test", "file": str(TESTS / f), "line": number},
                }
            )
        if "=#" in text:
            commented = False
    assert len(cases) == 28 and sum(c["run_in_oscar_ci"] for c in cases) == 10, (len(cases), sum(c["run_in_oscar_ci"] for c in cases))
    a7 = find(f, "L = root_lattice(:A, 7)", stop)
    bij1 = find(f, "@test is_bijective(GLinOq)", a7)
    a7_gram = matrix(ZZ, 10, 10, 0)
    for i in range(7):
        a7_gram[i, i] = 2
        if i < 6:
            a7_gram[i, i + 1] = a7_gram[i + 1, i] = -1
    a7_gram[7, 7], a7_gram[8, 8], a7_gram[9, 9] = 1, 1, -1
    cases.append(
        {
            "id": "A7_plus_diag_1_1_-1",
            "gram": rows(a7_gram),
            "signature": indefinite(a7_gram),
            "image_in_Oq_is_all_of_O_qL": True,
            "source": [at(f, a7, "root_lattice(:A, 7)"), at(f, a7 + 1, "QQ[1 0 0; 0 1 0; 0 0 -1]"), at(f, bij1, "is_bijective")],
        }
    )
    u3 = find(f, "L = integer_lattice(; gram = QQ[0 1 0;1 0 0; 0 0 -3])", bij1)
    u3_gram = matrix(ZZ, [[0, 1, 0], [1, 0, 0], [0, 0, -3]])
    cases.append(
        {
            "id": "U_plus_-3",
            "gram": rows(u3_gram),
            "signature": indefinite(u3_gram),
            "image_in_Oq_signed_order": 2,
            "image_in_Oq_signed_is_all_of_O_qL": True,
            "source": [at(f, u3, "0 0 -3"), at(f, u3 + 2, "@test order(GL) == 2"), at(f, u3 + 3, "is_bijective")],
        }
    )
    return cases


def isometry_group() -> list[dict[str, object]]:
    f = "Groups/isometry_group.jl"
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    h = find(f, "H = hyperbolic_plane_lattice()")
    cases = [
        {
            "id": "O_U_order",
            "gram": rows(hyperbolic),
            "signature": indefinite(hyperbolic),
            "orthogonal_group_order": 4,
            "source": [at(f, h, "hyperbolic_plane_lattice()"), at(f, h + 2, "@test order(G) == 4")],
        }
    ]
    s = find(f, "L = integer_lattice(gram=ZZ[4 1 1; 1 -2 0; 1 0 -2;])")
    gram = matrix(ZZ, [[4, 1, 1], [1, -2, 0], [1, 0, -2]])
    cases.append(
        {
            "id": "stabilizer_of_100",
            "gram": rows(gram),
            "signature": indefinite(gram),
            "vector": [1, 0, 0],
            "vector_stabilizer_order": 2,
            "source": [at(f, s, "[4 1 1; 1 -2 0; 1 0 -2;]"), at(f, s + 1, "QQ[1 0 0;]"), at(f, s + 3, "@test order(G)==2")],
        }
    )
    return cases


def main() -> None:
    data = {"vinberg": vinberg(), "discriminant_images": discriminant_images(), "isometry_groups": isometry_group()}
    TARGET.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}: " + ", ".join(f"{k} {len(v)}" for k, v in data.items()))


if __name__ == "__main__":
    main()
