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

The same file records facts about definite lattices, kept in their own keys because
definite isometry groups are the research preamble's concern; this port reaches them
only through the preamble:

- definite_orthogonal_group_orders: |O(L)| for every lattice whose order the file
  records, including all 24 entries of the list LL/orders (OSCAR's CI runs entries
  5-10; run_in_oscar_ci says which);
- definite_isometry_tests: one recorded non-isometric pair and one isometric pair;
- root_lattice_groups: for A_i, i = 2..5, the isomorphism types of the stable,
  special and special stable orthogonal groups, with their orders;
- hecke_reduced_automorphism_group_orders: Hecke's reduced_automorphism_group_order
  for A2 and E8, recorded under Hecke's name.

The script asserts each recorded value on its source line, that every lattice in the
indefinite keys is indefinite and nondegenerate, that every lattice in the definite
keys is definite, and that explicit roots give integral reflections.

Run from the repository root under Sage's Python:
    "$(dirname $(sage -c 'import sys; print(sys.executable)'))/python3" references/extract/oscar_lattice_tests.py
"""

import json
import re
from pathlib import Path
from typing import NotRequired, TypedDict

# Importing from sage.all runs Sage's session startup.
from sage.all import AA, ZZ, diagonal_matrix, factorial, matrix
from sage.matrix.matrix2 import Matrix
from sage.rings.integer import Integer

TESTS = Path("references/vendor/Oscar.jl@d135b70b/test")
TARGET = Path("tests/fixtures/oscar_lattice_oracles.json")

type Source = dict[str, str | int]


class VinbergCase(TypedDict):
    id: str
    gram: list[list[int]]
    signature: list[int]
    num_simple_roots: int
    a_simple_root_up_to_sign: NotRequired[list[int]]
    cusps: NotRequired[int]
    coxeter_diagram_edges: NotRequired[int]
    simple_roots_of_test_chamber: NotRequired[list[list[int]]]
    chamber: NotRequired[str]
    source: list[Source]


class DiscriminantCase(TypedDict):
    id: str
    gram: list[list[int]]
    signature: list[int]
    image_in_Oq_order: NotRequired[int]
    image_in_Oq_signed_order: NotRequired[int]
    O_qL_order: NotRequired[int]
    run_in_oscar_ci: NotRequired[bool]
    image_in_Oq_is_all_of_O_qL: NotRequired[bool]
    image_in_Oq_signed_is_all_of_O_qL: NotRequired[bool]
    source: Source | list[Source]


class IsometryGroupCase(TypedDict):
    id: str
    gram: list[list[int]]
    signature: list[int]
    orthogonal_group_order: NotRequired[int]
    vector: NotRequired[list[int]]
    vector_stabilizer_order: NotRequired[int]
    source: list[Source]


class OrderCase(TypedDict):
    id: str
    gram: list[list[int]]
    signature: list[int]
    orthogonal_group_order: int
    run_in_oscar_ci: NotRequired[bool]
    source: list[Source]


class IsometryTest(TypedDict):
    id: str
    gram_1: list[list[int]]
    gram_2: list[list[int]]
    signature_1: list[int]
    signature_2: list[int]
    isometric: bool
    source: list[Source]


class GroupType(TypedDict):
    isomorphism_type: str
    order: int


class RootLatticeGroups(TypedDict):
    id: str
    gram: list[list[int]]
    signature: list[int]
    stable_orthogonal_group: GroupType
    special_orthogonal_group: GroupType
    special_stable_orthogonal_group: GroupType
    source: list[Source]


class ReducedOrder(TypedDict):
    lattice: str
    hecke_reduced_automorphism_group_order: int
    source: Source


type OracleRecords = (
    list[VinbergCase] | list[DiscriminantCase] | list[IsometryGroupCase] | list[OrderCase] | list[IsometryTest] | list[RootLatticeGroups] | list[ReducedOrder]
)


def lines_of(relative: str) -> list[str]:
    return (TESTS / relative).read_text(encoding="utf-8").splitlines()


def at(relative: str, line: int, needle: str) -> Source:
    text = lines_of(relative)[line - 1]
    assert needle in text, f"{TESTS / relative}:{line} does not contain {needle!r}: {text.strip()!r}"
    return {"kind": "independent_implementation_test", "file": str(TESTS / relative), "line": line}


def find(relative: str, needle: str, start: int = 1) -> int:
    return next(i for i, line in enumerate(lines_of(relative), start=1) if i >= start and needle in line)


def indefinite(gram: Matrix[Integer]) -> list[int]:
    roots = gram.charpoly().roots(AA, multiplicities=True)
    positive = sum(m for r, m in roots if r > 0)
    negative = sum(m for r, m in roots if r < 0)
    assert positive + negative == gram.nrows() and positive > 0 and negative > 0, f"not indefinite: {gram}"
    return [int(positive), int(negative)]


def rows(m: Matrix[Integer]) -> list[list[int]]:
    return [[int(entry) for entry in row] for row in m.rows()]


def integral_reflection(gram: Matrix[Integer], root: list[int]) -> None:
    r = matrix(ZZ, [root])
    square = (r * gram * r.transpose())[0, 0]
    assert square != 0 and all((2 * entry) % square == 0 for entry in (gram * r.transpose()).list()), f"{root} gives no integral reflection"


def vinberg() -> list[VinbergCase]:
    f = "NumberTheory/vinberg.jl"
    cases: list[VinbergCase] = []
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
    roots: list[list[int]] = []
    sources: list[Source] = []
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


def discriminant_images() -> list[DiscriminantCase]:
    f = "Groups/spinor_norms.jl"
    lines = lines_of(f)
    start = find(f, "from_sage = [")
    stop = find(f, "for (g,ks,k,n) in from_sage", start)
    entry = re.compile(r"\(\[([-\d, ]+)\],\s*(\d+),\s*(\d+),\s*(\d+)\)")
    cases: list[DiscriminantCase] = []
    commented = False
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


def isometry_group() -> list[IsometryGroupCase]:
    f = "Groups/isometry_group.jl"
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    h = find(f, "H = hyperbolic_plane_lattice()")
    cases: list[IsometryGroupCase] = [
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


GRAM_LITERAL = re.compile(r"(?:ZZ|QQ)\[([^\]]*)\]")


def gram_on(relative: str, line: int) -> Matrix[Integer]:
    """The single Gram literal on a source line, as an integer matrix."""
    literals = GRAM_LITERAL.findall(lines_of(relative)[line - 1])
    assert len(literals) == 1, f"{TESTS / relative}:{line}: expected one Gram literal, found {len(literals)}"
    rows_ = [[int(entry) for entry in row.split()] for row in literals[0].split(";") if row.strip()]
    return matrix(ZZ, rows_)


def definite(gram: Matrix[Integer]) -> list[int]:
    assert gram.is_symmetric() and gram.det() != 0, f"degenerate or asymmetric: {gram}"
    if gram.is_positive_definite():
        return [gram.nrows(), 0]
    assert (-gram).is_positive_definite(), f"not definite: {gram}"
    return [0, gram.nrows()]


def evaluate_order(expression: str) -> int:
    """An order as written in the test: a product of integers and powers, e.g. 696729600^2*2*2."""
    value = 1
    for factor in expression.split("*"):
        base, _, exponent = factor.strip().partition("^")
        value *= int(base) ** int(exponent or 1)
    return value


def root_lattice_a(n: int) -> Matrix[Integer]:
    """OSCAR's root_lattice(:A, n): the positive definite Cartan Gram matrix."""
    return matrix(ZZ, n, n, lambda i, j: 2 if i == j else (-1 if abs(i - j) == 1 else 0))


def definite_orthogonal_groups() -> dict[str, OracleRecords]:
    f = "Groups/isometry_group.jl"
    order_test = re.compile(r"@test order\((.*)\)\s*==\s*([\d^*]+)\s*$")

    def recorded_order(line: int) -> int:
        match = order_test.search(lines_of(f)[line - 1])
        assert match, f"{TESTS / f}:{line}: no recorded order"
        return evaluate_order(match.group(2))

    def order_case(case_id: str, gram: Matrix[Integer], gram_source: list[Source], test_lines: list[int]) -> OrderCase:
        orders = {recorded_order(line) for line in test_lines}
        assert len(orders) == 1, f"{case_id}: the tests record different orders {orders}"
        return {
            "id": case_id,
            "gram": rows(gram),
            "signature": definite(gram),
            "orthogonal_group_order": orders.pop(),
            "source": gram_source + [at(f, line, "@test order(") for line in test_lines],
        }

    a2 = root_lattice_a(2)
    n_gram = a2.block_sum(4 * a2)
    cases: list[OrderCase] = [
        order_case(
            "A2_plus_A2(4)",
            n_gram,
            [at(f, 2, "root_lattice(:A, 2)"), at(f, 3, "rescale(N1, 4)"), at(f, 4, "direct_sum(N1,N2)")],
            [5],
        )
    ]
    assert gram_on(f, 7) == gram_on(f, 12)
    cases.append(order_case("isometry_group_line_7", gram_on(f, 7), [at(f, 7, "integer_lattice")], [9, 10, 11, 13]))
    cases.append(order_case("isometry_group_line_15", gram_on(f, 15), [at(f, 15, "gram = ZZ["), at(f, 16, "integer_lattice(; gram)")], [17]))
    for gram_line, test_line in ((48, 53), (59, 64), (72, 73), (77, 78), (81, 82), (84, 85), (88, 89)):
        cases.append(order_case(f"isometry_group_line_{gram_line}", gram_on(f, gram_line), [at(f, gram_line, "integer_lattice(gram=")], [test_line]))
    orders_line = find(f, "orders = ZZRingElem[")
    orders_match = re.search(r"ZZRingElem\[([\d, ]+)\]", lines_of(f)[orders_line - 1])
    assert orders_match, f"{TESTS / f}:{orders_line}: no ZZRingElem[...] order list"
    orders = [int(x) for x in orders_match.group(1).split(",")]
    ll_lines = list(range(94, 118))
    assert len(orders) == len(ll_lines) == 24
    loop = find(f, "for (L,ord) in zip(LL[5:10], orders[5:10])")
    at(f, loop + 1, "@test order(orthogonal_group(L))==ord")
    for position, (line, order) in enumerate(zip(ll_lines, orders, strict=True), start=1):
        gram = gram_on(f, line)
        cases.append(
            {
                "id": f"isometry_group_LL_{position}",
                "gram": rows(gram),
                "signature": definite(gram),
                "orthogonal_group_order": order,
                "run_in_oscar_ci": 5 <= position <= 10,
                "source": [at(f, line, "integer_lattice(gram="), at(f, orders_line, str(order)), at(f, loop, "zip(LL[5:10], orders[5:10])")],
            }
        )
    for case in cases:
        case.setdefault("run_in_oscar_ci", True)
    # Cross-check within the file: L2 (line 59, order written 696729600^2*2*2) is LL[1] (order from the list).
    by_id = {case["id"]: case for case in cases}
    assert by_id["isometry_group_line_59"]["gram"] == by_id["isometry_group_LL_1"]["gram"]
    assert by_id["isometry_group_line_59"]["orthogonal_group_order"] == by_id["isometry_group_LL_1"]["orthogonal_group_order"]

    def pair(case_id: str, first: int, second: int, test: int, needle: str, isometric: bool) -> IsometryTest:
        left, right = gram_on(f, first), gram_on(f, second)
        return {
            "id": case_id,
            "gram_1": rows(left),
            "gram_2": rows(right),
            "signature_1": definite(left),
            "signature_2": definite(right),
            "isometric": isometric,
            "source": [at(f, first, "integer_lattice(gram="), at(f, second, "integer_lattice(gram="), at(f, test, needle)],
        }

    isometry_tests = [
        pair("isometry_group_L1_L2_line_69", 48, 59, 69, "@test !is_isometric_with_isometry(L1,L2)[1]", False),
        pair("isometry_group_L1_L2_line_130", 127, 129, 130, "@test is_isometric(L1, L2)", True),
    ]

    at(f, 29, "for i in 2:5")
    at(f, 31, "S = symmetric_group(i+1)")
    at(f, 32, "A = alternating_group(i+1)")
    at(f, 33, "T = isodd(i) ? S : direct_product(A, cyclic_group(2))")
    groups: list[RootLatticeGroups] = []
    for i in range(2, 6):
        n = factorial(i + 1)
        groups.append(
            {
                "id": f"A{i}_orthogonal_subgroups",
                "gram": rows(root_lattice_a(i)),
                "signature": definite(root_lattice_a(i)),
                "stable_orthogonal_group": {"isomorphism_type": f"S_{i + 1}", "order": int(n)},
                "special_orthogonal_group": {
                    "isomorphism_type": f"S_{i + 1}" if i % 2 else f"A_{i + 1} x C_2",
                    "order": int(n),
                },
                "special_stable_orthogonal_group": {"isomorphism_type": f"A_{i + 1}", "order": int(n // 2)},
                "source": [at(f, 30, "root_lattice(:A, i)"), at(f, 35, "is_isomorphic(O_st, S)"), at(f, 38, "is_isomorphic(O_sp, T)"), at(f, 41, "is_isomorphic(O_spst, A)")],
            }
        )

    e8 = find(f, "reduced_automorphism_group_order(root_lattice(:E, 8)) == 1")
    reduced: list[ReducedOrder] = [
        {"lattice": "A2 (OSCAR root_lattice(:A, 2))", "hecke_reduced_automorphism_group_order": 2, "source": at(f, e8 - 1, "root_lattice(:A, 2)) == 2")},
        {"lattice": "E8 (OSCAR root_lattice(:E, 8))", "hecke_reduced_automorphism_group_order": 1, "source": at(f, e8, "root_lattice(:E, 8)) == 1")},
    ]
    return {
        "definite_orthogonal_group_orders": cases,
        "definite_isometry_tests": isometry_tests,
        "root_lattice_groups": groups,
        "hecke_reduced_automorphism_group_orders": reduced,
    }


def main() -> None:
    data: dict[str, OracleRecords] = {"vinberg": vinberg(), "discriminant_images": discriminant_images(), "isometry_groups": isometry_group(), **definite_orthogonal_groups()}
    TARGET.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}: " + ", ".join(f"{k} {len(v)}" for k, v in data.items()))


if __name__ == "__main__":
    main()
