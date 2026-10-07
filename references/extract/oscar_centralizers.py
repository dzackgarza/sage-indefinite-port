"""Extract tests/fixtures/isometry_centralizers.json from OSCAR's vendored tests.

Sources: references/vendor/Oscar.jl@d135b70b/test/NumberTheory/QuadFormAndIsom.
OSCAR is an independent implementation, so its recorded outputs are external
oracles for this port's centralizer computations. Each case is read mechanically
from the cited line range: the ambient basis matrix B, ambient Gram matrix G and
ambient isometry f (OSCAR acts on row vectors), and the recorded value from the
@test line in that range. The script computes the lattice's Gram matrix B G B^T
and the isometry in the lattice basis (the solution f_L of f_L B = B f), and
asserts that f_L is integral and preserves the Gram matrix. Cases are split by the
computed signature: indefinite ones into centralizer_cases, definite ones (the
research preamble's concern) into definite_centralizer_cases; the split is asserted.

enumeration.jl also records how many isomorphism classes of lattices with isometry
(L, f) OSCAR's enumeration returns. Those counts are kept with the OSCAR call verbatim,
since its filter arguments (char_poly, min_poly, pos_sigs, neg_sigs, fix_root, the
hermitian-type arguments) are OSCAR's parameters, together with every further @test
recorded on the result:
- lattice_class_counts: indefinite lattices given explicitly (U + E8, 5U, 4U, 2U, and
  the rank-6 lattice of the "Fix type condition" testset);
- hermitian_genus_class_counts: hermitian-type classes specified by signature pairs and
  determinant, with no explicit lattice;
- definite_class_counts: E6 and E8, the research preamble's concern.
The A4 loop at lines 15-18 compares OSCAR's enumeration with OSCAR's own conjugacy
classes and records no count, so it is not transcribed.

Run from the repository root under Sage's Python:
    "$(dirname $(sage -c 'import sys; print(sys.executable)'))/python3" references/extract/oscar_centralizers.py
"""

import json
import re
from fractions import Fraction
from pathlib import Path

# Importing from sage.all runs Sage's session startup.
from sage.all import AA, QQ, ZZ, CartanMatrix, block_diagonal_matrix, identity_matrix, matrix

TESTS = Path("references/vendor/Oscar.jl@d135b70b/test/NumberTheory/QuadFormAndIsom")
TARGET = Path("tests/fixtures/isometry_centralizers.json")
LITERAL = re.compile(r"(\w+)\s*=\s*matrix\(QQ,\s*(\d+),\s*(\d+)\s*,\s*\[(.*?)\]\)", re.DOTALL)

# (file, first line, last line, the @test pattern carrying the recorded value, record key)
CENTRALIZER_CASES = (
    ("lattices_with_isometry.jl", 80, 87, r"@test order\(GLf\) == (\d+)", "centralizer_image_order"),
    ("lattices_with_isometry.jl", 110, 116, r"@test order\(GL\) == (\d+)", "centralizer_image_order"),
    ("lattices_with_isometry.jl", 118, 124, r"@test order\(GL\) == (\d+)", "centralizer_image_order"),
    ("lattices_with_isometry.jl", 133, 139, r"@test (is_bijective)\(image_centralizer_in_Oq\(Lf\)\[2\]\)", "centralizer_image_is_all_of_O_qL"),
    ("enumeration.jl", 88, 94, r"@test order\(GLf\) == (\d+)", "centralizer_image_order"),
    ("enumeration.jl", 96, 102, r"@test order\(GLf\) == (\d+)", "centralizer_image_order"),
)
INVOLUTION_CLASSES = ("enumeration.jl", 52, 58)


def entries(body: str) -> list[Fraction]:
    tokens = [token for token in re.split(r"[\s,;]+", body.strip()) if token]
    return [Fraction(*map(int, token.split("//"))) if "//" in token else Fraction(int(token)) for token in tokens]


def literals(path: Path, first: int, last: int) -> tuple[dict[str, object], str]:
    text = "\n".join(path.read_text(encoding="utf-8").splitlines()[first - 1 : last])
    found = {}
    for name, rows, cols, body in LITERAL.findall(text):
        values = entries(body)
        assert len(values) == int(rows) * int(cols), f"{path}:{first}-{last}: {name} has {len(values)} entries"
        found[name] = matrix(QQ, int(rows), int(cols), [QQ(v.numerator) / v.denominator for v in values])
    return found, text


def signature(gram) -> list[int]:
    roots = gram.charpoly().roots(AA, multiplicities=True)
    positive = sum(m for r, m in roots if r > 0)
    negative = sum(m for r, m in roots if r < 0)
    assert positive + negative == gram.nrows(), "degenerate lattice"
    return [int(positive), int(negative)]


def integral(m) -> list[list[int]]:
    assert all(entry in ZZ for entry in m.list()), f"expected an integral matrix, got {m}"
    return [[int(entry) for entry in row] for row in m.rows()]


def source(relative: str, first: int, last: int) -> dict[str, str]:
    return {"kind": "independent_implementation_test", "file": str(TESTS / relative), "lines": f"{first}-{last}"}


def at(relative: str, line: int, needle: str) -> str:
    """The stripped text of a source line, asserted to contain needle."""
    text = (TESTS / relative).read_text(encoding="utf-8").splitlines()[line - 1].strip()
    assert needle in text, f"{TESTS / relative}:{line} does not contain {needle!r}: {text!r}"
    return text


def recorded_count(relative: str, line: int) -> int:
    match = re.fullmatch(r"@test length\((?:r|reps|r_global|r_local|.*\))\)\s*==\s*(\d+)", at(relative, line, "@test length("))
    assert match, f"{TESTS / relative}:{line}: no recorded count"
    return int(match.group(1))


def class_count(relative: str, case_id: str, call_line: int, count_line: int, property_lines: tuple[int, ...] = ()) -> dict[str, object]:
    return {
        "id": case_id,
        "oscar_call": at(relative, call_line, "("),
        "class_count": recorded_count(relative, count_line),
        "recorded_properties": [at(relative, line, "@test") for line in property_lines],
        "source": {"kind": "independent_implementation_test", "file": str(TESTS / relative), "lines": f"{call_line}-{max((count_line, *property_lines))}"},
    }


def class_counts() -> dict[str, list[dict[str, object]]]:
    e = "enumeration.jl"
    hyperbolic = matrix(ZZ, [[0, 1], [1, 0]])
    e8 = matrix(ZZ, CartanMatrix(["E", 8]))
    e6 = matrix(ZZ, CartanMatrix(["E", 6]))
    at(e, 39, "U = hyperbolic_plane_lattice()")
    at(e, 40, "E8 = root_lattice(:E, 8)")
    at(e, 41, "L = direct_sum(U, E8)[1]")
    at(e, 60, "M = direct_sum(U, U, U, U, U)[1]")
    at(e, 77, "L = direct_sum(U, U, U, U)[1]")
    at(e, 107, "L, _ = direct_sum(U, U)")
    found, _ = literals(TESTS / e, 119, 121)
    fix_type = found["B"] * found["G"] * found["B"].transpose()
    explicit = (
        (class_count(e, "U_plus_E8_order_30", 42, 43, (44,)), block_diagonal_matrix(hyperbolic, e8), 30),
        (class_count(e, "5U_order_4", 61, 62, (63,)), block_diagonal_matrix([hyperbolic] * 5), 4),
        (class_count(e, "4U_hermitian_order_5", 78, 78), block_diagonal_matrix([hyperbolic] * 4), 5),
        (class_count(e, "4U_hermitian_order_5_fix_5", 79, 79), block_diagonal_matrix([hyperbolic] * 4), 5),
        (class_count(e, "2U_hermitian_order_4_fix_4", 108, 109), block_diagonal_matrix([hyperbolic] * 2), 4),
        (class_count(e, "fix_type_condition_order_14", 122, 123), fix_type, 14),
    )
    lattice_counts = []
    for record, gram, order in explicit:
        sig = signature(gram)
        assert all(sig), f"{record['id']}: lattice is definite {sig}"
        lattice_counts.append({"id": record["id"], "gram": integral(gram), "signature": sig, "isometry_order": order} | {k: v for k, v in record.items() if k != "id"})
    genus_counts = [
        class_count(e, "hermitian_sig_2_4_det_3^5_order_9", 128, 129, (131, 132, 134, 135)),
        class_count(e, "hermitian_unimodular_sig_2_2_2_10_2_18_order_12", 137, 138, (139, 140, 141, 142)),
        class_count(e, "hermitian_sig_2_2_2_6_det_1_to_16_order_4", 144, 145, (146, 147, 148, 149)),
    ]
    at(e, 29, "E6 = root_lattice(:E, 6)")
    definite_counts = []
    for record, gram in (
        (class_count(e, "E6_order_20", 30, 30), e6),
        (class_count(e, "E6_order_9", 31, 31), e6),
        (class_count(e, "E6_genus_order_1", 32, 32), e6),
        (class_count(e, "E6_admissible_triples_2_IpA_2", 34, 34), e6),
        (class_count(e, "E6(2)_admissible_triples_2_IpB_4", 35, 35), 2 * e6),
        (class_count(e, "E6_admissible_triples_3_IpA_2_IpB_4", 36, 36), e6),
        (class_count(e, "E8_order_12_x8-x6+x2-1", 46, 47), e8),
        (class_count(e, "E8_order_12_x8+x6+x2+1", 48, 49), e8),
    ):
        assert not all(signature(gram))
        definite_counts.append({"id": record["id"], "gram": integral(gram), "signature": signature(gram)} | {k: v for k, v in record.items() if k != "id"})
    return {"lattice_class_counts": lattice_counts, "hermitian_genus_class_counts": genus_counts, "definite_class_counts": definite_counts}


def a3_identity_case() -> dict[str, object]:
    """lattices_with_isometry.jl 28-36: the image of O(A3) in O(q_A3) (centralizer of the identity) has order 2."""
    relative = "lattices_with_isometry.jl"
    at(relative, 2, "A3 = root_lattice(:A, 3)")
    at(relative, 28, "L = @inferred integer_lattice_with_isometry(A3)")
    assert at(relative, 36, "@test order(image_centralizer_in_Oq(L)[1]) == 2")
    gram = matrix(ZZ, CartanMatrix(["A", 3]))
    return {
        "id": f"oscar_{Path(relative).stem}_28",
        "gram": integral(gram),
        "signature": signature(gram),
        "isometry": integral(identity_matrix(ZZ, 3)),
        "isometry_order": 1,
        "isometry_convention": "row vectors: v -> v * isometry",
        "centralizer_image_order": 2,
        "source": source(relative, 28, 36),
    }


def main() -> None:
    cases = []
    for relative, first, last, pattern, key in CENTRALIZER_CASES:
        found, text = literals(TESTS / relative, first, last)
        b, g, f = found["B"], found["G"], found["f"]
        lattice_gram = b * g * b.transpose()
        restricted = b.solve_left(b * f)
        assert restricted * b == b * f, f"{relative}:{first}: f does not preserve the lattice"
        gram, isometry = integral(lattice_gram), integral(restricted)
        assert restricted * lattice_gram * restricted.transpose() == lattice_gram, f"{relative}:{first}: f is not an isometry"
        sig = signature(lattice_gram)
        match = re.search(pattern, text)
        assert match, f"{relative}:{first}-{last}: no recorded value matching {pattern}"
        value = True if match.group(1) == "is_bijective" else int(match.group(1))
        order = next(n for n in range(1, 1000) if (restricted**n) == 1)
        cases.append(
            {
                "id": f"oscar_{Path(relative).stem}_{first}",
                "gram": gram,
                "signature": sig,
                "isometry": isometry,
                "isometry_order": order,
                "isometry_convention": "row vectors: v -> v * isometry",
                key: value,
                "source": source(relative, first, last),
            }
        )
    indefinite = [case for case in cases if all(case["signature"])]
    definite = [case for case in cases if not all(case["signature"])]
    assert [case["id"] for case in definite] == ["oscar_lattices_with_isometry_80"], [case["id"] for case in definite]
    relative, first, last = INVOLUTION_CLASSES
    found, text = literals(TESTS / relative, first, last)
    m_gram = found["B"] * found["G"] * found["B"].transpose()
    assert "char_poly=(x-1)^4*(x+1)^6" in text
    local = int(re.search(r"@test length\(r_local\)==(\d+)", text).group(1))
    global_ = int(re.search(r"@test length\(r_global\)==(\d+)", text).group(1))
    involution_classes = {
        "id": f"oscar_{Path(relative).stem}_{first}",
        "genus_representative_gram": integral(m_gram),
        "signature": signature(m_gram),
        "characteristic_polynomial": "(x-1)^4 (x+1)^6",
        "isometry_order": 2,
        "classes_in_genus": global_,
        "local_classes": local,
        "meaning": "isomorphism classes of pairs (L, f) with L in the genus of M and f of the given characteristic polynomial",
        "source": source(relative, first, last),
    }
    definite.append(a3_identity_case())
    counts = class_counts()
    data = {"centralizer_cases": indefinite, "definite_centralizer_cases": definite, "involution_classes": [involution_classes], **counts}
    TARGET.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    print(
        f"wrote {TARGET}: {len(indefinite)} indefinite and {len(definite)} definite centralizer cases, "
        f"involution classes {global_} (local {local}), " + ", ".join(f"{key} {len(value)}" for key, value in counts.items())
    )


if __name__ == "__main__":
    main()
