"""Extract tests/fixtures/isometry_centralizers.json from OSCAR's vendored tests.

Sources: references/vendor/Oscar.jl@d135b70b/test/NumberTheory/QuadFormAndIsom.
OSCAR is an independent implementation, so its recorded outputs are external
oracles for this port's centralizer computations. Each case is read mechanically
from the cited line range: the ambient basis matrix B, ambient Gram matrix G and
ambient isometry f (OSCAR acts on row vectors), and the recorded value from the
@test line in that range. The script computes the lattice's Gram matrix B G B^T
and the isometry in the lattice basis (the solution f_L of f_L B = B f), and
asserts that f_L is integral, preserves the Gram matrix, and that the lattice is
indefinite; definite cases belong to the research preamble.

Run from the repository root under Sage's Python:
    "$(dirname $(sage -c 'import sys; print(sys.executable)'))/python3" references/extract/oscar_centralizers.py
"""

import json
import re
from fractions import Fraction
from pathlib import Path

# Importing from sage.all runs Sage's session startup.
from sage.all import AA, QQ, ZZ, matrix

TESTS = Path("references/vendor/Oscar.jl@d135b70b/test/NumberTheory/QuadFormAndIsom")
TARGET = Path("tests/fixtures/isometry_centralizers.json")
LITERAL = re.compile(r"(\w+)\s*=\s*matrix\(QQ,\s*(\d+),\s*(\d+)\s*,\s*\[(.*?)\]\)", re.S)

# (file, first line, last line, the @test pattern carrying the recorded value, record key)
CENTRALIZER_CASES = (
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
        assert sig[0] > 0 and sig[1] > 0, f"{relative}:{first}: lattice is definite {sig}"
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
    TARGET.write_text(json.dumps({"centralizer_cases": cases, "involution_classes": [involution_classes]}, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}: {len(cases)} centralizer cases, involution classes {global_} (local {local})")


if __name__ == "__main__":
    main()
