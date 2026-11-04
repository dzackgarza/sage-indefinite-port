#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract K3 orbit and boundary strata from the committed source TeX.

The primary Scattone memoir is not distributed. The published boundary
classifications are reproduced in the vendored arXiv manuscripts and extracted
from those cited mathematical statements instead.

- K3 lattice: Eichler's criterion, Gritsenko--Hulek--Sankaran arXiv:1012.4155,
  Lemma lem:eichler. For unimodular L the discriminant group is trivial, so the
  orbit of a primitive vector is determined by its norm.
- Degree 2 and degree 4: Scattone, On the compactification of moduli spaces for
  algebraic K3 surfaces (Mem. AMS 374, 1987). Section 6: the Baily--Borel boundary
  of D_k/O(L_k) adds one point and h(k) curves meeting there, with h(1) = 4 and
  h(2) = 9; (6.2.1) gives the root types of the four degree-2 Type II components,
  and Section 6.3 the generalised types of the nine degree-4 components.
- The group: Scattone Proposition 5.4.7 and Corollary 5.4.8(2). For squarefree k
  the O(L_k)-orbits of primitive isotropic planes E are in bijection with G(k) by
  E -> E^perp/E, so there are h(k) of them.

The published Scattone classifications are cross-checked in the committed
secondary papers; no copyrighted local reading copy is required.

Run from the repository root: uv run references/extract/k3_strata.py
"""

import json
import re
from pathlib import Path

ARXIV = Path("references/vendor/arxiv")
SCATTONE = "F. Scattone, On the compactification of moduli spaces for algebraic K3 surfaces, Mem. Amer. Math. Soc. 70 (1987), no. 374"
TARGET = Path("tests/fixtures/k3_modular_strata.json")
DEGREE_TWO = ARXIV / "1205.3144/ksba.tex"
DEGREE_FOUR = ARXIV / "2502.04301/main.tex"


def cite(path: Path, first: int, last: int, *needles: str) -> dict[str, str]:
    text = "\n".join(path.read_text(encoding="utf-8").splitlines()[first - 1 : last])
    for needle in needles:
        assert needle in text, f"{path}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "file": str(path), "lines": f"{first}-{last}"}


def root_systems(source: str) -> list[str]:
    """Read root-system factors from Scattone's displayed equations."""
    components = re.findall(r"(?:\\mathrm\{([ADE])\}|(?<![a-zA-Z])([ADE]))\s*_\s*(?:\{\s*(\d+(?:\s*\d)*)\s*\}|(\d+))(?:\s*\^\s*\{\s*(\d+)\s*\})?", source)
    assert components, "no ADE factors recognized"
    return [f"{a or b}{int((n or bare).replace(' ', ''))}" for a, b, n, bare, exponent in components for _ in range(int(exponent or "1"))]


def main() -> None:
    eichler = cite(ARXIV / "1012.4155/main.tex", 3036, 3041, r"\label{lem:eichler}", "two orthogonal isotropic planes", "determined by two invariants")
    two_source = DEGREE_TWO.read_text(encoding="utf-8").splitlines()
    four_source = DEGREE_FOUR.read_text(encoding="utf-8").splitlines()
    two_statement = two_source[262]  # degree-two Type II labels in the arithmetic discussion
    assert "For degree $2$" in two_statement and "four Type II components" in two_statement
    two_types = []
    for label in re.findall(r"\$([^$]+)\$", two_statement.split("correspond to the root lattices", 1)[1]):
        if not re.search(r"[ADE]_", label):
            continue
        multiplicity = re.search(r"(\d+)([ADE])_(?:\{(\d+)\}|(\d+))", label)
        if multiplicity:
            count, kind, subscript, bare = multiplicity.groups()
            label = label.replace(multiplicity.group(), "")
            factors = [f"{kind}{subscript or bare}"] * int(count)
        else:
            factors = []
        two_types.append("+".join(factors + root_systems(label)))
    assert len(two_types) == 4, two_types
    two_types = sorted("+".join(sorted(t.split("+"))) for t in two_types)
    assert len(set(two_types)) == 4

    two_boundary = cite(DEGREE_TWO, 252, 256, "four curves", "single point")
    degree_two_types = cite(DEGREE_TWO, 262, 264, "four Type II components", "A_{17}")
    four_statement = four_source[292]  # Scattone's nine generalized types reproduced verbatim
    assert "A_{11}+E_6" in four_statement and "D_{17}" in four_statement
    four_types = []
    for label in re.split(r",\s*|,\\\\", four_statement):
        factors = root_systems(label) if re.search(r"[ADE]_", label) else []
        if re.search(r"\\langle\s*-4\s*\\rangle", label):
            factors.insert(0, "<-4>")
        if factors:
            four_types.append("+".join(sorted(factors)))
    four_types = sorted(four_types)
    assert len(four_types) == 9 and len(set(four_types)) == 9, four_types
    degree_four = cite(DEGREE_FOUR, 289, 296, "A_{11}+E_6", "9 boundary curves")
    correspondence = cite(DEGREE_FOUR, 273, 282, "5.4.7", "I^\\perp/I", "boundary curve")
    h_counts = {1: len(two_types), 2: len(four_types)}
    record = {
        "type_ii_correspondence": {
            "statement": (
                "for squarefree k, E -> E^perp/E is a bijection from the O(Lambda_{2k})-orbits of primitive isotropic planes E of Lambda_{2k} "
                "onto the genus G(k) of <-2k> + 2E8(-1), so there are h(k) orbits"
            ),
            "group": "O(Lambda_{2k})",
            "source": correspondence,
        },
        "k3_unimodular_lattice": {
            "id": "k3_lattice_lambda",
            "construction": "3U + 2E8(-1)",
            "signature": [3, 19],
            "rank": 22,
            "det": -1,
            "is_even": True,
            "is_unimodular": True,
            "primitive_vector_orbits": "one orbit per represented norm (trivial discriminant group)",
            "tested_represented_norms": [0, -2, 2, 4, -4],
            "source": eichler,
        },
        "degree_two_polarized_k3": {
            "id": "k3_degree_two_modular_variety",
            "construction": "2U + 2E8(-1) + <-2>",
            "signature": [2, 19],
            "rank": 21,
            "det": 2,
            "baily_borel_boundary": {"type_iii_points": 1, "type_ii_curves": h_counts[1], "point_curve_incidences": h_counts[1]},
            "type_ii_root_types": two_types,
            "source": {"boundary": two_boundary, "root_types": degree_two_types},
        },
        "degree_four_polarized_k3": {
            "id": "k3_degree_four_modular_variety",
            "construction": "2U + 2E8(-1) + <-4>",
            "signature": [2, 19],
            "rank": 21,
            "det": 4,
            "baily_borel_boundary": {"type_ii_curves": h_counts[2]},
            "type_ii_generalised_types": four_types,
            "source": {"boundary": degree_four, "generalised_types": degree_four},
        },
    }
    TARGET.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}")


if __name__ == "__main__":
    main()
