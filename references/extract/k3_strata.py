#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/k3_modular_strata.json from local primary sources.

The Scattone boundary counts and types are parsed from a local licensed primary memoir;
Eichler's criterion is checked against the vendored arXiv TeX.

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

Scattone's values are read from a locally supplied memoir copy, not distributed in Git.

Run from the repository root: uv run references/extract/k3_strata.py
"""

import json
import re
from pathlib import Path

ARXIV = Path("references/vendor/arxiv")
SCATTONE = "F. Scattone, On the compactification of moduli spaces for algebraic K3 surfaces, Mem. Amer. Math. Soc. 70 (1987), no. 374"
TARGET = Path("tests/fixtures/k3_modular_strata.json")
MEMOIR = Path("references/vendor/zotero/SF7T3C8G-scattone-1987-memoir/extracted.md")


def cite(path: Path, first: int, last: int, *needles: str) -> dict[str, str]:
    text = "\n".join(path.read_text(encoding="utf-8").splitlines()[first - 1 : last])
    for needle in needles:
        assert needle in text, f"{path}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "file": str(path), "lines": f"{first}-{last}"}


def source_excerpt(lines: list[str], first: int, last: int, *needles: str) -> dict[str, str]:
    excerpt = "\n".join(lines[first - 1 : last])
    for needle in needles:
        assert needle in excerpt, f"{MEMOIR}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "file": str(MEMOIR), "lines": f"{first}-{last}"}


def root_systems(source: str) -> list[str]:
    """Read root-system factors from Scattone's displayed equations."""
    components = re.findall(r"(?:\\mathrm\{([ADE])\}|(?<![a-zA-Z])([ADE]))\s*_\s*\{\s*(\d(?:\s*\d)*)\s*\}(?:\s*\^\s*\{\s*(\d+)\s*\})?", source)
    assert components, "no ADE factors recognized"
    return [f"{a or b}{int(n.replace(' ', ''))}" for a, b, n, exponent in components for _ in range(int(exponent or "1"))]


def main() -> None:
    eichler = cite(ARXIV / "1012.4155/main.tex", 3036, 3041, r"\label{lem:eichler}", "two orthogonal isotropic planes", "determined by two invariants")
    if not MEMOIR.is_file():
        raise FileNotFoundError(f"Scattone primary source is not distributed: supply your local reading copy at {MEMOIR}")
    lines = MEMOIR.read_text(encoding="utf-8").splitlines()
    boundary = source_excerpt(lines, 3740, 3770, "adding one point p and h(k) curves", "boundary curves")
    degree_two_types = source_excerpt(lines, 3850, 3860, r"\tag{6.2.1}", r"\mathrm{A} _ {1}")
    degree_four = source_excerpt(lines, 3932, 3948, "Degree four", "generalized) types")
    correspondence = source_excerpt(lines, 3224, 3250, "Proposition 5.4.7", "Corollary 5.4.8", "if k is squarefree")
    h_count_text = "\n".join(lines[3757:3773])
    h_counts = {int(k): int(n) for k, n in re.findall(r"h\((\d+)\)\s*=\s*(\d+)", h_count_text)}
    assert h_counts == {1: 4, 2: 9}, f"unexpected Scattone h(k): {h_counts}"
    degree_two_display = "\n".join(lines[3852:3857])
    two_types = ["+".join(root_systems(item)) for item in re.findall(r"(.*?)(?:,|\\tag\{6\.2\.1\})", degree_two_display) if "\\mathrm" in item]
    assert len(two_types) == h_counts[1], two_types
    degree_four_display = "\n".join(lines[3939:3947])
    display_match = re.search(r"\\begin\{array\}\{l\}(.*?)\\end\{array\}", degree_four_display, re.S)
    assert display_match is not None, "missing degree-four type array"
    four_types = []
    for label in re.split(r",|\\\\", display_match.group(1)):
        factors = root_systems(label) if re.search(r"(?:\\mathrm|[ADE] _)", label) else []
        minus_four = "<-4>" if re.search(r"<\s*-\s*4\s*>", label) else None
        if minus_four:
            factors.insert(0, minus_four)
        if factors:
            four_types.append("+".join(factors))
    assert len(four_types) == h_counts[2], four_types
    # Canonical form of each ADE multiset, independent of printed factor order.
    two_types = sorted("+".join(sorted(t.split("+"))) for t in two_types)
    four_types = sorted("+".join(sorted(t.split("+"))) for t in four_types)
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
            "source": {"boundary": boundary, "root_types": degree_two_types},
        },
        "degree_four_polarized_k3": {
            "id": "k3_degree_four_modular_variety",
            "construction": "2U + 2E8(-1) + <-4>",
            "signature": [2, 19],
            "rank": 21,
            "det": 4,
            "baily_borel_boundary": {"type_ii_curves": h_counts[2]},
            "type_ii_generalised_types": four_types,
            "source": {"boundary": boundary, "generalised_types": degree_four},
        },
    }
    TARGET.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}")


if __name__ == "__main__":
    main()
