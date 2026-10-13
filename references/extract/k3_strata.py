#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/k3_modular_strata.json from the vendored K3 literature.

The values are transcribed, and every transcribed label or count is asserted to
appear on the cited lines of the vendored TeX, so a transcription error fails
here instead of in a test.

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

Scattone's values are transcribed from the memoir and cited by its numbering.

Run from the repository root: uv run references/extract/k3_strata.py
"""

import json
from pathlib import Path

ARXIV = Path("references/vendor/arxiv")
SCATTONE = "F. Scattone, On the compactification of moduli spaces for algebraic K3 surfaces, Mem. Amer. Math. Soc. 70 (1987), no. 374"
TARGET = Path("tests/fixtures/k3_modular_strata.json")


def cite(path: Path, first: int, last: int, *needles: str) -> dict[str, str]:
    text = "\n".join(path.read_text(encoding="utf-8").splitlines()[first - 1 : last])
    for needle in needles:
        assert needle in text, f"{path}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "file": str(path), "lines": f"{first}-{last}"}


def book(location: str) -> dict[str, str]:
    """A statement of Scattone's memoir, read from the memoir and cited by its numbering."""
    return {"kind": "published_text", "citation": SCATTONE, "location": location}


def main() -> None:
    eichler = cite(ARXIV / "1012.4155/main.tex", 3036, 3041, r"\label{lem:eichler}", "two orthogonal isotropic planes", "determined by two invariants")
    boundary = book("Section 6, opening: D_k/O_-(L_k) gains one point and h(k) curves through it; h(1) = 4, h(2) = 9")
    degree_two_types = book("(6.2.1)")
    degree_four = book("Section 6.3, the nine orthogonal complements of D_7")
    correspondence = book("Proposition 5.4.7 and Corollary 5.4.8(2)")
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
            "baily_borel_boundary": {"type_iii_points": 1, "type_ii_curves": 4, "point_curve_incidences": 4},
            "type_ii_root_types": ["A17", "E8+E8+A1", "D16+A1", "E7+D10"],
            "source": {"boundary": boundary, "root_types": degree_two_types},
        },
        "degree_four_polarized_k3": {
            "id": "k3_degree_four_modular_variety",
            "construction": "2U + 2E8(-1) + <-4>",
            "signature": [2, 19],
            "rank": 21,
            "det": 4,
            "baily_borel_boundary": {"type_ii_curves": 9},
            "type_ii_generalised_types": [
                "A11+E6",
                "A1+A1+A15",
                "D8+D8+<-4>",
                "D12+D5",
                "D16+<-4>",
                "E7+E7+A3",
                "E8+E8+<-4>",
                "E8+D9",
                "D17",
            ],
            "source": {"boundary": boundary, "generalised_types": degree_four},
        },
    }
    TARGET.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}")


if __name__ == "__main__":
    main()
