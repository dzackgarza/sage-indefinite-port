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
- Degree 2: Laza arXiv:1205.3144, Theorem (citing Scattone 6.2): the Baily--Borel
  boundary is four Type II curves meeting in one Type III point; Shah's table
  there gives the root types of the four Type II components.
- Degree 4: Jones arXiv:2502.04301, Theorem thm:scattone lattices (citing
  Scattone 6.3): nine Type II boundary curves with their generalised types.

Run from the repository root: uv run references/extract/k3_strata.py
"""

import json
from pathlib import Path

ARXIV = Path("references/vendor/arxiv")
TARGET = Path("tests/fixtures/k3_modular_strata.json")


def cite(relative: str, first: int, last: int, *needles: str) -> dict[str, str]:
    path = ARXIV / relative
    text = "\n".join(path.read_text(encoding="utf-8").splitlines()[first - 1 : last])
    for needle in needles:
        assert needle in text, f"{path}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "file": str(path), "lines": f"{first}-{last}"}


def main() -> None:
    eichler = cite("1012.4155/main.tex", 3036, 3041, r"\label{lem:eichler}", "two orthogonal isotropic planes", "determined by two invariants")
    degree_two_boundary = cite("1205.3144/ksba.tex", 256, 256, "four curves", "single point")
    degree_two_types = cite("1205.3144/ksba.tex", 194, 205, "$A_{17}$", "$E_8^{2}+A_1$", "$D_{16}+A_1$", "$E_7+D_{10}$")
    degree_four = cite(
        "2502.04301/main.tex",
        289,
        296,
        r"\cite[Section 6.3]{Sca87}",
        "A_{11}+E_6",
        "A_1+A_1+A_{15}",
        r"D_8+D_8+\langle -4 \rangle",
        "D_{12}+D_5",
        r"D_{16}+\langle -4 \rangle",
        "E_7+E_7+A_3",
        r"E_8+E_8+\langle -4 \rangle",
        "E_8+D_9",
        "D_{17}",
        "9 boundary curves",
    )
    correspondence = cite(
        "2502.04301/main.tex",
        276,
        281,
        r"\cite[5.4.7]{Sca87} Let $I \subset \Lambda_{2k}$ be a rank 2 isotropic sublattice",
        r"I \text{ } \operatorname{modulo} \text{ }O(\Lambda_{2k}) \leftrightarrow C_I",
    )
    record = {
        "type_ii_correspondence": {
            "statement": (
                "boundary curves C_I of F_{2k} correspond bijectively to rank-2 isotropic sublattices I of Lambda_{2k} modulo O(Lambda_{2k}), and to I^perp/I in the genus G(k)"
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
            "source": {"boundary": degree_two_boundary, "root_types": degree_two_types},
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
            "source": degree_four,
        },
    }
    TARGET.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}")


if __name__ == "__main__":
    main()
