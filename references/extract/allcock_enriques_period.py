#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/allcock_i_2_10_orbits.json from Allcock math/9905166.

Source (vendored TeX references/vendor/arxiv/9905166/main.tex): the Enriques period
lattice K-hat is I_{2,10} (lines 177-179), and Gamma = Aut(K-hat) (line 185).
Corollary 3: Aut(I_{2,10}) is transitive on norm -1 vectors. Corollary 4: there are
two orbits of primitive isotropic vectors v, distinguished by v^perp / v being I_{1,9}
or II_{1,9}, and two orbits of primitive isotropic planes V, distinguished by
V^perp / V being E8(-1) or I_{0,8}. Every transcribed statement is asserted on its
cited lines.

Run from the repository root: uv run references/extract/allcock_enriques_period.py
"""

import json
from pathlib import Path

SOURCE = Path("references/vendor/arxiv/9905166/main.tex")
TARGET = Path("tests/fixtures/allcock_i_2_10_orbits.json")


def cite(first: int, last: int, *needles: str) -> dict[str, str]:
    text = "\n".join(SOURCE.read_text(encoding="utf-8").splitlines()[first - 1 : last])
    for needle in needles:
        assert needle in text, f"{SOURCE}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "file": str(SOURCE), "lines": f"{first}-{last}"}


def main() -> None:
    lattice = cite(177, 179, r"\Khat\isomorphism", r"\II1,9\oplus\I1,1\isomorphism \I2,{10}")
    group = cite(185, 185, r"\G=\aut\Khat")
    norm_minus_one = cite(214, 220, "Corollary 3", r"transitivity of $\aut\Khat$", r"on the norm $-1$ vectors of $\Khat$")
    isotropic = cite(
        232,
        241,
        "Corollary 4",
        "There are two orbits of primitive isotropic vectors $v$",
        r"$v^\perp/\spanof{v}\isomorphism\I1,9$",
        r"$v^\perp/\spanof{v}\isomorphism\II1,9$",
        "There are two orbits of $2$-dimensional",
        r"$V^\perp/V\isomorphism E_8(-1)$",
        r"$V^\perp/V\isomorphism\I0,8$",
    )
    gram = [[(1 if i < 2 else -1) if i == j else 0 for j in range(12)] for i in range(12)]
    record = {
        "id": "allcock_enriques_period_lattice",
        "lattice": "I_{2,10}",
        "gram": gram,
        "signature": [2, 10],
        "group": "O(I_{2,10})",
        "norm_minus_one_vector_orbits": 1,
        "primitive_isotropic_vector_orbits": [
            {"invariant": "v^perp / v", "isomorphism_class": "I_{1,9}"},
            {"invariant": "v^perp / v", "isomorphism_class": "II_{1,9}"},
        ],
        "isotropic_plane_orbits": [
            {"invariant": "V^perp / V", "isomorphism_class": "E8(-1)"},
            {"invariant": "V^perp / V", "isomorphism_class": "I_{0,8}"},
        ],
        "source": {"lattice": lattice, "group": group, "norm_minus_one": norm_minus_one, "isotropic": isotropic},
    }
    TARGET.write_text(json.dumps(record, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}")


if __name__ == "__main__":
    main()
