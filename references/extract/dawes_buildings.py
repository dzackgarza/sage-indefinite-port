#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/dawes_buildings.json from Dawes's vendored papers.

Each transcribed fact is asserted to appear on the cited lines of the vendored TeX.
Dawes assumes all root lattices negative definite (arXiv:2205.10601, line 158), so
A2 below is the negative definite A2 and every lattice has signature (2, n).

- arXiv:2205.10601, section "Examples": the Tits building of the stable
  orthogonal group of 2U + A2 has one point, one curve and one edge; the buildings
  of O^+ and of the stable group of 2U + <-2> + <-6> are the same path with two
  points, two curves and three edges; and the stable groups of 2U(2) + A2,
  U + U(2) + A2, 2U + A2 and O^+(2U + A2) form a chain with indices 20, 27 and 2.
- arXiv:2108.06236, Theorem L2boundarythm: for L_2 = 2U + <-2> + <-6> and
  Gamma_2 = {g in O^+(L_2) : g v* = v* mod L_2}, v generating <-2>, the boundary
  has two curves, three points and four point-curve incidences.

Run from the repository root: uv run references/extract/dawes_buildings.py
"""

import json
from pathlib import Path

ARXIV = Path("references/vendor/arxiv")
TARGET = Path("tests/fixtures/dawes_buildings.json")
U = [[0, 1], [1, 0]]


def cite(relative: str, first: int, last: int, *needles: str) -> dict[str, str]:
    path = ARXIV / relative
    text = "\n".join(path.read_text(encoding="utf-8").splitlines()[first - 1 : last])
    for needle in needles:
        assert needle in text, f"{path}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "file": str(path), "lines": f"{first}-{last}"}


def block_sum(*blocks: list[list[int]]) -> list[list[int]]:
    size = sum(len(block) for block in blocks)
    gram = [[0] * size for _ in range(size)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, entry in enumerate(row):
                gram[offset + i][offset + j] = entry
        offset += len(block)
    return gram


def main() -> None:
    paper = "2205.10601/orbits_in_lattices.tex"
    convention = cite(paper, 158, 158, "we will assume all roots lattices are negative definite")
    a2 = [[-2, 1], [1, -2]]
    records = [
        {
            "id": "dawes_2U_A2_stable",
            "lattice": "2U + A2",
            "gram": block_sum(U, U, a2),
            "signature": [2, 4],
            "group": "stable orthogonal group (discriminant kernel of O^+)",
            "building": {"points": 1, "curves": 1, "edges": 1},
            "source": {
                "convention": convention,
                "statement": cite(paper, 1310, 1310, r"Figure \ref{stable_2U_A2_tits}"),
                "figure": cite(paper, 1349, 1364, "(-3,0)--(3,0)", r"\label{stable_2U_A2_tits}"),
            },
        },
        {
            "id": "dawes_2U_minus2_minus6",
            "lattice": "2U + <-2> + <-6>",
            "gram": block_sum(U, U, [[-2]], [[-6]]),
            "signature": [2, 4],
            "group": "O^+ and the stable orthogonal group (same building)",
            "building": {"points": 2, "curves": 2, "edges": 3},
            "source": {
                "statement": cite(paper, 1392, 1392, r"\widetilde{\opn{O}}^+(L))$ and $\bc(\opn{O}^+(L))$ are given by Figure \ref{gkbuilding}"),
                "figure": cite(paper, 1426, 1444, "(-6,0)--(6,0)", r"\label{gkbuilding}"),
            },
        },
        {
            "id": "dawes_L2_gamma2",
            "lattice": "L_2 = 2U + <-2> + <-6>",
            "gram": block_sum(U, U, [[-2]], [[-6]]),
            "signature": [2, 4],
            "group": "Gamma_2 = {g in O^+(L_2) : g v* = v* mod L_2}, v generating <-2>",
            "building": {"points": 3, "curves": 2, "edges": 4},
            "source": {
                "definition": cite("2108.06236/main.tex", 204, 212, r"L_{2d} = 2U \op \la -2d \ra \op \la -6 \ra", r"\Gamma_{2d} = \{ g \in \opn{O}^+(L) \mid g \underline{v}^* \equiv \underline{v}^* \bmod{ L} \}"),
                "theorem": cite("2108.06236/main.tex", 1044, 1047, r"\label{L2boundarythm}", "curves $\\cc_1$ and $\\cc_2$", "points $P_1$, $P_2$, $P_3$", r"$\overline{\cc}_1 \cap P_1$, $\overline{\cc}_1 \cap P_2$, $\overline{\cc}_2 \cap P_2$ and $\overline{\cc}_2 \cap P_3$"),
            },
        },
    ]
    index_chain = {
        "id": "dawes_2U2_A2_index_chain",
        "chain": [
            "stable orthogonal group of 2U(2) + A2",
            "stable orthogonal group of U + U(2) + A2",
            "stable orthogonal group of 2U + A2",
            "O^+(2U + A2)",
        ],
        "indices": [20, 27, 2],
        "total_index": 1080,
        "source": cite(paper, 1501, 1513, "$\\vert G_2:G_1 \\vert = 1080$", "20, 27 and 2, respectively"),
    }
    assert index_chain["indices"][0] * index_chain["indices"][1] * index_chain["indices"][2] == index_chain["total_index"]
    TARGET.write_text(json.dumps({"buildings": records, "index_chains": [index_chain]}, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}")


if __name__ == "__main__":
    main()
