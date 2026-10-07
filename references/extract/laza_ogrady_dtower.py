#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# ///
"""Extract tests/fixtures/dtower_boundaries.json from Laza--O'Grady arXiv:1801.04845.

Source (vendored TeX references/vendor/arxiv/1801.04845/LOG3-2018-08-24.tex):

- Lambda_N = U^2 + D_{N-2}, with D_n negative definite (line 625);
- F(N) = Gamma(N) \\ D^+(N) with Gamma(N) = O^+(Lambda_N) unless n = 6 mod 8 (line
  642; the text leaves n undefined, but for N = 9, 10, 11 neither reading n = N nor
  n = N - 2 is 6 mod 8, so Gamma(N) = O^+(Lambda_N) for the cases recorded here);
- the Baily--Borel boundaries of F(9), F(10), F(11) (diagram bbpicture, lines
  4629-4631): Type III points (isotropic lines) and Type II curves (isotropic planes)
  with their labels and incidences.

Every transcribed label is asserted on the cited lines.

Run from the repository root: uv run references/extract/laza_ogrady_dtower.py
"""

import json
from pathlib import Path

SOURCE = Path("references/vendor/arxiv/1801.04845/LOG3-2018-08-24.tex")
TARGET = Path("tests/fixtures/dtower_boundaries.json")


def cite(first: int, last: int, *needles: str) -> dict[str, str]:
    text = "\n".join(SOURCE.read_text(encoding="utf-8").splitlines()[first - 1 : last])
    for needle in needles:
        assert needle in text, f"{SOURCE}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "file": str(SOURCE), "lines": f"{first}-{last}"}


def d_lattice(n: int) -> list[list[int]]:
    """The negative definite D_n: chain 0-1-...-(n-2), node n-1 attached to node n-3."""
    gram = [[-2 if i == j else 0 for j in range(n)] for i in range(n)]
    edges = [(i, i + 1) for i in range(n - 2)] + [(n - 3, n - 1)]
    for i, j in edges:
        gram[i][j] = gram[j][i] = 1
    return gram


def lambda_n(big_n: int) -> list[list[int]]:
    blocks = [[[0, 1], [1, 0]], [[0, 1], [1, 0]], d_lattice(big_n - 2)]
    size = sum(len(b) for b in blocks)
    gram = [[0] * size for _ in range(size)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, entry in enumerate(row):
                gram[offset + i][offset + j] = entry
        offset += len(block)
    return gram


def main() -> None:
    definition = cite(625, 625, r"\Lambda_N$ be the lattice $U^2\oplus D_{N-2}", "negative definite")
    group = cite(642, 642, r"\Gamma(N)=O^{+}(\Lambda_N)$ if $n\not\equiv 6\pmod{8}")
    pictures = {
        9: (cite(4629, 4629, r"N=9: & &&\bullet^{III_a}\ar@{-}[r]&\circ^{II(D_7)}"), ["III_a"], ["D7"], [("III_a", "D7")]),
        10: (
            cite(4630, 4630, r"N=10: & \bullet^{III_b}\ar@{-}[r]&\circ^{II(E_8)}\ar@{-}[r]&\bullet^{III_a}\ar@{-}[r]&\circ^{II(D_8)}"),
            ["III_a", "III_b"],
            ["E8", "D8"],
            [("III_b", "E8"), ("III_a", "E8"), ("III_a", "D8")],
        ),
        11: (
            cite(4631, 4631, r"N=11: &&\circ^{II(E_8\oplus D_1)}\ar@{-}[r]&\bullet^{III_a}\ar@{-}[r]&\circ^{II(D_9)}"),
            ["III_a"],
            ["E8+D1", "D9"],
            [("III_a", "E8+D1"), ("III_a", "D9")],
        ),
    }
    records = []
    for big_n, (picture, type_iii, type_ii, incidences) in pictures.items():
        assert big_n % 8 != 6 and (big_n - 2) % 8 != 6
        records.append(
            {
                "id": f"lambda_{big_n}",
                "N": big_n,
                "lattice": f"U^2 + D_{big_n - 2} (D negative definite)",
                "gram": lambda_n(big_n),
                "signature": [2, big_n],
                "group": "O^+(Lambda_N)",
                "type_iii_points": len(type_iii),
                "type_ii_curves": len(type_ii),
                "type_ii_labels": type_ii,
                "incidences": [list(edge) for edge in incidences],
                "source": {"lattice": definition, "group": group, "boundary": picture},
            }
        )
    TARGET.write_text(json.dumps(records, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {len(records)} boundaries to {TARGET}")


if __name__ == "__main__":
    main()
