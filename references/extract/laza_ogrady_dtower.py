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

Gamma(N) is defined in Laza--O'Grady arXiv:1607.01324 (vendored; LOG3 line 642 cites
its Prop. 1.2.3): Gamma_xi = {phi in O^+(Lambda) : phi(xi) = xi} for a decoration xi in
A_Lambda of square 1 mod 2Z, recorded under "gamma" with every clause asserted on its line.

Run from the repository root: uv run references/extract/laza_ogrady_dtower.py
"""

import json
from pathlib import Path
from typing import TypedDict

SOURCE = Path("references/vendor/arxiv/1801.04845/LOG3-2018-08-24.tex")
TARGET = Path("tests/fixtures/dtower_boundaries.json")


class TypeCount(TypedDict):
    N: int
    type_ii_components: int
    type_ii_from_genus_of_d: int
    type_ii_from_even_unimodular: int
    type_iii_components: int
    group_note: str
    derivation: str
    source: dict[str, dict[str, str]]


class FHatRow(TypedDict):
    label: str
    type: str
    dimension_in_F_hat: int
    geometric_meaning: str
    quartic_case: str
    source: dict[str, str | int]


LOG1 = Path("references/vendor/arxiv/1607.01324/bir-quartics-arxiv.tex")


def cite_log1(first: int, last: int, *needles: str) -> dict[str, str]:
    text = "\n".join(LOG1.read_text(encoding="utf-8").splitlines()[first - 1 : last])
    for needle in needles:
        assert needle in text, f"{LOG1}:{first}-{last} does not contain {needle!r}"
    return {"kind": "published_text", "file": str(LOG1), "lines": f"{first}-{last}"}


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


def genus_of_d() -> dict[int, list[str]]:
    """Theorem thm:classifydn: root sublattices classifying the genus of D_n, 1 <= n <= 18."""
    source_text = "\n".join(SOURCE.read_text(encoding="utf-8").splitlines()[4463:4515])
    stated = {
        1: ["empty (the single class D1 = <-4>, whose root sublattice is empty)"],
        9: ["D9", "E8"],
        13: ["D13", "D5+E8", "D12"],
        14: ["D14", "D6+E8", "D12+D2"],
        15: ["D15", "D7+E8", "D12+D3", "A15", "E7+E7"],
        16: ["D16", "D8+E8", "D12+D4", "D2+E7+E7", "D8+D8", "A15"],
        17: ["D17", "D9+E8", "D12+D5", "D3+E7+E7", "A15+D2", "A11+E6", "D8+D8", "D16", "E8+E8"],
        18: ["D18", "D10+E8", "D12+D6", "D4+E7+E7", "A15+D3", "D8+D8+D2", "D16+D2", "D2+E8+E8", "A11+E6", "D6+D6+D6", "A9+A9", "D10+E7+A1", "A17+A1"],
    }
    for n in range(2, 9):
        stated[n] = [f"D{n}"]
    for n in range(10, 13):
        stated[n] = [f"D{n}", f"D{n - 8}+E8"]
    for needle in (
        r"$n=1$: $\es$",
        r"$2\le n\le 8$: $D_n$",
        r"$n=9$: $D_9$, $E_8$",
        r"$10\le n \le 12$: $D_n$, $D_{n-8}\oplus E_8$",
        r"$n=13$: $D_{13}$, $D_{5}\oplus E_8$, $D_{12}$",
        r"$n=14$: $D_{14}$, $D_{6}\oplus E_8$, $D_{12}\oplus D_{2}$",
        r"$n=15$: $D_{15}$, $D_7\oplus E_8$, $D_{12}\oplus D_3$,  $A_{15}$, $(E_7)^2$",
        r"$n=16$: $D_{16}$, $D_8\oplus E_8$, $D_{12}\oplus D_4$, $D_2\oplus (E_7)^2$, $(D_8)^2$, $A_{15}$",
        r"$A_{11}\oplus E_6$,  $(D_8)^2$, $D_{16}$, $(E_8)^2$",
        r"$(A_9)^2$, $D_{10}\oplus E_7\oplus A_1$, $A_{17}\oplus A_1$",
    ):
        assert needle in source_text, needle
    assert sorted(stated) == list(range(1, 19))
    return stated


def type_counts(genus: dict[int, list[str]], definition: dict[str, str], group: dict[str, str]) -> list[TypeCount]:
    """Theorem thm:bbtype2comp with the unimodular lattices of line 4460, for 3 <= N <= 20."""
    theorem = cite(
        4604,
        4616,
        r"Let $3\le N\le 20$. The number of Type II  components",
        "with the exception of $N=14$, in this case  there are three Type II components labeled by $D_{12}$",
        r"If $N\equiv 2 \pmod 8$, there are two kinds of Type II components",
    )
    unimodular = cite(4460, 4460, "namely  $E_8$", r"namely $E_8\oplus E_8$, and the unique (up to isomorphism) unimodular overlattice of $D_{16}$")
    type_iii = cite(4391, 4391, r"the number of Type III boundary components is $1$ if    $N\not\equiv 2\pmod{8}$, and $2$ if $N\equiv 2\pmod{8}$")
    even_unimodular = {8: ["E8"], 16: ["E8+E8", "D16+"]}
    counts: list[TypeCount] = []
    for big_n in range(3, 21):
        n = big_n - 2
        first_kind = len(genus[n]) + (2 if big_n == 14 else 0)
        second_kind = len(even_unimodular[n]) if big_n % 8 == 2 else 0
        counts.append(
            {
                "N": big_n,
                "type_ii_components": first_kind + second_kind,
                "type_ii_from_genus_of_d": first_kind,
                "type_ii_from_even_unimodular": second_kind,
                "type_iii_components": 2 if big_n % 8 == 2 else 1,
                "group_note": "Gamma(N) = O^+(Lambda_N) unless n = 6 mod 8, then an index-3 subgroup (line 642)",
                "derivation": (
                    "Type II = classes in the genus of D_{N-2} (two extra components for D12 when N = 14), plus the even unimodular lattices of rank N-2 when N = 2 mod 8"
                ),
                "source": {"lattice": definition, "group": group, "theorem": theorem, "unimodular": unimodular, "type_iii": type_iii},
            }
        )
    return counts


def f_hat_table() -> list[FHatRow]:
    """Table tabletype2: the eight Type II components of F(18)* and their pre-images in F-hat."""
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    label = next(i for i, line in enumerate(lines) if r"\label{tabletype2}" in line)
    begin = max(i for i in range(label) if r"\begin{tabular}" in lines[i])
    rows: list[FHatRow] = []
    for index in range(begin + 1, label):
        line = lines[index].strip()
        if not line.endswith(r"\\") or line.startswith("Label"):
            continue
        cells = [cell.strip() for cell in line[:-2].split("&")]
        assert len(cells) == 5, f"{SOURCE}:{index + 1}"
        rows.append(
            {
                "label": cells[0],
                "type": cells[1],
                "dimension_in_F_hat": int(cells[2].strip("$")),
                "geometric_meaning": cells[3],
                "quartic_case": cells[4],
                "source": {"kind": "published_table", "file": str(SOURCE), "table": "tabletype2", "line": index + 1},
            }
        )
    return rows


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
    genus = genus_of_d()
    counts = type_counts(genus, definition, group)
    for record in records:
        derived = next(c for c in counts if c["N"] == record["N"])
        assert (derived["type_ii_components"], derived["type_iii_components"]) == (record["type_ii_curves"], record["type_iii_points"]), record["N"]
    k3 = json.loads(Path("tests/fixtures/k3_modular_strata.json").read_text())
    n19 = next(c for c in counts if c["N"] == 19)
    assert n19["type_ii_components"] == k3["degree_four_polarized_k3"]["baily_borel_boundary"]["type_ii_curves"] == 9, "N = 19 must match the quartic K3 count"
    f_hat = f_hat_table()
    n18 = next(c for c in counts if c["N"] == 18)
    assert n18["type_ii_components"] == len(f_hat) == 8
    labels_a = sorted(row["label"] for row in f_hat if row["type"] == "a")
    assert len(labels_a) == len(genus[16]) == 6
    gamma = {
        "group": "Gamma(N) = Gamma_xi = {phi in O^+(Lambda_N) : phi(xi) = xi}",
        "decoration": "xi in A_{Lambda_N} with q(xi) = 1 mod 2Z; unique unless N = 6 mod 8",
        "index_in_O_plus": "1, or 3 when N = 6 mod 8",
        "source": {
            "gamma_n": cite(642, 642, r"$\Gamma(N)=O^{+}(\Lambda_N)$ if $n\not\equiv 6\pmod{8}$", "Prop.~1.2.3 ibid"),
            "decoration": cite_log1(619, 619, r"a \emph{decoration} of  $\Lambda$  is an element $\xi\in A_{\Lambda}$  of square $1$ (modulo $2\ZZ$)"),
            "uniqueness": cite_log1(626, 626, r"it is unique unless $N\equiv 6 \pmod 8$"),
            "definition": cite_log1(698, 698, r"\Gamma_\xi:=\{\phi\in O^+(\Lambda)\mid \phi(\xi)=\xi\}"),
            "index": cite_log1(758, 758, r"$\Gamma_\xi= O^+(\Lambda)$ unless $N\equiv 6\pmod 8$, in which case $\Gamma_\xi< O^+(\Lambda)$ is of index $3$"),
        },
    }
    data = {
        "gamma": gamma,
        "boundary_pictures": records,
        "type_counts": counts,
        "genus_of_d": [{"n": n, "root_sublattices": genus[n]} for n in sorted(genus)],
        "f18_hat_type_ii": f_hat,
    }
    TARGET.write_text(json.dumps(data, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}: {len(records)} pictures, Type II/III counts for N = 3..20 (match the pictures and the 8 rows of the N = 18 table)")


if __name__ == "__main__":
    main()
