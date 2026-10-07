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
  points, two curves and three edges; the building of the stable orthogonal group of
  2U(2) + A2 (figure 2u2a2building) is transcribed from its TikZ source; and the
  stable groups of 2U(2) + A2,
  U + U(2) + A2, 2U + A2 and O^+(2U + A2) form a chain with indices 20, 27 and 2.
- arXiv:2108.06236, Theorem L2boundarythm: for L_2 = 2U + <-2> + <-6> and
  Gamma_2 = {g in O^+(L_2) : g v* = v* mod L_2}, v generating <-2>, the boundary
  has two curves, three points and four point-curve incidences.

For figure 2u2a2building, black nodes are curves and white nodes are points (as in the
other building figures of the paper). The outer circle joins consecutive nodes at
k * 30 degrees; each chord "edge [bend right=30]" joins two white nodes and passes
through one inner black node, which the script determines by evaluating the midpoint
of the TikZ Bezier curve (control distance 0.3915 * chord length) and asserting it lies
within 0.01 cm of exactly one inner node. The chord contributes the two point-curve
edges at that inner node. Labels are the figure's node labels.

Run from the repository root: uv run references/extract/dawes_buildings.py
"""

import json
import math
import re
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


def polar(degrees: float, radius: float) -> tuple[float, float]:
    return (radius * math.cos(math.radians(degrees)), radius * math.sin(math.radians(degrees)))


def bend_right_midpoint(start: tuple[float, float], end: tuple[float, float], bend: float) -> tuple[float, float]:
    """Midpoint of TikZ's 'edge [bend right=bend]' cubic Bezier from start to end."""
    length = math.dist(start, end)
    direction = math.atan2(end[1] - start[1], end[0] - start[0])
    out_angle, in_angle = direction - math.radians(bend), direction + math.pi + math.radians(bend)
    reach = 0.3915 * length
    first = (start[0] + reach * math.cos(out_angle), start[1] + reach * math.sin(out_angle))
    second = (end[0] + reach * math.cos(in_angle), end[1] + reach * math.sin(in_angle))
    return ((start[0] + 3 * first[0] + 3 * second[0] + end[0]) / 8, (start[1] + 3 * first[1] + 3 * second[1] + end[1]) / 8)


def circular_building(relative: str, label: str) -> tuple[dict[str, object], int, int]:
    """Transcribe the circular building figure with the given label from its TikZ source."""
    path = ARXIV / relative
    lines = path.read_text(encoding="utf-8").splitlines()
    last = next(i for i, line in enumerate(lines) if rf"\label{{{label}}}" in line)
    first = max(i for i in range(last) if r"\begin{tikzpicture}" in lines[i])
    body = "\n".join(lines[first : last + 1])
    angle = r"\(([\d*]+)(?:\*360/12)?:([\d.]+)cm\)"
    nodes: dict[tuple[float, float], dict[str, object]] = {}
    for fill, step, radius in re.findall(r"\\draw\[black, fill=(black|white)\] " + angle + r" circle", body):
        key = (int(step.split("*")[0]) * 30.0, float(radius))
        nodes[key] = {"kind": "curve" if fill == "black" else "point", "angle": key[0], "radius": key[1]}
    for step, radius, text in re.findall(r"\\node at " + angle + r" \{(\d+)\}", body):
        degrees, label_radius = int(step.split("*")[0]) * 30.0, float(radius)
        owner = min(nodes, key=lambda key: math.dist(polar(*key), polar(degrees, label_radius)))
        assert math.dist(polar(*owner), polar(degrees, label_radius)) < 0.31, f"{path}: label {text} at {degrees} is not next to a node"
        assert "label" not in nodes[owner], f"{path}: node at {owner} has two labels"
        nodes[owner]["label"] = int(text)
    assert all("label" in node for node in nodes.values()), f"{path}: unlabelled node"
    outer = sorted(key for key in nodes if key[1] == 4.0)
    inner = [key for key in nodes if key[1] != 4.0]
    edges = {frozenset((outer[k], outer[(k + 1) % len(outer)])) for k in range(len(outer))}
    for start, end in re.findall(r"\\path \((\d+)\*360/12:4cm\) edge \[bend right=30\] \((\d+)\*360/12:4cm\)", body):
        a, b = (int(start) * 30.0, 4.0), (int(end) * 30.0, 4.0)
        midpoint = bend_right_midpoint(polar(*a), polar(*b), 30)
        hit = [key for key in inner if math.dist(midpoint, polar(*key)) < 0.01]
        assert len(hit) == 1, f"{path}: chord {start}-{end} passes through {hit}"
        edges |= {frozenset((a, hit[0])), frozenset((b, hit[0]))}
    for edge in edges:
        assert {nodes[key]["kind"] for key in edge} == {"point", "curve"}, f"{path}: edge {edge} is not a point-curve incidence"

    building = {
        "points": sorted(node["label"] for node in nodes.values() if node["kind"] == "point"),
        "curves": sorted(node["label"] for node in nodes.values() if node["kind"] == "curve"),
        "incidences": sorted(
            [next(nodes[k]["label"] for k in edge if nodes[k]["kind"] == "point"), next(nodes[k]["label"] for k in edge if nodes[k]["kind"] == "curve")] for edge in edges
        ),
    }
    return building, first + 1, last + 1


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
                "definition": cite(
                    "2108.06236/main.tex",
                    204,
                    212,
                    r"L_{2d} = 2U \op \la -2d \ra \op \la -6 \ra",
                    r"\Gamma_{2d} = \{ g \in \opn{O}^+(L) \mid g \underline{v}^* \equiv \underline{v}^* \bmod{ L} \}",
                ),
                "theorem": cite(
                    "2108.06236/main.tex",
                    1044,
                    1047,
                    r"\label{L2boundarythm}",
                    "curves $\\cc_1$ and $\\cc_2$",
                    "points $P_1$, $P_2$, $P_3$",
                    r"$\overline{\cc}_1 \cap P_1$, $\overline{\cc}_1 \cap P_2$, $\overline{\cc}_2 \cap P_2$ and $\overline{\cc}_2 \cap P_3$",
                ),
            },
        },
    ]
    building, first, last = circular_building(paper, "2u2a2building")
    assert (len(building["points"]), len(building["curves"]), len(building["incidences"])) == (6, 9, 18), building
    assert building["points"] == list(range(6)) and building["curves"] == list(range(9))
    assert all(sum(p == point for p, _ in building["incidences"]) == 3 for point in building["points"])
    assert all(sum(c == curve for _, c in building["incidences"]) == 2 for curve in building["curves"])
    records.append(
        {
            "id": "dawes_2U2_A2_stable",
            "lattice": "2U(2) + A2",
            "gram": block_sum([[0, 2], [2, 0]], [[0, 2], [2, 0]], a2),
            "signature": [2, 4],
            "group": "stable orthogonal group (discriminant kernel of O^+)",
            "building": {"points": 6, "curves": 9, "edges": 18},
            "source": {
                "convention": convention,
                "statement": cite(paper, 1496, 1498, r"The building $\bc(\otl)$ is given in Figure \ref{2u2a2building}"),
                "figure": cite(paper, first, last, r"\caption{$\bc(\widetilde{\opn{O}}^+(2U(2) \op A_2))$}"),
            },
        }
    )
    incidence_graph = {
        "building_id": "dawes_2U2_A2_stable",
        "points": building["points"],
        "curves": building["curves"],
        "point_curve_incidences": building["incidences"],
        "source": {"kind": "published_figure", "file": str(ARXIV / paper), "lines": f"{first}-{last}", "label": "2u2a2building"},
    }
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
    TARGET.write_text(json.dumps({"buildings": records, "incidence_graphs": [incidence_graph], "index_chains": [index_chain]}, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {TARGET}")


if __name__ == "__main__":
    main()
