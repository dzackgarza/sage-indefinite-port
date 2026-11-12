"""Diagnostic witness for the dawes_2U2_A2_stable incidence graph (#64), on U(2) ⊕ U(2) ⊕ A2.

The frozen test reports only that the canonical graphs differ. This witness prints the
computed records, the distinct (line class, plane class) pairs and the degree sequences, then
asserts the published structure: 18 distinct pairs, every point on 3 curves and every curve
through 2 points (Dawes, arXiv:2205.10601, Figure 2u2a2building).
"""
from collections import Counter
import json

from tests.acceptance.consumer import lattice


def test_u2u2a2_stable_building_records():
    case = next(c for c in json.load(open("tests/fixtures/dawes_buildings.json"))["buildings"] if c["id"] == "dawes_2U2_A2_stable")
    ambient = lattice(case["gram"])
    group = ambient.stable_orthogonal_group().intersection(ambient.O_plus())
    lines, planes = list(group.cusps(1)), list(group.cusps(2))
    pairs = [(lines.index(r.line_cusp()), planes.index(r.plane_cusp())) for r in group.tits_building_incidence()]
    print("lines", len(lines), "planes", len(planes), "records", len(pairs), flush=True)
    print("pairs", sorted(pairs), flush=True)
    print("point degrees", sorted(Counter(p for p, _c in set(pairs)).values()), flush=True)
    print("curve degrees", sorted(Counter(c for _p, c in set(pairs)).values()), flush=True)
    assert len(set(pairs)) == 18
    assert sorted(Counter(p for p, _c in set(pairs)).values()) == [3] * 6
    assert sorted(Counter(c for _p, c in set(pairs)).values()) == [2] * 9
