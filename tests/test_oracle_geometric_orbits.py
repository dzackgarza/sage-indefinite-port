"""Independent geometric building oracles from Dawes, arXiv:2205.10601.

Pinned primary source: references/vendor/arxiv/2205.10601/orbits_in_lattices.tex.
The building is the bipartite incidence graph of group-orbits of rational
isotropic lines and planes (definitions at source lines 800-812); its edges
are distinct incident pairs of line and plane orbits (lines 803-809),
not necessarily the set of orbits of complete flags.

Each expected triple below is counted directly from the indicated published
TikZ figure, not from the port, its existing fixture JSON, or an old test.
The stable group is O^+(L) intersect ker[O(L)->O(D(L))], as in the paper.
"""

import pytest

from tests.acceptance.consumer import (
    MARKED_CLASSES,
    lattice,
    require,
)

_U = ((0, 1), (1, 0))
_A2_NEG = ((-2, 1), (1, -2))


def _block_diagonal(*blocks):
    size = sum(len(block) for block in blocks)
    result = [[0] * size for _ in range(size)]
    start = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, entry in enumerate(row):
                result[start + i][start + j] = entry
        start += len(block)
    return result


def _stable_positive_group(L):
    return L.stable_orthogonal_group().intersection(L.O_plus())


@pytest.mark.parametrize(
    "gram,stable,expected",
    [
        pytest.param(
            _block_diagonal(_U, _U, _A2_NEG),
            True,
            (1, 1, 1),
            id="dawes-stable-2U-A2",
        ),
        pytest.param(
            _block_diagonal(_U, _U, ((-2,),), ((-6,),)),
            False,
            (2, 2, 3),
            id="dawes-Oplus-2U-minus2-minus6",
        ),
        pytest.param(
            _block_diagonal(((0, 2), (2, 0)), ((0, 2), (2, 0)), _A2_NEG),
            True,
            (9, 6, 18),
            id="dawes-stable-2U2-A2",
        ),
    ],
)
def test_published_tits_building_orbit_counts(gram, stable, expected):
    """Published figure counts (points, curves, incidences).

    Source correspondence:
      * 2U+A2: Example at 1310-1347, Fig. at 1349-1364: 1 white,
        1 black, 1 connecting edge.
      * 2U+<-2>+<-6>: Example at 1391-1424, Fig. at 1426-1444:
        2 white, 2 black, 3 consecutive edges; group explicitly O^+(L).
      * 2U(2)+A2: Example at 1497-1500, Fig. at 1450-1495:
        6 white, 9 black (6 outer + 3 inner), 18 edges
        (12 perimeter segments + 3 curved chords, each split into two
        distinct incidences by the black vertex lying on that chord).

    These are distinct arithmetic subgroups; do not substitute O(L).
    """
    require(MARKED_CLASSES)
    L = lattice(gram)
    group = _stable_positive_group(L) if stable else L.O_plus()
    lines, planes = tuple(group.cusps(1)), tuple(group.cusps(2))
    incidences = {(lines.index(edge.line_cusp()), planes.index(edge.plane_cusp()))
                  for edge in group.tits_building_incidence()}
    computed = (len(lines), len(planes), len(incidences))
    assert computed == expected
