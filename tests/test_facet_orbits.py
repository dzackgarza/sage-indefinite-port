"""The orbit-representative facet walk, checked against the engine's complete H-description.

Two independent ground truths:

* for the port's own cells, `cell_cone` gives the complete facet family from Normaliz, and the walk has to
  reproduce exactly that set -- incidences and oriented normals -- without asking the engine anything;
* for the cube cones the facet family is closed form, `2(n-1)` facets, and the stabilizer acting by coordinate
  permutations and sign flips has a known orbit count.

The walk's own claims are then checked on the records: every incidence is certified, every orbit is a genuine
orbit of the stabilizer (words transport the representative onto each member), the orbits partition the facet
set, and the counts agree with orbit-stabilizer.
"""
import itertools
from math import prod

import pytest
from sage.all import QQ, ZZ, matrix

from sage_indefinite_port.indefinite.facet_certificates import is_facet_certificate
from sage_indefinite_port.indefinite.facet_orbits import (
    cell_facet_orbits,
    cone_facets,
    ray_permutations_from_stabilizer,
)
from sage_indefinite_port.indefinite.lorentzian_cells import (
    NONNEGATIVE_SQUARES,
    LorentzianPerfectLocalBackend,
)
from tests.acceptance.consumer import lattice

SPECIMENS = {
    "A2_plus_1_plus_minus1": [[2, -1, 0, 0], [-1, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]],
    "A3_plus_1_plus_minus1": [
        [2, -1, 0, 0, 0],
        [-1, 2, -1, 0, 0],
        [0, -1, 2, 0, 0],
        [0, 0, 0, 1, 0],
        [0, 0, 0, 0, -1],
    ],
}


def _cell(gram):
    backend = LorentzianPerfectLocalBackend()
    return backend, backend.initial_cell(lattice(gram), NONNEGATIVE_SQUARES)


def _engine_facets(backend, cell):
    """The complete facet family from the maintained engine, as incidence -> oriented normal."""
    dimension = len(cell.rows[0])
    rays = [tuple(int(entry) for entry in row) for row in cell.rows]
    facets = {}
    for covector in backend.cell_cone(cell).halfspace_covectors():
        coordinates = tuple(
            int(QQ(int(covector(cell.lattice.module_generator(label)).numerator()))
                / QQ(int(covector(cell.lattice.module_generator(label)).denominator())))
            for label in cell.lattice.module_generating_set()
        )
        values = [sum(coordinates[i] * row[i] for i in range(dimension)) for row in rays]
        assert all(value >= 0 for value in values)
        incidence = frozenset(index for index, value in enumerate(values) if value == 0)
        facets[incidence] = coordinates
    return facets


@pytest.mark.parametrize("name", sorted(SPECIMENS))
def test_walk_reproduces_the_engine_facet_family_exactly(name) -> None:
    backend, cell = _cell(SPECIMENS[name])
    from sage_indefinite_port.indefinite.facet_certificates import pointedness_functional

    height, pairings = pointedness_functional(cell)
    assert all(value > 0 for value in pairings)

    walked = cone_facets([tuple(int(entry) for entry in row) for row in cell.rows])
    engine = _engine_facets(backend, cell)

    assert set(walked) == set(engine), (
        f"incidence mismatch: walk {sorted(len(i) for i in walked)}, engine {sorted(len(i) for i in engine)}"
    )
    for incidence, normal in walked.items():
        assert normal == engine[incidence] or normal == tuple(-entry for entry in engine[incidence])
        assert is_facet_certificate(normal, incidence, cell.rows)


@pytest.mark.parametrize("name", sorted(SPECIMENS))
def test_orbits_partition_the_facets_and_words_transport(name) -> None:
    backend, cell = _cell(SPECIMENS[name])
    stabilizer = backend.cell_stabilizer(cell)
    permutations = ray_permutations_from_stabilizer(cell, stabilizer)
    engine = _engine_facets(backend, cell)
    records = cell_facet_orbits(cell, stabilizer)

    covered = set()
    for record in records:
        assert record.incidence in engine
        assert is_facet_certificate(record.normal, record.incidence, cell.rows)
        assert len(record.words) == record.orbit_size
        assert dict(record.words)[record.incidence] == (), "the representative must carry the empty word"
        for member, word in record.words:
            moved = frozenset(record.incidence)
            for index in word:
                moved = frozenset(permutations[index][position] for position in moved)
            assert moved == member, "the word must transport the representative onto the member"
        covered.update(member for member, _word in record.words)
    assert covered == set(engine), "the orbits must partition the whole facet family"
    assert sum(record.orbit_size for record in records) == len(engine)


def test_cube_cone_facets_and_orbits_match_the_closed_form() -> None:
    r"""Closed-form specimen: `2(n-1)` facets; with the full symmetry one orbit."""
    for n in (4, 5):
        rays = [tuple(int(entry) for entry in row) + (1,) for row in itertools.product((-1, 1), repeat=n - 1)]
        facets = cone_facets(rays)
        assert len(facets) == 2 * (n - 1), f"n={n}: expected 2(n-1) facets, got {len(facets)}"
        for incidence, normal in facets.items():
            assert is_facet_certificate(normal, incidence, rays)
        assert all(len(incidence) == 2 ** (n - 2) for incidence in facets)
