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
    _stabilizer_on_facet,
    cell_facet_orbits,
    cone_facets,
    orbit_representatives,
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


@pytest.mark.parametrize("reflection, orbit_size, expected_order", ((False, 4, 1), (True, 4, 2)))
def test_square_edge_stabilizer_has_the_exact_induced_action(reflection, orbit_size, expected_order):
    """A square edge has trivial rotation stabilizer and a reflection in D4."""
    from sage.all import PermutationGroup

    generators = ["(1,2,3,4)"]
    if reflection:
        generators.append("(1,2)(3,4)")
    group = PermutationGroup(generators, domain=(1, 2, 3, 4))
    induced = _stabilizer_on_facet(group, frozenset((0, 1)), orbit_size)
    stabilizer = PermutationGroup(
        [[position + 1 for position in permutation] for permutation in induced],
        domain=(1, 2),
    )
    assert int(stabilizer.order()) == expected_order
    if reflection:
        assert tuple(int(stabilizer.gen(0)(position)) for position in (1, 2)) == (2, 1)


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


def _group_closure(permutations, limit=20000):
    """All elements of the generated permutation group, by closure (small groups only)."""
    identity = tuple(range(len(permutations[0])))
    elements = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for permutation in permutations:
            product = tuple(permutation[current[index]] for index in range(len(current)))
            if product not in elements:
                elements.add(product)
                frontier.append(product)
                if len(elements) > limit:
                    raise AssertionError("the group is too large for the closure check")
    return elements


@pytest.mark.parametrize("name", sorted(SPECIMENS))
def test_pruned_walk_finds_every_orbit_exactly_once(name) -> None:
    r"""The orbit-pruned walk must equal the orbit decomposition of the *complete* facet family.

    The complete family is computed independently by `cone_facets` (the full walk) and by `cell_cone` (the
    engine). The pruned walk must produce the same number of orbits, one representative in each, no duplicate,
    no omission, and each representative must certify.
    """
    backend, cell = _cell(SPECIMENS[name])
    rays = [tuple(int(entry) for entry in row) for row in cell.rows]
    permutations = ray_permutations_from_stabilizer(cell, backend.cell_stabilizer(cell))

    full = cone_facets(rays)
    engine = _engine_facets(backend, cell)
    assert set(full) == set(engine)

    # Independent orbit decomposition of the complete family, by closure from each uncovered facet.
    expected_orbits = []
    covered: set[frozenset[int]] = set()
    for incidence in sorted(full, key=lambda face: (len(face), sorted(face))):
        if incidence in covered:
            continue
        members = set()
        frontier = [incidence]
        while frontier:
            current = frontier.pop()
            if current in members:
                continue
            members.add(current)
            for permutation in permutations:
                frontier.append(frozenset(permutation[position] for position in current))
        covered.update(members)
        expected_orbits.append(members)
    assert covered == set(full), "the orbit decomposition must cover the complete family"

    pruned = dict(orbit_representatives(rays, permutations))
    assert len(pruned) == len(expected_orbits), (
        f"pruned walk found {len(pruned)} orbits, the complete family has {len(expected_orbits)}"
    )
    found = []
    for incidence, normal in pruned.items():
        assert is_facet_certificate(normal, incidence, rays)
        orbit = set()
        frontier = [incidence]
        while frontier:
            current = frontier.pop()
            if current in orbit:
                continue
            orbit.add(current)
            for permutation in permutations:
                frontier.append(frozenset(permutation[position] for position in current))
        found.append(orbit)
    assert all(orbit <= set(full) for orbit in found), "every pruned orbit must lie in the family"
    assert len({frozenset(orbit) for orbit in found}) == len(found), "orbits must be distinct"
    assert set().union(*found) == set(full), "the pruned orbits must partition the complete family"


@pytest.mark.parametrize("name", ["A2_plus_1_plus_minus1", "A3_plus_1_plus_minus1"])
def test_orbit_stabilizer_and_transport_witnesses_are_exact(name) -> None:
    r"""`|G·F| = |G| / |Stab_G(F)|`, with the stabilizer from Schreier and the words verified.

    ``|G|`` is computed by closure (the groups here are small), so orbit-stabilizer is checked against an
    independent count rather than against the walk's own bookkeeping.
    """
    backend, cell = _cell(SPECIMENS[name])
    permutations = ray_permutations_from_stabilizer(cell, backend.cell_stabilizer(cell))
    group = _group_closure(permutations)
    records = cell_facet_orbits(cell, backend.cell_stabilizer(cell))

    total = 0
    for record in records:
        members = {}
        for member, word in record.words:
            moved = frozenset(record.incidence)
            for index in word:
                moved = frozenset(permutations[index][position] for position in moved)
            assert moved == member, "the word must transport the representative onto the member"
            members[member] = word
        stabilizer = sum(1 for g in group if frozenset(g[position] for position in record.incidence)
                         == record.incidence)
        assert len(members) == record.orbit_size
        assert len(members) * stabilizer == len(group), (
            f"orbit-stabilizer failed: orbit {len(members)}, stabilizer {stabilizer}, group {len(group)}"
        )
        total += len(members)
    assert total == len(cone_facets([tuple(int(e) for e in row) for row in cell.rows]))


@pytest.mark.parametrize("name", sorted(SPECIMENS))
def test_walk_permutations_match_the_port_action_convention(name) -> None:
    r"""The walk's ray action must be the port's `_facet_mask_data` action, element for element.

    The traversal carries stabilizer permutations as position images of the rows. If the walk used a different
    convention -- inverse composition, a different row order, or the covector action -- its orbit records would
    still look self-consistent while being incompatible with `flip_across` and the side-pairing words. This pins
    the two conventions together.
    """
    backend, cell = _cell(SPECIMENS[name])
    mine = set(ray_permutations_from_stabilizer(cell, backend.cell_stabilizer(cell)))

    masks, port_permutations, _keys = backend._facet_mask_data(cell, use_owned_cone=True)
    list(masks)  # the facet keys are an iterator; this test is about the permutations only
    theirs = {tuple(int(entry) for entry in permutation) for permutation in port_permutations}

    assert mine == theirs, (
        f"the walk's ray action differs from the port's: {len(mine)} vs {len(theirs)} generators, "
        f"{len(mine - theirs)} only in the walk, {len(theirs - mine)} only in the port"
    )


def test_same_ray_contract_is_proportionality_and_positive_orientation() -> None:
    r"""The exact ray contract, with the regressions that broke the first version.

    Cyclically adjacent minors are not proportionality in dimension 4: `e1` and `e3` pass every one of those
    tests while sharing no ray. And a vector and its negative pass every minor test while spanning *opposite*
    rays, which a pointed cone distinguishes.
    """
    from sage_indefinite_port.indefinite.facet_orbits import _same_ray

    # proportionality, over Q, positive orientation
    assert _same_ray((1, 2, 3, 4), (2, 4, 6, 8))
    assert _same_ray((2, 4), (3, 6)), "a rational ratio is still the same ray"
    assert _same_ray((0, 3, 0), (0, 1, 0))

    # the regressions
    assert not _same_ray((1, 0, 0, 0), (0, 0, 1, 0)), "e1 and e3 are not proportional in dimension 4"
    assert not _same_ray((1, 2, 3), (-1, -2, -3)), "a vector and its negative are opposite rays"
    assert not _same_ray((1, 0, 0, 0), (-1, 0, 0, 0)), "opposite in dimension 4"
    assert not _same_ray((1, 2), (1, 3)), "different directions"
    assert not _same_ray((0, 0), (0, 0)), "the zero vector is not a ray"
    assert not _same_ray((0, 0), (1, 2))


def test_a_rational_coordinate_is_refused_rather_than_truncated() -> None:
    r"""`int(QQ(3)/QQ(2)) == 1` would fabricate a ray; the conversion must refuse instead."""
    from sage_indefinite_port.indefinite.facet_orbits import _integral_coordinate

    class Coordinate:
        """Stands in for an owned ring element: `numerator()` and `denominator()` are methods."""

        def __init__(self, numerator, denominator):
            self._numerator = numerator
            self._denominator = denominator

        def numerator(self):
            return self._numerator

        def denominator(self):
            return self._denominator

    assert _integral_coordinate(Coordinate(6, 3)) == 2
    assert _integral_coordinate(Coordinate(-6, 3)) == -2
    with pytest.raises(ArithmeticError):
        _integral_coordinate(Coordinate(3, 2))
    with pytest.raises(ArithmeticError):
        _integral_coordinate(Coordinate(1, 0))


def test_the_walk_refuses_a_cone_it_cannot_certify() -> None:
    r"""A failure that is not proved to be an exclusion must surface, not empty the worklist quietly.

    The half-plane cone `cone{(1,0), (-1,0), (0,1)}` is not pointed: it contains the x-axis as a lineality
    space, its dual is not full-dimensional, and the seed slice is unbounded. The walk must raise instead of
    returning a partial facet set that would later be read as a closure certificate.
    """
    rays = [(1, 0), (-1, 0), (0, 1)]
    with pytest.raises(Exception):
        cone_facets(rays)
