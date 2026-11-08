"""The traversal seam: the walk's facet orbits must be the ones the engine route produces.

`_expand_traversal_once` takes its facet orbits from the pinned split/direct dispatch. An independent engine
route computes the complete H-description and groups it.
What the traversal consumes is identical either way: per orbit the representative's `FacetMaskKey`, and -- when
side-pairing witnesses are retained -- a `(key, Schreier word)` entry for every member, the word carrying the
representative onto that member.

These tests compare the two routes on the bounded specimens, where both can run, on exactly those objects:
representative incidences, orbit membership, the neighbour obtained by `flip_across` on each representative, and
the action of every word. They also compare complete traversals through the source dispatch and the independent
engine route.
"""
import pytest
from sage.all import QQ, ZZ, matrix

from sage_indefinite_port.indefinite import lorentzian_cells
from sage_indefinite_port.indefinite.lorentzian_cells import (
    NONNEGATIVE_SQUARES,
    LorentzianPerfectLocalBackend,
    _facet_incidence,
    _schreier_generator_indices,
)
from tests.acceptance.consumer import lattice

SPECIMENS = {
    "A2_plus_1_plus_minus1": [[2, -1, 0, 0], [-1, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]],
    "A3_plus_1_plus_minus1": [
        [2, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 2, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, -1]
    ],
}


def _cell(gram):
    backend = LorentzianPerfectLocalBackend()
    return backend, backend.initial_cell(lattice(gram), NONNEGATIVE_SQUARES)


def _engine_orbits(backend, cell):
    facet_masks, permutations, canonical = backend._facet_mask_data(cell, use_owned_cone=False)
    return permutations, list(
        backend._iter_facet_mask_orbits(facet_masks, permutations, retain_words=True, canonical_keys=canonical)
    )


def _apply(word, permutations, incidence):
    moved = frozenset(incidence)
    for generator_index in _schreier_generator_indices(word):
        moved = frozenset(permutations[generator_index][position] for position in moved)
    return moved


@pytest.mark.parametrize("name", sorted(SPECIMENS))
def test_walk_orbits_match_the_engine_route_at_the_seam(name) -> None:
    backend, cell = _cell(SPECIMENS[name])
    ray_count = len(cell.rows)

    engine_permutations, engine_orbits = _engine_orbits(backend, cell)
    walk_permutations, walk_orbits = backend._facet_mask_orbits_via_walk(cell, retain_words=True)

    assert {tuple(int(entry) for entry in permutation) for permutation in engine_permutations} == {
        tuple(int(entry) for entry in permutation) for permutation in walk_permutations
    }, "the two routes must use the same ray permutations"

    def incidence_sets(orbits):
        groups = []
        for orbit in orbits:
            groups.append(frozenset(_facet_incidence(key, ray_count) for key, _word in orbit))
        return set(groups)

    assert incidence_sets(engine_orbits) == incidence_sets(walk_orbits), (
        "the two routes must produce the same orbits as sets of facet incidences"
    )


@pytest.mark.parametrize("name", sorted(SPECIMENS))
def test_walk_words_transport_and_neighbours_agree(name) -> None:
    backend, cell = _cell(SPECIMENS[name])
    ray_count = len(cell.rows)
    engine_permutations, engine_orbits = _engine_orbits(backend, cell)
    walk_permutations, walk_orbits = backend._facet_mask_orbits_via_walk(cell, retain_words=True)

    # Every word must carry the representative onto the member it is paired with, in both encodings.
    for permutations, orbits in ((engine_permutations, engine_orbits), (walk_permutations, walk_orbits)):
        for orbit in orbits:
            representative = _facet_incidence(orbit[0][0], ray_count)
            for key, word in orbit:
                member = _facet_incidence(key, ray_count)
                assert _apply(word, permutations, representative) == member, (
                    "the Schreier word must transport the representative onto its member"
                )

    # The neighbour across each representative, and the cell it produces, must agree between routes.
    def neighbours(orbits):
        produced = []
        for orbit in orbits:
            representative = _facet_incidence(orbit[0][0], ray_count)
            neighbour = backend.flip_across(cell, representative)
            if neighbour is None:
                produced.append(None)
            else:
                produced.append(tuple(sorted(tuple(int(entry) for entry in row) for row in neighbour.rows)))
        return sorted(produced, key=lambda item: (item is not None, item))

    assert neighbours(engine_orbits) == neighbours(walk_orbits), (
        "flipping across the representative must give the same neighbouring cell on both routes"
    )


@pytest.mark.parametrize("name", sorted(SPECIMENS))
def test_source_dispatch_traverses_to_the_same_generators_as_engine(name, monkeypatch) -> None:
    r"""Compare the source dispatch's complete traversal with direct facet enumeration.

    `perfect_complex` is cached per lattice, so the cache is cleared between the two runs; otherwise the second
    call would return the first route's complex and the comparison would be vacuous.
    """
    from sage_indefinite_port.indefinite.lorentzian_cells import perfect_complex

    gram = lattice(SPECIMENS[name])

    def generators():
        complex_ = perfect_complex(gram, NONNEGATIVE_SQUARES, retain_facet_orbits=True)
        return tuple(
            tuple(tuple(int(entry) for entry in row) for row in _rows_of(generator))
            for generator in complex_.component_preserving_generators()
        )

    original = LorentzianPerfectLocalBackend._facet_mask_orbits_via_walk

    def direct(self, cell, *, retain_words=True):
        masks, permutations, canonical = self._facet_mask_data(cell, use_owned_cone=False)
        return permutations, self._iter_facet_mask_orbits(
            masks, permutations, retain_words=retain_words, canonical_keys=canonical
        )

    lorentzian_cells._cached_perfect_complex.cache_clear()
    monkeypatch.setattr(LorentzianPerfectLocalBackend, "_facet_mask_orbits_via_walk", direct)
    engine = sorted(generators())

    lorentzian_cells._cached_perfect_complex.cache_clear()
    monkeypatch.setattr(LorentzianPerfectLocalBackend, "_facet_mask_orbits_via_walk", original)
    walked = sorted(generators())

    assert engine == walked, (
        f"the source dispatch produced {len(walked)} generators, the engine route {len(engine)}, "
        "and they must agree as exact row matrices"
    )


def _rows_of(isometry):
    """The row matrix of an owned isometry on its lattice's declared generators."""
    lattice_ = isometry.domain()
    labels = tuple(lattice_.module_generating_set())
    return [
        [int(QQ(int(isometry(lattice_.module_generator(label)).to_vector()(target).numerator()))
             / QQ(int(isometry(lattice_.module_generator(label)).to_vector()(target).denominator())))
         for target in labels]
        for label in labels
    ]
