"""Integral transport of perfect cells in U + E8(-2).

The explicit unimodular matrix certifies that the two complete configurations
are isometric, independently of the transporter search.
"""

from array import array

import pytest
from sage.all import ZZ, identity_matrix, matrix, vector

from sage_indefinite_port.indefinite.lorentzian_cells import (
    _CellFaces,
    _bit_positions,
    _facet_complement_key,
    _facet_incidence,
    _iter_owned_row_blocks,
    _permuted_bitset,
    _permuted_facet_key,
    _perfect_form_hash_key,
    _schreier_generator_indices,
    LorentzianPerfectCell,
    LorentzianPerfectLocalBackend,
    NONNEGATIVE_SQUARES,
    ZERO_SQUARES,
    _extend_configuration_signatures,
)
import sage_indefinite_port.indefinite.lorentzian_cells as lorentzian_cells
from tests.acceptance.consumer import lattice, rows_of


def test_owned_row_blocks_release_only_consumed_list_entries() -> None:
    rows = [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]]
    blocks = _iter_owned_row_blocks(rows, 2)

    first = next(blocks)
    assert first == ((1, 2), (3, 4))
    assert rows == [None, None, [5, 6], [7, 8], [9, 10]]
    second = next(blocks)
    assert second == ((5, 6), (7, 8))
    assert rows == [None, None, None, None, [9, 10]]
    assert next(blocks) == ((9, 10),)
    assert rows == [None] * 5
    with pytest.raises(StopIteration):
        next(blocks)


def test_owned_row_blocks_leave_immutable_sequences_unchanged() -> None:
    rows = ((1, 2), (3, 4), (5, 6))
    assert tuple(_iter_owned_row_blocks(rows, 2)) == (
        ((1, 2), (3, 4)),
        ((5, 6),),
    )
    assert rows == ((1, 2), (3, 4), (5, 6))

    with pytest.raises(ValueError, match="block size must be positive"):
        tuple(_iter_owned_row_blocks(rows, 0))


def test_schreier_generator_indices_streams_nested_word_in_action_order() -> None:
    word = (((None, 3), 1), 4)
    indices = _schreier_generator_indices(word)

    assert iter(indices) is indices
    assert tuple(indices) == (3, 1, 4)


def test_facet_orbit_materializers_release_consumed_bfs_entries(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    backend = LorentzianPerfectLocalBackend()
    original = backend._iter_facet_mask_orbits
    emitted = []

    def recording_orbits(facet_keys, permutations):
        for orbit in original(facet_keys, permutations):
            emitted.append(orbit)
            yield orbit

    monkeypatch.setattr(backend, "_iter_facet_mask_orbits", recording_orbits)
    permutation_orbits = backend.facet_orbit_permutations(cell)
    assert permutation_orbits
    assert emitted and all(not any(orbit) for orbit in emitted)

    emitted.clear()
    isometry_orbits = backend.facet_orbits(cell)
    assert isometry_orbits
    assert emitted and all(not any(orbit) for orbit in emitted)

GRAM = ((-4, 0, 0, 0, 0, 2, 0, 0, 0, 0),
 (0, -4, 2, 0, 0, 0, 0, 0, 0, 0),
 (0, 2, -4, 2, 0, 0, 0, 0, 0, 0),
 (0, 0, 2, -4, 2, 0, 0, 0, 0, 0),
 (0, 0, 0, 2, -4, 2, 0, 0, 0, 0),
 (2, 0, 0, 0, 2, -4, 2, 0, 0, 0),
 (0, 0, 0, 0, 0, 2, -4, 2, 0, 0),
 (0, 0, 0, 0, 0, 0, 2, -4, 0, 0),
 (0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
 (0, 0, 0, 0, 0, 0, 0, 0, 1, 0))

SOURCE_ROWS = ((-8, -3, -6, -9, -12, -15, -10, -5, 2, 4),
 (0, 0, 0, 0, 0, 0, 0, 0, 0, 1),
 (-5, -1, -3, -5, -7, -9, -6, -3, 2, 2),
 (-5, -2, -3, -5, -7, -9, -6, -3, 2, 2),
 (-5, -2, -4, -5, -7, -9, -6, -3, 2, 2),
 (-5, -2, -4, -6, -7, -9, -6, -3, 2, 2),
 (-5, -2, -4, -6, -8, -10, -6, -3, 2, 2),
 (-5, -2, -4, -6, -8, -9, -6, -3, 2, 2),
 (-5, -2, -4, -6, -8, -10, -7, -3, 2, 2),
 (-5, -2, -4, -6, -8, -10, -7, -4, 2, 2))

TARGET_ROWS = ((-5, -2, -4, -6, -8, -10, -7, -3, 2, 2),
 (-11, -4, -8, -12, -16, -20, -14, -7, 4, 4),
 (-11, -5, -9, -13, -17, -21, -14, -7, 4, 4),
 (-5, -2, -4, -6, -8, -9, -6, -3, 2, 2),
 (-5, -2, -3, -5, -7, -9, -6, -3, 2, 2),
 (-5, -2, -4, -5, -7, -9, -6, -3, 2, 2),
 (-5, -2, -4, -6, -7, -9, -6, -3, 2, 2),
 (-11, -4, -8, -12, -16, -20, -13, -6, 4, 4),
 (-8, -3, -6, -9, -12, -15, -10, -5, 4, 2),
 (-9, -4, -7, -10, -13, -16, -11, -5, 4, 3))

CERTIFICATE = ((-8, -3, -6, -8, -11, -14, -9, -4, 4, 2),
 (0, 1, 1, 1, 1, 1, 0, 0, 0, 0),
 (-6, -3, -5, -7, -9, -12, -8, -4, 2, 2),
 (0, 0, -1, -1, -1, 0, 0, 0, 0, 0),
 (0, 0, 1, 1, 0, 0, 0, 0, 0, 0),
 (0, 0, 0, -1, 0, 0, 0, 0, 0, 0),
 (6, 2, 4, 7, 9, 11, 7, 3, -2, -2),
 (-3, -1, -2, -3, -4, -5, -3, -1, 0, 2),
 (-12, -5, -9, -13, -17, -21, -14, -7, 5, 4),
 (-9, -4, -7, -10, -13, -16, -11, -5, 4, 3))


def test_integral_transporter_of_spanning_perfect_cells() -> None:
    gram = matrix(ZZ, GRAM)
    certificate = matrix(ZZ, CERTIFICATE)
    source_rows = matrix(ZZ, SOURCE_ROWS)
    target_rows = matrix(ZZ, TARGET_ROWS)
    assert abs(certificate.det()) == 1
    assert certificate * gram * certificate.transpose() == gram
    assert set(map(tuple, source_rows * certificate)) == set(TARGET_ROWS)

    ambient = lattice(GRAM)
    source = LorentzianPerfectCell(ambient, SOURCE_ROWS, NONNEGATIVE_SQUARES)
    target = LorentzianPerfectCell(ambient, TARGET_ROWS, NONNEGATIVE_SQUARES)
    witness = LorentzianPerfectLocalBackend().cell_transporter(source, target)

    assert witness is not None
    transform = matrix(ZZ, rows_of(witness))
    assert abs(transform.det()) == 1
    assert transform * gram * transform.transpose() == gram
    assert set(map(tuple, source_rows * transform)) == set(map(tuple, target_rows))


def test_all_integral_transporters_extend_exact_anchor_signatures() -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    transporters = LorentzianPerfectLocalBackend().cell_transporter(cell, cell, all_images=True)
    assert len(transporters) == 2
    for witness in transporters:
        transform = matrix(ZZ, rows_of(witness))
        assert abs(transform.det()) == 1
        assert transform * gram * transform.transpose() == gram
        assert set(map(tuple, matrix(ZZ, cell.rows) * transform)) == set(cell.rows)


def test_coarse_large_configuration_colours_preserve_exact_stabilizer(monkeypatch) -> None:
    """Skipping quadratic graph refinement must not weaken witness checks."""
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)
    monkeypatch.setattr(
        lorentzian_cells._ConfigurationPairings,
        "blocks",
        lambda _self: pytest.fail("coarse invariant path must not construct the quadratic pairing graph"),
    )

    transporters = LorentzianPerfectLocalBackend().cell_transporter(
        cell, cell, all_images=True, generators_only=True
    )

    assert len(transporters) == 1
    for witness in transporters:
        transform = matrix(ZZ, rows_of(witness))
        assert transform != identity_matrix(ZZ, 2)
        assert abs(transform.det()) == 1
        assert transform * gram * transform.transpose() == gram
        assert set(map(tuple, matrix(ZZ, cell.rows) * transform)) == set(cell.rows)


def test_coarse_colours_preserve_changed_configuration_transporters(monkeypatch) -> None:
    """Coarse invariant colours cannot reject a valid nonidentity transporter."""
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    ambient = lattice(gram.rows())
    source = LorentzianPerfectCell(ambient, ((1, 0), (0, 1)), ZERO_SQUARES)
    target = LorentzianPerfectCell(ambient, ((-1, 0), (0, -1)), ZERO_SQUARES)
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)
    monkeypatch.setattr(
        lorentzian_cells._ConfigurationPairings,
        "blocks",
        lambda _self: pytest.fail("coarse invariant path must not construct the quadratic pairing graph"),
    )
    monkeypatch.setattr(
        lorentzian_cells,
        "Sets",
        lambda: pytest.fail("large single-transporter search must not use the preamble pairing-graph transporter"),
    )
    assert _perfect_form_hash_key(source) == _perfect_form_hash_key(target)

    backend = LorentzianPerfectLocalBackend()
    witness = backend.cell_transporter(source, target)
    assert witness is not None
    witness_matrix = matrix(ZZ, rows_of(witness))
    assert abs(witness_matrix.det()) == 1
    assert witness_matrix * gram * witness_matrix.transpose() == gram
    assert set(map(tuple, matrix(ZZ, source.rows) * witness_matrix)) == set(target.rows)

    transporters = backend.cell_transporter(source, target, all_images=True)

    assert len(transporters) == 2
    for witness in transporters:
        transform = matrix(ZZ, rows_of(witness))
        assert abs(transform.det()) == 1
        assert transform * gram * transform.transpose() == gram
        assert set(map(tuple, matrix(ZZ, source.rows) * transform)) == set(target.rows)


def test_anchor_signature_extension_appends_in_common_assignment_order() -> None:
    signatures = ((0, (10, 30)), (1, (11, 31)), (0, (12, 32)))
    column = (20, 21, 22)
    assert _extend_configuration_signatures(signatures, column) == (
        (0, (10, 30, 20)),
        (1, (11, 31, 21)),
        (0, (12, 32, 22)),
    )


def test_empty_stable_configuration_has_empty_permutation_for_each_generator() -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    backend = LorentzianPerfectLocalBackend()
    generators = backend.cell_stabilizer(cell)

    assert backend.point_permutations_from_rows(cell, ()) == tuple(() for _ in generators)


def test_permuted_facet_bitsets_preserve_exact_positions() -> None:
    permutation = (2, 0, 3, 1, 5, 4, 6, 7, 9, 8)
    all_positions = (1 << len(permutation)) - 1
    masks = (0, 1, 1 << 9, (1 << 0) | (1 << 3) | (1 << 8) | (1 << 9), (1 << 10) - 1)
    for mask in masks:
        expected = sum(1 << permutation[position] for position in range(len(permutation)) if mask & (1 << position))
        image = _permuted_bitset(mask, permutation)
        assert image == expected
        assert _permuted_bitset(mask, array("H", permutation)) == expected
        # The production facet traversal stores nonincident-ray complements;
        # a permutation must commute with complement in the fixed ray universe.
        assert _permuted_bitset(all_positions ^ mask, permutation) == all_positions ^ image

    # Match the A7 configuration scale and the production sparse-complement
    # representation: only four of 49,861 ray positions are nonincident.
    ray_count = 49_861
    reverse = tuple(reversed(range(ray_count)))
    sparse_complement = sum(1 << position for position in (0, 1, ray_count // 2, ray_count - 1))
    sparse_key = _facet_complement_key(sparse_complement, ray_count)
    assert sparse_key == (True, (0, 1, ray_count // 2, ray_count - 1))
    sparse_image = _permuted_facet_key(sparse_key, reverse)
    assert sparse_image == (True, tuple(sorted(ray_count - 1 - position for position in sparse_key[1])))
    assert _facet_incidence(sparse_key, ray_count) == frozenset(range(ray_count)).difference(sparse_key[1])
    sparse_incident = _facet_complement_key(((1 << ray_count) - 1) ^ sparse_complement, ray_count)
    assert sparse_incident == (False, sparse_key[1])
    assert _facet_incidence(sparse_incident, ray_count) == frozenset(sparse_key[1])
    assert _permuted_facet_key(sparse_incident, reverse) == (
        False, tuple(sorted(ray_count - 1 - position for position in sparse_incident[1]))
    )
    # For small ray universes, an integer bitset remains compact.
    small_key = _facet_complement_key((1 << 1) | (1 << 3), 10)
    assert isinstance(small_key, int)

    dense_complement = sum(1 << position for position in range(200))
    dense_key = _facet_complement_key(dense_complement, ray_count)
    assert isinstance(dense_key, int)
    dense_image = _permuted_facet_key(dense_key, reverse)
    assert dense_image == sum(1 << (ray_count - 1 - position) for position in range(200))


def test_normaliz_pairing_rows_match_owned_covectors() -> None:
    """Large-cell incidence reuses the owned H-family without changing zeros."""
    gram = matrix(ZZ, [[0, 1], [1, 0]])

    class LargeInput:
        pass

    cell = LargeInput()
    cell.lattice = lattice(gram.rows())
    # Trigger the large-cell Normaliz path without an expensive V-to-H case:
    # Normaliz removes these repeated input rays exactly.
    cell.rows = ((1, 0), (0, 1)) * 5_000
    backend = LorentzianPerfectLocalBackend()
    labels = tuple(cell.lattice.module_generating_set())
    pairing_rows = tuple(backend._iter_cell_halfspace_pairing_rows(cell, labels))
    assert cell not in backend.__dict__.get("_normaliz_cell_cones", {})
    owned_rows = tuple(
        tuple(int(covector(label)) for label in labels)
        for covector in backend._iter_cell_halfspace_covectors(cell)
    )
    assert cell not in backend.__dict__.get("_normaliz_cell_cones", {})
    assert tuple(tuple(value == 0 for value in row) for row in pairing_rows) == tuple(
        tuple(value == 0 for value in row) for row in owned_rows
    )


def test_large_perfect_form_hash_skips_quadratic_pairing_graph(monkeypatch) -> None:
    """Large-cell bucketing uses exact coarse invariants, not all pairings."""
    gram = matrix(ZZ, [[0, 1], [1, 0]])

    class LargeInput:
        pass

    cell = LargeInput()
    cell.lattice = lattice(gram.rows())
    cell.rows = ((1, 0), (0, 1)) * 5_000
    monkeypatch.setattr(
        lorentzian_cells._ConfigurationPairings,
        "blocks",
        lambda _self: pytest.fail("large perfect-form hashing must not build the quadratic pairing graph"),
    )

    key = _perfect_form_hash_key(cell)

    assert key[1] == ()
    assert sum(count for _value, count in key[0]) == len(cell.rows)
    assert len(key[2]) == 2
    assert sum(count for _content, count in key[3]) == len(cell.rows)
    assert sum(count for _divisibility, count in key[4]) == len(cell.rows)


def test_configuration_invariant_cache_is_bounded_and_identity_keyed() -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    ambient = lattice(gram.rows())
    cells = (
        LorentzianPerfectCell(ambient, ((1, 0), (0, 1)), ZERO_SQUARES),
        LorentzianPerfectCell(ambient, ((0, 1), (1, 0)), ZERO_SQUARES),
        LorentzianPerfectCell(ambient, ((-1, 0), (0, -1)), ZERO_SQUARES),
    )
    backend = LorentzianPerfectLocalBackend()
    first = backend._integral_configuration_invariants(cells[0])

    assert backend._integral_configuration_invariants(cells[0]) is first
    backend._integral_configuration_invariants(cells[1])
    assert backend._integral_configuration_invariants(cells[0]) is first
    backend._integral_configuration_invariants(cells[2])

    cache = backend.__dict__["_configuration_invariant_cache"]
    assert len(cache) == 2
    assert id(cells[0]) in cache
    assert id(cells[1]) not in cache
    assert cache[id(cells[0])][0] is cells[0]
    assert cache[id(cells[2])][0] is cells[2]


def test_owned_large_cell_cone_releases_normaliz_engine_cache() -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])

    class LargeInput:
        pass

    cell = LargeInput()
    cell.lattice = lattice(gram.rows())
    cell.rows = ((1, 0), (0, 1)) * 5_000
    backend = LorentzianPerfectLocalBackend()
    cone = backend.cell_cone(cell)
    assert cone.primitive_rays()
    assert cell not in backend.__dict__.get("_normaliz_cell_cones", {})


def test_large_cell_owned_cone_bypasses_normaliz_polyhedron_wrapper(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])

    class LargeInput:
        pass

    cell = LargeInput()
    cell.lattice = lattice(gram.rows())
    cell.rows = ((1, 0), (0, 1)) * 5_000
    assert len(cell.rows) >= 10_000
    backend = LorentzianPerfectLocalBackend()
    monkeypatch.setattr(
        backend,
        "_normaliz_cell_polyhedron",
        lambda _cell: pytest.fail("large owned-cone conversion materialized a Sage polyhedron"),
    )
    cone = backend.cell_cone(cell)
    assert {tuple(ray) for ray in cone.primitive_rays()} == {(1, 0), (0, 1)}
    assert len(cone.halfspace_covectors()) == 2
    assert cell not in backend.__dict__.get("_normaliz_cell_cones", {})


def test_large_cell_owned_cone_preserves_normaliz_equations(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])

    class LargeInput:
        pass

    cell = LargeInput()
    cell.lattice = lattice(gram.rows())
    # A one-dimensional cone embedded in rank two has a nontrivial equation;
    # support hyperplanes alone describe a larger half-plane.
    cell.rows = ((1, 0),) * 10_000
    backend = LorentzianPerfectLocalBackend()
    original_h_presentation = backend._normaliz_cell_h_presentation
    h_presentation_calls = []

    def record_h_presentation(requested_cell, *, include_equations):
        h_presentation_calls.append((requested_cell, include_equations))
        return original_h_presentation(
            requested_cell, include_equations=include_equations
        )

    monkeypatch.setattr(backend, "_normaliz_cell_h_presentation", record_h_presentation)
    monkeypatch.setattr(
        backend,
        "_normaliz_cell_polyhedron",
        lambda _cell: pytest.fail("large owned-cone conversion materialized a Sage polyhedron"),
    )
    cone = backend.cell_cone(cell)
    rays = tuple(cone.primitive_rays())
    equations = cone.equation_covectors()
    assert {tuple(ray) for ray in rays} == {(1, 0)}
    assert equations
    labels = tuple(cell.lattice.module_generating_set())
    assert all(
        sum(equation(label) * ray(label) for label in labels) == 0
        for equation in equations
        for ray in rays
    )
    assert h_presentation_calls == [(cell, True)]


def test_large_cell_ray_extraction_bypasses_normaliz_polyhedron_wrapper(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])

    class LargeInput:
        pass

    cell = LargeInput()
    cell.lattice = lattice(gram.rows())
    cell.rows = ((1, 0), (0, 1)) * 5_000
    backend = LorentzianPerfectLocalBackend()
    monkeypatch.setattr(
        backend,
        "_normaliz_cell_polyhedron",
        lambda _cell: pytest.fail("large ray extraction materialized a Sage polyhedron"),
    )
    rays = backend.cell_primitive_rays(cell)
    assert {tuple(ray) for ray in rays} == {(1, 0), (0, 1)}
    assert cell not in backend.__dict__.get("_normaliz_cell_cones", {})


def test_facet_incidence_streaming_preserves_exact_orbits() -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    backend = LorentzianPerfectLocalBackend()
    labels = tuple(cell.lattice.module_generating_set())
    covectors = tuple(backend._iter_cell_halfspace_covectors(cell))
    pairing_rows = tuple(backend._iter_cell_halfspace_pairing_rows(cell, labels))
    assert pairing_rows == tuple(tuple(int(covector(label)) for label in labels) for covector in covectors)
    permutation_orbits = backend.facet_orbit_permutations(cell)
    permutation_orbits_again = backend.facet_orbit_permutations(cell)
    assert permutation_orbits_again == permutation_orbits
    assert permutation_orbits_again is not permutation_orbits
    assert all(facet != frozenset(range(len(cell.rows))) for orbit in permutation_orbits for facet, _ in orbit)
    isometry_orbits = backend.facet_orbits(cell)
    isometry_orbits_again = backend.facet_orbits(cell)
    assert isometry_orbits_again == isometry_orbits
    assert isometry_orbits_again is not isometry_orbits
    word_orbits, generator_permutations = backend._facet_orbit_transversal_words(cell)
    mask_orbit_iterator, mask_generators = backend._facet_orbit_transversal_masks(cell)
    mask_orbits = tuple(mask_orbit_iterator)
    all_ray_mask = (1 << len(cell.rows)) - 1
    assert all(isinstance(orbit, list) for orbit in mask_orbits)
    assert mask_generators == generator_permutations
    assert all(isinstance(permutation, array) for permutation in mask_generators)
    assert tuple(
        tuple((_facet_incidence(mask, len(cell.rows)), word) for mask, word in orbit)
        for orbit in mask_orbits
    ) == word_orbits

    assert tuple(tuple(facet for facet, _ in orbit) for orbit in permutation_orbits) == tuple(
        tuple(facet for facet, _ in orbit) for orbit in isometry_orbits
    )
    assert tuple(tuple(facet for facet, _ in orbit) for orbit in word_orbits) == tuple(
        tuple(facet for facet, _ in orbit) for orbit in permutation_orbits
    )
    for word_orbit, permutation_orbit in zip(word_orbits, permutation_orbits, strict=True):
        for (facet, word), (permutation_facet, permutation) in zip(word_orbit, permutation_orbit, strict=True):
            assert facet == permutation_facet
            composed = tuple(range(len(cell.vector_configuration)))
            generator_indices = []
            node = word
            while node is not None:
                node, generator_index = node
                generator_indices.append(generator_index)
            for generator_index in reversed(generator_indices):
                generator = generator_permutations[generator_index]
                composed = tuple(generator[position] for position in composed)
            assert composed == permutation
    for permutation_orbit, isometry_orbit in zip(permutation_orbits, isometry_orbits, strict=True):
        for (facet, permutation), (isometry_facet, isometry) in zip(
            permutation_orbit, isometry_orbit, strict=True
        ):
            assert facet == isometry_facet
            for position, image_position in enumerate(permutation):
                assert isometry(cell.vector_configuration[position]) == cell.vector_configuration[image_position]


def test_cell_faces_reuse_owned_halfspace_covectors_exactly() -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    lattice_ = lattice(gram.rows())
    cell = LorentzianPerfectCell(lattice_, ((1, 0), (0, 1)), ZERO_SQUARES)
    backend = LorentzianPerfectLocalBackend()
    cone = backend.cell_cone(cell)
    structure = _CellFaces(backend, cell)
    assert structure.gram == gram
    labels = tuple(lattice_.module_generating_set())
    assert backend.point_permutations(cell, cell.vector_configuration) == backend.point_permutations_from_rows(cell, cell.rows)

    def coordinates(covector):
        return tuple(int(covector.to_vector()(label)) for label in lattice_.module_generating_set())

    assert {tuple(int(entry) for entry in row) for row in structure.rows} == {
        tuple(int(entry) for entry in row)
        for row in (tuple(int(ray.to_vector()(label)) for label in labels) for ray in cone.primitive_rays())
    }
    assert set(backend.cell_primitive_rays(cell)) == set(cone.primitive_rays())
    assert set(backend.cell_primitive_rays(cell, square=0)) == {
        ray for ray in cone.primitive_rays() if ray.q() == 0
    }
    assert {tuple(int(entry) for entry in row) for row in structure.covectors} == {
        coordinates(covector) for covector in cone.facet_covectors()
    }
    assert {
        tuple(int(covector.to_vector()(label)) for label in labels)
        for covector in backend.cell_halfspace_covectors(cell)
    } == {coordinates(covector) for covector in cone.halfspace_covectors()}
    assert structure.incidences == frozenset(
        frozenset(index for index, row in enumerate(structure.rows) if covector.dot_product(row) == 0)
        for covector in structure.covectors
    )
    interior_mark = structure.mark(vector(ZZ, (1, 1)))
    assert interior_mark.chain == (structure.full_face,)
    assert interior_mark.coefficients == (1,)
    boundary_mark = structure.mark(vector(ZZ, (1, 0)))
    assert boundary_mark.chain == (frozenset({0}),)
    assert boundary_mark.coefficients == (1,)
    mixed_mark = structure.mark(vector(ZZ, (2, 1)))
    assert mixed_mark.chain == (structure.full_face, frozenset({0}))
    assert mixed_mark.coefficients == (1, 1)
