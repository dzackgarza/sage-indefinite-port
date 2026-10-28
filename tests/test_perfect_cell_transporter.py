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
    _class_sum_pairing_signatures,
    _configuration_signature_multisets_match,
    _facet_complement_key,
    _facet_incidence,
    _gram_component_sign_flips,
    _integral_basis_reflection,
    _iter_owned_row_blocks,
    _permuted_bitset,
    _permuted_facet_key,
    _packed_pairing_fibres,
    _perfect_form_hash_key,
    _schreier_generator_indices,
    LorentzianPerfectComplex,
    LorentzianPerfectCell,
    LorentzianPerfectLocalBackend,
    NONNEGATIVE_SQUARES,
    ZERO_SQUARES,
    lorentzian_orthogonal_group_generators,
    perfect_complex,
)
import sage_indefinite_port.indefinite.lorentzian_cells as lorentzian_cells
from tests.acceptance.consumer import lattice, rows_of


def test_group_generator_complex_does_not_alias_full_domain_cache() -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    ambient = lattice(gram.rows())

    full_domain = perfect_complex(ambient, ZERO_SQUARES)
    group_only = perfect_complex(ambient, ZERO_SQUARES, retain_facet_orbits=False)

    assert full_domain is perfect_complex(ambient, ZERO_SQUARES)
    assert full_domain is perfect_complex(ambient, ZERO_SQUARES, retain_facet_orbits=True)
    assert group_only is perfect_complex(ambient, ZERO_SQUARES, retain_facet_orbits=False)
    assert full_domain is not group_only
    assert full_domain._retain_facet_orbits
    assert not group_only._retain_facet_orbits


def test_group_only_expansion_discards_facet_orbit_records(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    for retain in (True, False):
        ambient = lattice(gram.rows())
        cell = LorentzianPerfectCell(ambient, ((1, 0), (0, 1)), ZERO_SQUARES)
        complex_ = LorentzianPerfectComplex(ambient, ZERO_SQUARES, retain_facet_orbits=retain)
        backend = complex_.local_backend()
        key = _perfect_form_hash_key(cell)
        complex_._traversal_representatives = [cell]
        complex_._traversal_representative_buckets = {key: [cell]}
        monkeypatch.setattr(
            backend,
            "_facet_mask_data",
            lambda _cell, *, use_owned_cone=True: (iter(()), (), None),
        )
        monkeypatch.setattr(
            backend,
            "_iter_facet_mask_orbits",
            lambda _facets, _permutations, *, retain_words=True, canonical_keys=None: iter(
                [[((True, (0,)), None)]]
            ),
        )
        monkeypatch.setattr(backend, "flip_across", lambda _cell, _facet: cell)
        monkeypatch.setattr(backend, "cell_transporter", lambda _source, _target: ambient.Aut().identity())

        complex_._expand_traversal_once()

        assert len(complex_._traversal_adjacencies) == 1
        if retain:
            assert complex_._traversal_facet_orbits[id(cell)]
            assert id(cell) in complex_._traversal_facet_generators
        else:
            assert complex_._traversal_facet_orbits == {}
            assert complex_._traversal_facet_generators == {}


def test_orthogonal_group_generators_use_raw_h_rows(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1, 0], [1, 0, 0], [0, 0, -2]])
    ambient = lattice(gram.rows())
    monkeypatch.setattr(
        LorentzianPerfectLocalBackend,
        "cell_cone",
        lambda _self, _cell: pytest.fail("group-only orthogonal traversal must not construct an owned cone"),
    )

    generators = lorentzian_orthogonal_group_generators(ambient)

    assert generators
    for generator in generators:
        action = matrix(ZZ, rows_of(generator))
        assert abs(action.det()) == 1
        assert action * gram * action.transpose() == gram
    complex_ = perfect_complex(ambient, NONNEGATIVE_SQUARES, retain_facet_orbits=False)
    assert complex_._traversal_facet_orbits == {}
    assert complex_._traversal_facet_generators == {}


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


def test_facet_orbit_enumerator_streams_representatives_without_words() -> None:
    swap_first_two = array("H", [1, 0, 2])
    facets = ((False, (0,)), (False, (1,)), (False, (2,)))
    backend = LorentzianPerfectLocalBackend()

    retained = tuple(backend._iter_facet_mask_orbits(iter(facets), (swap_first_two,)))
    canonical = {key: key for key in facets}
    canonical_retained = tuple(
        backend._iter_facet_mask_orbits(
            iter(facets), (swap_first_two,), canonical_keys=canonical
        )
    )
    streamed = tuple(
        backend._iter_facet_mask_orbits(iter(facets), (swap_first_two,), retain_words=False)
    )

    assert tuple(len(orbit) for orbit in retained) == (2, 1)
    assert tuple(tuple(entry[0] for entry in orbit) for orbit in retained) == (
        ((False, (0,)), (False, (1,))),
        ((False, (2,)),),
    )
    assert canonical_retained[0][1][0] is facets[1]
    assert all(value is None for value in canonical.values())
    with pytest.raises(ArithmeticError, match="outside the exact H-presentation"):
        tuple(
            backend._iter_facet_mask_orbits(
                iter((facets[0],)),
                (swap_first_two,),
                canonical_keys={facets[0]: facets[0]},
            )
        )
    assert tuple(orbit[0][0] for orbit in streamed) == ((False, (0,)), (False, (2,)))
    assert all(len(orbit) == 1 for orbit in streamed)
    assert _permuted_facet_key(facets[0], swap_first_two) == facets[1]


def test_facet_orbit_materializers_release_consumed_bfs_entries(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    backend = LorentzianPerfectLocalBackend()
    original = backend._iter_facet_mask_orbits
    emitted = []

    def recording_orbits(facet_keys, permutations, *, retain_words=True, canonical_keys=None):
        for orbit in original(
            facet_keys,
            permutations,
            retain_words=retain_words,
            canonical_keys=canonical_keys,
        ):
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


def test_cell_stabilizer_seeds_certified_gram_automorphisms() -> None:
    gram = matrix(ZZ, [[2, -1, 0], [-1, 2, 0], [0, 0, -2]])
    roots = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1))
    rows = tuple((x, y, 1) for x, y in roots)
    cell = LorentzianPerfectCell(lattice(gram.rows()), rows, ZERO_SQUARES)

    generators = LorentzianPerfectLocalBackend().cell_stabilizer(cell)
    actions = {tuple(tuple(int(entry) for entry in row) for row in rows_of(g)) for g in generators}
    expected = {
        tuple(tuple(int(entry) for entry in row) for row in _integral_basis_reflection(gram, index).rows())
        for index in range(2)
    }
    expected.add(tuple(tuple(int(entry) for entry in row) for row in matrix(ZZ, [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]).rows()))

    assert expected <= actions
    for action_rows in actions:
        action = matrix(ZZ, action_rows)
        assert abs(action.det()) == 1
        assert action * gram * action.transpose() == gram
        assert set(map(tuple, matrix(ZZ, cell.rows) * action)) == set(cell.rows)


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
    original_counter = lorentzian_cells.Counter

    def reject_row_multiset_counter(iterable=(), *args, **kwargs):
        code = getattr(iterable, "gi_code", None)
        if iterable is cell.rows or (
            code is not None
            and "_configuration_signature_multisets_match" not in code.co_qualname
            and "backtrack" not in code.co_qualname
        ):
            pytest.fail("large stabilizer certification must use its exact row-position map")
        return original_counter(iterable, *args, **kwargs)

    monkeypatch.setattr(lorentzian_cells, "Counter", reject_row_multiset_counter)

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


def test_class_sum_pairing_colours_are_isometry_invariants() -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    rows = matrix(ZZ, [[1, 0], [0, 1], [1, 1], [-1, -1]])
    signatures = ((0,), (0,), (1,), (1,))
    action = matrix(ZZ, [[0, 1], [1, 0]])

    source = _class_sum_pairing_signatures(rows, gram, signatures)
    target = _class_sum_pairing_signatures(rows * action, gram, signatures)

    assert len(source[0]) == len(signatures[0]) + 2
    assert lorentzian_cells.Counter(source) == lorentzian_cells.Counter(target)
    assert _class_sum_pairing_signatures(rows, gram, signatures, max_classes=1) == signatures


def test_cell_stabilizer_seeds_basis_transpositions() -> None:
    gram = matrix(ZZ, [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
    rows = ((2, 0, 0, 2), (-2, 0, 0, 2), (0, 2, 0, 2), (0, -2, 0, 2), (0, 0, 2, 2), (0, 0, -2, 2))
    cell = LorentzianPerfectCell(lattice(gram.rows()), rows, ZERO_SQUARES)

    generators = LorentzianPerfectLocalBackend().cell_stabilizer(cell)
    expected = matrix(ZZ, [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    actions = {tuple(tuple(int(entry) for entry in row) for row in rows_of(generator)) for generator in generators}

    assert tuple(tuple(int(entry) for entry in row) for row in expected.rows()) in actions
    assert expected * gram * expected.transpose() == gram
    assert set(map(tuple, matrix(ZZ, rows) * expected)) == set(rows)


def test_gram_component_sign_flips_are_exact_isometries() -> None:
    gram = matrix(ZZ, [[2, -1, 0, 0], [-1, 2, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])
    actions = _gram_component_sign_flips(gram)

    assert len(actions) == 3
    assert all(abs(action.det()) == 1 and action * gram * action.transpose() == gram for action in actions)
    assert {tuple(action.diagonal()) for action in actions} == {
        (-1, -1, 1, 1),
        (1, 1, -1, 1),
        (1, 1, 1, -1),
    }


def test_integral_basis_reflection_is_exact_or_rejected() -> None:
    root_gram = matrix(ZZ, [[2, -1], [-1, 2]])
    reflection = _integral_basis_reflection(root_gram, 0)
    assert reflection == matrix(ZZ, [[-1, 0], [1, 1]])
    assert abs(reflection.det()) == 1
    assert reflection * root_gram * reflection.transpose() == root_gram

    nonintegral_gram = matrix(ZZ, [[4, 1], [1, 4]])
    assert _integral_basis_reflection(nonintegral_gram, 0) is None
    null_gram = matrix(ZZ, [[0, 1], [1, 0]])
    assert _integral_basis_reflection(null_gram, 0) is None


def test_configuration_pairing_diagonal_matches_exact_row_products() -> None:
    rows = matrix(ZZ, [[1, 2], [3, 4], [5, 6], [-2, 7]])
    gram = matrix(ZZ, [[2, 1], [1, -3]])
    pairings = lorentzian_cells._ConfigurationPairings(rows, gram)

    expected = tuple((rows.row(i) * gram).dot_product(rows.row(i)) for i in range(rows.nrows()))

    assert pairings.diagonal() == expected


def test_configuration_pairing_diagonal_crosses_temporary_block_boundary() -> None:
    rows = matrix(ZZ, [[index, index * index - 3] for index in range(1025)])
    gram = matrix(ZZ, [[2, 1], [1, -3]])
    pairings = lorentzian_cells._ConfigurationPairings(rows, gram)

    expected = tuple((rows.row(i) * gram).dot_product(rows.row(i)) for i in range(rows.nrows()))

    assert len(pairings.diagonal()) == rows.nrows()
    assert pairings.diagonal() == expected


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


def test_configuration_signature_multisets_match_exact_columns() -> None:
    colours = (0, 1, 0)
    columns = ((10, 11, 12), (30, 31, 32), (20, 21, 22))
    source_counts = lorentzian_cells.Counter(
        (colour, tuple(values[index] for values in columns))
        for index, colour in enumerate(colours)
    )
    original_counts = source_counts.copy()
    assert _configuration_signature_multisets_match(
        colours, columns, colours, columns, source_signature_counts=source_counts
    )
    assert source_counts == original_counts
    changed_columns = (*columns[:-1], (20, 21, 99))
    assert not _configuration_signature_multisets_match(
        colours, columns, colours, changed_columns, source_signature_counts=source_counts
    )
    assert source_counts == original_counts


def test_packed_pairing_fibres_preserve_exact_position_partition() -> None:
    fibres = _packed_pairing_fibres((0, 1, 0, 0, 1), (7, -2, 7, 9, -2), "I")

    assert {key: tuple(indices) for key, indices in fibres.items()} == {
        (0, 7): (0, 2),
        (1, -2): (1, 4),
        (0, 9): (3,),
    }
    assert all(isinstance(indices, array) and indices.typecode == "I" for indices in fibres.values())
    with pytest.raises(ValueError, match="one colour per target row"):
        _packed_pairing_fibres((0,), (7, 8), "I")


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

    # When the incident side is the smaller dense side, retain its bitset so
    # Schreier traversal permutes 2,000 bits rather than 47,861.
    incident_mask = sum(1 << position for position in range(2_000))
    incident_key = _facet_complement_key(((1 << ray_count) - 1) ^ incident_mask, ray_count)
    assert incident_key == (False, incident_mask)
    assert _facet_incidence(incident_key, ray_count) == frozenset(range(2_000))
    incident_image = _permuted_facet_key(incident_key, reverse)
    assert incident_image == (False, sum(1 << (ray_count - 1 - position) for position in range(2_000)))


def test_normaliz_pairing_rows_match_owned_covectors(monkeypatch) -> None:
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
    original_h_presentation = backend._normaliz_cell_h_presentation
    original_normaliz_cone = backend._normaliz_cell_cone
    h_presentation_calls = 0
    cached_normaliz_cone_calls = 0

    def count_h_presentation(*args, **kwargs):
        nonlocal h_presentation_calls
        h_presentation_calls += 1
        return original_h_presentation(*args, **kwargs)

    def count_normaliz_cone(*args, **kwargs):
        nonlocal cached_normaliz_cone_calls
        cached_normaliz_cone_calls += 1
        return original_normaliz_cone(*args, **kwargs)

    monkeypatch.setattr(backend, "_normaliz_cell_h_presentation", count_h_presentation)
    monkeypatch.setattr(backend, "_normaliz_cell_cone", count_normaliz_cone)
    labels = tuple(cell.lattice.module_generating_set())
    pairing_rows = tuple(backend._iter_cell_halfspace_pairing_rows(cell, labels))
    raw_pairing_rows = tuple(
        backend._iter_cell_halfspace_pairing_rows(cell, labels, use_owned_cone=False)
    )
    assert cell not in backend.__dict__.get("_normaliz_cell_cones", {})
    owned_rows = tuple(
        tuple(int(covector(label)) for label in labels)
        for covector in backend._iter_cell_halfspace_covectors(cell)
    )
    # The primitive-covector view now independently consumes the raw H-data
    # rather than forcing the owned cone; the owned, raw-pairing and primitive
    # views therefore make three exact H requests in this explicit comparison.
    assert h_presentation_calls == 3
    assert cached_normaliz_cone_calls == 1
    assert cell not in backend.__dict__.get("_normaliz_cell_cones", {})
    def zero_pattern(rows):
        return tuple(tuple(value == 0 for value in row) for row in rows)

    assert zero_pattern(pairing_rows) == zero_pattern(owned_rows)
    assert zero_pattern(raw_pairing_rows) == zero_pattern(owned_rows)


def test_group_only_facet_keys_do_not_construct_owned_cone(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    backend = LorentzianPerfectLocalBackend()
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)
    monkeypatch.setattr(
        backend,
        "cell_cone",
        lambda _cell: pytest.fail("group-only facet enumeration must consume raw Normaliz H-rows"),
    )
    original_h_presentation = backend._normaliz_cell_h_presentation
    raw_h_rows = []

    def record_h_rows(_cell, *, include_equations):
        assert id(_cell) not in backend.__dict__.get("_configuration_invariant_cache", {})
        hyperplanes, equations = original_h_presentation(_cell, include_equations=include_equations)
        raw_h_rows.append(hyperplanes)
        return hyperplanes, equations

    monkeypatch.setattr(backend, "_normaliz_cell_h_presentation", record_h_rows)

    facet_keys, _permutations, canonical_keys = backend._facet_mask_data(cell, use_owned_cone=False)

    assert raw_h_rows and all(row is None for row in raw_h_rows[0])
    assert id(cell) not in backend.__dict__.get("_configuration_invariant_cache", {})
    assert canonical_keys is not None
    assert all(canonical_keys[key] is key for key in canonical_keys)
    assert set(facet_keys) == {1, 2}


def test_public_facet_orbit_permutations_keep_owned_cone_lazy(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    backend = LorentzianPerfectLocalBackend()
    monkeypatch.setattr(
        backend,
        "cell_cone",
        lambda _cell: pytest.fail("facet orbit permutations need only raw H incidence keys"),
    )

    orbits = backend.facet_orbit_permutations(cell)

    assert len(orbits) == 1
    assert {facet for facet, _permutation in orbits[0]} == {frozenset({0}), frozenset({1})}
    assert {permutation for _facet, permutation in orbits[0]} == {(0, 1), (1, 0)}


def test_perfect_traversal_keeps_owned_cone_lazy(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    complex_ = perfect_complex(lattice(gram.rows()), ZERO_SQUARES)
    backend = complex_.local_backend()
    monkeypatch.setattr(
        backend,
        "cell_cone",
        lambda _cell: pytest.fail("facet traversal needs incidence rows, not the owned cone"),
    )

    complex_._ensure_traversal_started()
    assert complex_._expand_traversal_once() == ()
    assert backend.__dict__.get("_normaliz_cell_cones", {}) == {}


def test_raw_normaliz_facet_route_matches_owned_cone_orbits() -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    backend = LorentzianPerfectLocalBackend()

    owned_keys, owned_generators, _ = backend._facet_mask_data(cell, use_owned_cone=True)
    raw_keys, raw_generators, raw_canonical = backend._facet_mask_data(cell, use_owned_cone=False)
    owned_keys = tuple(owned_keys)
    raw_keys = tuple(raw_keys)
    assert raw_canonical is not None
    assert tuple(tuple(generator) for generator in owned_generators) == tuple(tuple(generator) for generator in raw_generators)
    owned_facets = {_facet_incidence(key, len(cell.rows)) for key in owned_keys}
    raw_facets = {_facet_incidence(key, len(cell.rows)) for key in raw_keys}
    assert raw_facets == owned_facets

    owned_orbits = tuple(
        tuple(sorted(_facet_incidence(key, len(cell.rows)) for key, _word in orbit))
        for orbit in backend._iter_facet_mask_orbits(iter(owned_keys), owned_generators)
    )
    raw_orbits = tuple(
        tuple(sorted(_facet_incidence(key, len(cell.rows)) for key, _word in orbit))
        for orbit in backend._iter_facet_mask_orbits(
            iter(raw_keys), raw_generators, canonical_keys=raw_canonical
        )
    )
    assert {frozenset(orbit) for orbit in raw_orbits} == {frozenset(orbit) for orbit in owned_orbits}


def test_large_facet_keys_use_the_smaller_dense_incidence_side(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    rows = tuple((index + 1, 0) for index in range(1_000)) + tuple(
        (0, index + 1) for index in range(9_000)
    )

    class LargeInput:
        pass

    cell = LargeInput()
    cell.lattice = lattice(gram.rows())
    cell.rows = rows
    backend = LorentzianPerfectLocalBackend()
    monkeypatch.setattr(backend, "cell_stabilizer", lambda _cell: ())

    facet_keys, _permutations, canonical_keys = backend._facet_mask_data(
        cell, use_owned_cone=False
    )
    keys = tuple(facet_keys)

    assert canonical_keys is not None
    assert len(keys) == 2
    dense_incident_keys = [
        key for key in keys if isinstance(key, tuple) and isinstance(key[1], int)
    ]
    assert len(dense_incident_keys) == 1
    assert len(_facet_incidence(dense_incident_keys[0], len(cell.rows))) == 1_000


def test_large_perfect_form_hash_skips_quadratic_pairing_graph(monkeypatch) -> None:
    """Large-cell bucketing uses exact linear-cost integral invariants."""
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
    monkeypatch.setattr(
        lorentzian_cells._ConfigurationPairings,
        "diagonal",
        lambda _self: pytest.fail("large perfect-form hashing must stream diagonal blocks"),
    )
    key = _perfect_form_hash_key(cell)
    backend = LorentzianPerfectLocalBackend()
    refined_key = _perfect_form_hash_key(cell, backend)
    # The backend has the cached linear-cost invariant data and may safely
    # append stronger exact colours; the uncached route remains a coarser key.
    assert refined_key[:5] == key
    assert len(refined_key) == 6
    assert id(cell) in backend.__dict__["_configuration_invariant_cache"]

    assert key[1] == ()
    assert sum(count for _value, count in key[0]) == len(cell.rows)
    assert len(key[2]) == 2
    assert sum(count for _value, count in key[3]) == len(cell.rows)
    assert sum(count for _value, count in key[4]) == len(cell.rows)
    assert sum(count for _value, count in refined_key[5]) == len(cell.rows)


def test_complex_isometry_uses_cached_coarse_key_for_large_route(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    ambient = lattice(gram.rows())
    source = LorentzianPerfectComplex(ambient, ZERO_SQUARES)
    target = LorentzianPerfectComplex(ambient, ZERO_SQUARES)
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)
    monkeypatch.setattr(
        lorentzian_cells._ConfigurationPairings,
        "blocks",
        lambda _self: pytest.fail("large-route complex isometry must not build the quadratic pairing graph"),
    )
    monkeypatch.setattr(
        lorentzian_cells._ConfigurationPairings,
        "diagonal",
        lambda _self: pytest.fail("large-route complex isometry must consume diagonal blocks"),
    )

    witness = source.isometry_to(target)

    assert witness is not None
    assert abs(matrix(ZZ, rows_of(witness)).det()) == 1


def test_large_distinct_row_transporter_uses_exact_target_membership(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    source = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    target_gram = -gram
    target = LorentzianPerfectCell(lattice(target_gram.rows()), ((1, 0), (0, -1)), ZERO_SQUARES)
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)
    original_counter = lorentzian_cells.Counter

    def reject_configuration_counter(iterable=(), *args, **kwargs):
        if iterable is source.rows or iterable is target.rows:
            pytest.fail("large distinct-row transporters must not build row-multiset Counters")
        return original_counter(iterable, *args, **kwargs)

    monkeypatch.setattr(lorentzian_cells, "Counter", reject_configuration_counter)
    witness = LorentzianPerfectLocalBackend().cell_transporter(source, target)

    assert witness is not None
    assert abs(matrix(ZZ, rows_of(witness)).det()) == 1


def test_large_transporter_matches_joint_gcd_colours_after_basis_change(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    shear = matrix(ZZ, [[1, 2**80 + 1], [0, 1]])
    target_gram = shear.inverse() * gram * shear.inverse().transpose()
    source_rows = matrix(ZZ, [[2, 0], [0, 3]])
    source = LorentzianPerfectCell(
        lattice(gram.rows()), tuple(tuple(row) for row in source_rows), ZERO_SQUARES
    )
    target = LorentzianPerfectCell(
        lattice(target_gram.rows()), tuple(tuple(row) for row in source_rows * shear), ZERO_SQUARES
    )
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)

    backend = LorentzianPerfectLocalBackend()
    source_spectra = backend._integral_configuration_invariants(source)[2]
    target_spectra = backend._integral_configuration_invariants(target)[2]
    assert lorentzian_cells.Counter(source_spectra) == lorentzian_cells.Counter(target_spectra)
    assert _perfect_form_hash_key(source, backend) == _perfect_form_hash_key(target, backend)
    witness = backend.cell_transporter(source, target)

    assert witness is not None
    action = matrix(ZZ, rows_of(witness))
    assert action == shear
    assert abs(action.det()) == 1
    assert action * target_gram * action.transpose() == gram


def test_large_same_lattice_transporter_skips_counter_identity_probe(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    ambient = lattice(gram.rows())
    source = LorentzianPerfectCell(ambient, ((1, 0), (0, 1)), ZERO_SQUARES)
    target = LorentzianPerfectCell(ambient, ((0, 1), (1, 0)), ZERO_SQUARES)
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)
    original_counter = lorentzian_cells.Counter

    def reject_row_counter(iterable=(), *args, **kwargs):
        if iterable is source.rows or iterable is target.rows:
            pytest.fail("large same-lattice transporters must not build row-multiset Counters")
        return original_counter(iterable, *args, **kwargs)

    monkeypatch.setattr(lorentzian_cells, "Counter", reject_row_counter)
    witness = LorentzianPerfectLocalBackend().cell_transporter(source, target)

    assert witness is not None
    action = matrix(ZZ, rows_of(witness))
    assert abs(action.det()) == 1
    assert action * gram * action.transpose() == gram


def test_repeated_row_transporter_retains_multiset_certificate(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    ambient = lattice(gram.rows())
    source = LorentzianPerfectCell(ambient, ((1, 0), (0, 1), (1, 0)), ZERO_SQUARES)
    target = LorentzianPerfectCell(ambient, ((0, 1), (1, 0), (0, 1)), ZERO_SQUARES)
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)
    backend = LorentzianPerfectLocalBackend()
    assert not backend._integral_configuration_invariants(source)[6]

    witness = backend.cell_transporter(source, target)

    assert witness is not None
    action = matrix(ZZ, rows_of(witness))
    assert abs(action.det()) == 1
    assert action * gram * action.transpose() == gram


def test_large_configuration_invariants_include_exact_gcd_colours(monkeypatch) -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    monkeypatch.setattr(lorentzian_cells, "_PAIRING_GRAPH_REFINEMENT_LIMIT", 0)
    backend = LorentzianPerfectLocalBackend()

    _rows, _pairings, spectra, _smith, row_sum, row_moment, distinct_rows, row_positions = backend._integral_configuration_invariants(cell)

    assert spectra == ((0, 1, 1, 1, 1, 0), (0, 1, 1, 1, 1, 0))
    assert tuple(row_sum) == (1, 1)
    assert row_moment == matrix(ZZ, [[1, 0], [0, 1]])
    assert distinct_rows
    assert row_positions == {(1, 0): 0, (0, 1): 1}
    assert len(backend._integral_configuration_invariants(cell)) == 8


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
    assert backend.cell_halfspace_covectors(cell) is cone.halfspace_covectors()
    structure = _CellFaces(backend, cell)
    assert structure.gram == gram
    labels = tuple(lattice_.module_generating_set())
    assert backend.point_permutations(cell, cell.vector_configuration) == backend.point_permutations_from_rows(cell, cell.rows)
    assert not hasattr(structure, "covectors")

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
    assert {tuple(int(entry) for entry in row) for row in structure._covector_matrix.rows()} == {
        coordinates(covector) for covector in cone.facet_covectors()
    }
    assert {
        tuple(int(covector.to_vector()(label)) for label in labels)
        for covector in backend.cell_halfspace_covectors(cell)
    } == {coordinates(covector) for covector in cone.halfspace_covectors()}
    assert structure.incidences == frozenset(
        frozenset(
            index
            for index, row in enumerate(structure.rows)
            if sum(a * b for a, b in zip(covector, row, strict=True)) == 0
        )
        for covector in structure._covector_matrix.rows()
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
