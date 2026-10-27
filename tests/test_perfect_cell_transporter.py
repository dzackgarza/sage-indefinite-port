"""Integral transport of perfect cells in U + E8(-2).

The explicit unimodular matrix certifies that the two complete configurations
are isometric, independently of the transporter search.
"""

from sage.all import ZZ, matrix, vector

from sage_indefinite_port.indefinite.lorentzian_cells import (
    _CellFaces,
    LorentzianPerfectCell,
    LorentzianPerfectLocalBackend,
    NONNEGATIVE_SQUARES,
    ZERO_SQUARES,
    _extend_configuration_signatures,
)
from tests.acceptance.consumer import lattice, rows_of

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


def test_anchor_signature_extension_appends_in_common_assignment_order() -> None:
    signatures = ((0, (10, 30)), (1, (11, 31)), (0, (12, 32)))
    column = (20, 21, 22)
    assert _extend_configuration_signatures(signatures, column) == (
        (0, (10, 30, 20)),
        (1, (11, 31, 21)),
        (0, (12, 32, 22)),
    )


def test_facet_incidence_streaming_preserves_exact_orbits() -> None:
    gram = matrix(ZZ, [[0, 1], [1, 0]])
    cell = LorentzianPerfectCell(lattice(gram.rows()), ((1, 0), (0, 1)), ZERO_SQUARES)
    backend = LorentzianPerfectLocalBackend()
    permutation_orbits = backend.facet_orbit_permutations(cell)
    isometry_orbits = backend.facet_orbits(cell)
    word_orbits, generator_permutations = backend._facet_orbit_transversal_words(cell)

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
    labels = tuple(lattice_.module_generating_set())

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
