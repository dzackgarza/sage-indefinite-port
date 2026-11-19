"""Complete, certified perfect-cell stabilizers from definite second moments."""
from math import factorial
import json
from pathlib import Path

from sage.matrix.constructor import matrix
from sage.rings.integer_ring import ZZ as SageZZ
from sage.groups.perm_gps.permgroup_named import SymmetricGroup
from sage_extensions.all import Lattices, ZZ
from sage_indefinite_port.indefinite.lorentzian_cells import (
    LorentzianPerfectCell,
    LorentzianPerfectLocalBackend,
    NONNEGATIVE_SQUARES,
    _gram_matrix,
)


def test_full_diagonal_moment_stabilizer_has_exact_hyperoctahedral_order():
    # The vectors (2, +/-e1), (2, +/-e2) lie in one affine plane.
    # Their moment is diag(16,2,2); their full stabilizer has order 8.
    ambient = Lattices(ZZ)([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])
    rays = ((2, 1, 0), (2, -1, 0), (2, 0, 1), (2, 0, -1))
    cell = LorentzianPerfectCell(ambient, rays, NONNEGATIVE_SQUARES)
    backend = LorentzianPerfectLocalBackend()
    generators = backend._diagonal_moment_stabilizer(cell)
    assert generators is not None and len(generators) == 3
    assert tuple(backend.cell_stabilizer(cell)) == generators
    permutations = backend.cell_stabilizer_permutations(cell)
    assert len(permutations) == 3
    finite = SymmetricGroup(len(rays)).subgroup(
        tuple(SymmetricGroup(len(rays))([i+1 for i in images]) for images in permutations)
    )
    assert int(finite.order()) == 8
    assert all(all(g(ambient.module_generator(label)).q() == ambient.module_generator(label).q()
                   for label in ambient.module_generating_set())
               for g in generators)
    # Reversing the time coordinate is in O(M), but does not
    # preserve the rays. The tested signed spatial group is exhaustive.
    assert all(row[0] > 0 for row in rays)


def test_frozen_rank_seventeen_initial_cell_stabilizer_is_certified_exactly():
    cases = json.loads(
        (Path(__file__).parent / "fixtures" / "isotropic_cases_8821.json").read_text()
    )
    data = next(item for item in cases if item["id"] == "IsotropicCases_8819")
    ambient = Lattices(ZZ)(data["gram"])
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(ambient, NONNEGATIVE_SQUARES)
    assert len(cell.rows) == 33

    rows = matrix(SageZZ, cell.rows)
    assert rows.transpose() * rows == matrix(
        SageZZ, 17, 17,
        lambda i, j: (33 if i == j == 0 else 2 if i == j else 0),
    )
    generators = backend.cell_stabilizer(cell)
    assert len(generators) == 31  # 16 independent signs and 15 adjacent swaps
    assert len(backend.cell_stabilizer_permutations(cell)) == 31
    positions = set(cell.rows)
    form = _gram_matrix(ambient)
    for generator, action in zip(
        generators, backend._stabilizer_actions[cell], strict=True
    ):
        assert abs(action.det()) == 1
        assert action * form * action.transpose() == form
        assert all(tuple(row * action) in positions for row in rows.rows())
        assert generator.domain() is ambient and generator.codomain() is ambient


def test_frozen_rank_thirteen_initial_cell_has_complete_integral_stabilizer():
    """The non-diagonal second moment is integrally congruent to diag(2049,486^12).

    Its 24 minimum vectors give twelve orthogonal primitive directions.
    Their saturated complement supplies one timelike direction, proving
    that the full moment group is C2 x (C2 wr S12). The ray-membership
    certificates exclude time reversal and retain the complete spatial
    subgroup, with no enumeration of the finite group's enormous order.
    """
    cases = json.loads(
        (Path(__file__).parent / "fixtures" / "isotropic_cases_8821.json").read_text()
    )
    data = next(item for item in cases if item["id"] == "IsotropicCases_8812")
    ambient = Lattices(ZZ)(data["gram"])
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(ambient, NONNEGATIVE_SQUARES)
    assert len(cell.rows) == 2049
    moment = matrix(SageZZ, cell.rows).transpose() * matrix(SageZZ, cell.rows)
    assert moment.det() == 2049 * 486 ** 12
    assert sum(int(moment[i,j] != 0)
               for i in range(moment.nrows()) for j in range(i+1,moment.ncols())) > 0
    generators = backend.cell_stabilizer(cell)
    assert len(generators) == 23  # 12 independent signs, 11 swaps
    actions = backend._stabilizer_actions[cell]
    assert len(actions) == len(generators)
    rows = matrix(SageZZ, cell.rows)
    positions = set(cell.rows)
    form = _gram_matrix(ambient)
    for generator, action in zip(generators, actions, strict=True):
        assert abs(action.det()) == 1
        assert action * form * action.transpose() == form
        assert action.transpose() * moment * action == moment
        assert all(tuple(int(value) for value in row) in positions
                   for row in (rows * action).rows())
        assert generator.domain() is ambient and generator.codomain() is ambient


def test_frozen_public_signature_twist_preserves_complete_moment_certificate():
    """The public dispatcher uses AttackProfile's opposite-signature lattice.

    The exact cell certificate must work there, not merely on the
    originally displayed indefinite Gram matrices of the frozen fixtures.
    """
    from sage_indefinite_port.invariants import AttackProfile
    from sage_indefinite_port.indefinite.indefinite_reduction import indefinite_reduced_presentation

    cases = json.loads(
        (Path(__file__).parent / "fixtures" / "isotropic_cases_8821.json").read_text()
    )
    by_id = {case["id"]: case for case in cases if case["id"] in (
        "IsotropicCases_8812", "IsotropicCases_8819",
    )}
    for case_id, expected in (("IsotropicCases_8812", 23), ("IsotropicCases_8819", 31)):
        original = Lattices(ZZ)(by_id[case_id]["gram"])
        profile = AttackProfile.from_lattice(original)
        reduced = indefinite_reduced_presentation(profile.signed_view).reduced
        assert int(tuple(reduced.signature_pair())[0]) == 1
        backend = LorentzianPerfectLocalBackend()
        cell = backend.initial_cell(reduced, NONNEGATIVE_SQUARES)
        generators = backend.cell_stabilizer(cell)
        assert len(generators) == expected
        assert len(backend.cell_stabilizer_permutations(cell)) == expected
        assert all(generator.domain() is reduced and generator.codomain() is reduced
                   for generator in generators)
        if case_id == "IsotropicCases_8812":
            # The actual signed-Weyl-chamber enumeration of this 2,049-ray
            # source cell has precisely two complete facet orbits.
            from sage_indefinite_port.indefinite.facet_orbits import _facet_orbit_records
            records = tuple(_facet_orbit_records(
                cell.rows, backend.cell_stabilizer_permutations(cell), retain_words=False
            ))
            assert len(records) == 2
            assert sorted((len(item.incidence), item.orbit_size) for item in records) == [
                (220, 4096), (243, 24),
            ]


def test_rank_seventeen_crosspolytope_facets_one_exact_stabilizer_orbit():
    """The 33-row cusp cell is a 16-dimensional crosspolytope.

    Its 32 extreme rays pair as (e0 +/- ei), i=1,...,16; the additional
    e0 ray lies in their positive interior. The facets are in bijection
    with choices of one sign for each coordinate, so there are 2^16.
    Independent coordinate signs act transitively. The recursive
    ridge walk must return the one orbit without first computing
    all 65,536 hyperplanes by direct dual description.
    """
    from sage_indefinite_port.indefinite.facet_orbits import _facet_orbit_records

    ambient = Lattices(ZZ)([[-1 if i == j == 0 else int(i == j)
                            for j in range(17)] for i in range(17)])
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(ambient, NONNEGATIVE_SQUARES)
    assert len(cell.rows) == 33
    generators = backend.cell_stabilizer(cell)
    assert len(generators) == 31
    permutations = backend.cell_stabilizer_permutations(cell)
    records = tuple(_facet_orbit_records(cell.rows, permutations, retain_words=True))
    assert len(records) == 1
    (facet,) = records
    assert facet.orbit_size == 2 ** 16
    assert len(facet.words) == 2 ** 16
    assert len(facet.incidence) == 16
    # Every stored Schreier word explicitly transports the seed facet
    # to the member it certifies.
    for member, word in (facet.words[0], facet.words[len(facet.words)//2], facet.words[-1]):
        moved = facet.incidence
        for step in word:
            moved = frozenset(permutations[step][position] for position in moved)
        assert moved == member


def test_packed_schreier_orbit_preserves_exact_words_and_membership():
    """Packing incidences changes neither group orbits nor Schreier transport words."""
    from collections import deque
    from itertools import combinations
    from sage_indefinite_port.indefinite.facet_orbits import _orbit

    degree = 10
    # Nine adjacent transpositions generate S_10 and satisfy the packed
    # path's many-generator criterion. The orbit has C(10,5)=252 facets.
    actions = []
    for index in range(degree - 1):
        permutation = list(range(degree))
        permutation[index], permutation[index+1] = permutation[index+1], permutation[index]
        actions.append(tuple(permutation))
    initial = frozenset((0, 2, 4, 6, 8))
    result = _orbit(initial, tuple(actions))
    expected = {initial: ()}
    pending = deque([initial])
    while pending:
        current = pending.popleft()
        for index, action in enumerate(actions):
            moved = frozenset(action[position] for position in current)
            if moved not in expected:
                expected[moved] = expected[current] + (index,)
                pending.append(moved)
    assert len(result) == len(list(combinations(range(degree), 5))) == 252
    assert result == expected


def test_large_ray_universe_packed_schreier_actions_are_exact():
    """Sparse-support permutations on 205 rays keep exact BFS words.

    The support has 10 active positions; 195 other rays are fixed and
    belong to every incidence. Generator-wise BFS must agree, even
    though the facet incidence has 102 positions and the ambient
    stabilizer acts on a 205-point set.
    """
    from collections import deque
    from sage_indefinite_port.indefinite.facet_orbits import _orbit

    degree = 205
    actions = []
    for index in range(9):
        permutation = list(range(degree))
        permutation[index], permutation[index + 1] = permutation[index + 1], permutation[index]
        actions.append(tuple(permutation))
    initial = frozenset((*range(5), *range(100, 197)))
    packed = _orbit(initial, tuple(actions))
    baseline = {initial: ()}
    pending = deque([initial])
    while pending:
        current = pending.popleft()
        for generator_index, permutation in enumerate(actions):
            moved = frozenset(permutation[position] for position in current)
            if moved not in baseline:
                baseline[moved] = baseline[current] + (generator_index,)
                pending.append(moved)
    assert len(packed) == 252
    assert packed == baseline


def test_frozen_rank_seventeen_second_large_cell_all_facet_orbits():
    r"""The source's 205-ray cell has exactly 65,533 facets in 15 orbits.

    Construct this cell by actual certified perfect-cell flips from the
    frozen I_(16,1) initial cell, rather than copying its rows into a test
    fixture. Its fourteen exchangeable coordinate positions generate
    S_14 inside the actual cell stabilizer. The exact 57-ray
    fundamental-chamber dual description, filtered by facet-rank
    certificates, is therefore complete under the original arithmetic
    group. All 15 orbit counts are certified by the supplied ray
    permutations. Independent general-orbit BFS in the source order
    gives these same sizes on the frozen cell.
    """
    from sage_indefinite_port.indefinite.lorentzian_cells import _facet_incidence
    from sage_indefinite_port.indefinite.facet_orbits import _facet_orbit_records

    ambient = Lattices(ZZ)([
        [int(i == j) * (1 if i == 0 else -1) for j in range(17)]
        for i in range(17)
    ])
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(ambient, NONNEGATIVE_SQUARES)
    for target_count in (17, 18, 205):
        _permutations, facet_orbits = backend._facet_mask_orbits_via_walk(cell)
        next_cell = None
        for orbit in facet_orbits:
            facet = _facet_incidence(orbit[0][0], len(cell.rows))
            neighbor = backend.flip_across(cell, facet)
            if neighbor is not None and len(neighbor.rows) == target_count:
                next_cell = neighbor
                break
        assert next_cell is not None, f"the exact cell walk lost its {target_count}-ray neighbor"
        cell = next_cell
    assert len(cell.rows) == 205

    # Complete source generator family before normalization. Every
    # discarded generator must be a swap whose endpoints are connected
    # by retained exchange generators, so group equality is proved
    # without a costly black-box permutation-group order comparison.
    original = tuple(backend.cell_transporter(
        cell, cell, all_images=True, generators_only=True,
    ))
    complete = tuple(backend.cell_stabilizer(cell))
    assert len(original) >= len(complete)
    assert len(complete) == 16
    kept = set(complete)
    labels = tuple(ambient.module_generating_set())
    parent = list(range(17))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for g in complete:
        matrix_rows = matrix(SageZZ, [
            [int(g(ambient.module_generator(label)).to_vector()(target))
             for target in labels] for label in labels
        ])
        moved = [i for i in range(17)
                 if any(matrix_rows[i, j] != int(i == j) for j in range(17))]
        if len(moved) == 2 and all(
            matrix_rows[moved[i], j] == int(j == moved[1-i])
            for i in range(2) for j in range(17)
        ):
            parent[root(moved[0])] = root(moved[1])
    for g in original:
        if g in kept:
            continue
        matrix_rows = matrix(SageZZ, [
            [int(g(ambient.module_generator(label)).to_vector()(target))
             for target in labels] for label in labels
        ])
        moved = [i for i in range(17)
                 if any(matrix_rows[i, j] != int(i == j) for j in range(17))]
        assert len(moved) == 2 and all(
            matrix_rows[moved[i], j] == int(j == moved[1-i])
            for i in range(2) for j in range(17)
        )
        assert root(moved[0]) == root(moved[1])

    permutations = backend.cell_stabilizer_permutations(cell)
    assert len(permutations) == len(complete)
    # A later standalone exhaustive stabilizer query must not overwrite
    # the cached subgroup's (reduced) row-action certificate.
    backend.cell_transporter(cell, cell, all_images=True, generators_only=True)
    assert len(backend._stabilizer_actions[cell]) == len(complete)
    assert backend.cell_stabilizer_permutations(cell) == permutations
    records = tuple(_facet_orbit_records(cell.rows, permutations, retain_words=True))
    assert len(records) == 15
    assert sum(record.orbit_size for record in records) == 65533
    assert sorted((len(record.incidence), record.orbit_size) for record in records) == sorted([
        (17,4),(29,56),(39,364),(47,1456),(59,12012),
        (39,1456),(57,8008),(53,4004),(59,13728),(53,8008),
        (17,56),(57,12012),(191,1),(47,4004),(29,364),
    ])
    assert all(len(record.words) == record.orbit_size for record in records)


def test_marked_null_orbits_are_identical_with_or_without_full_facet_transversals():
    r"""Null-mark gluing uses complete quotient adjacencies, not every facet word.

    The two traversal representations retain identical exact cell and
    adjacency data; only the richer public traversal materializes the
    per-facet Schreier transversals. The null-mark algorithm must give
    exactly the same partition and vectors under either representation.
    """
    from sage_indefinite_port.indefinite.lorentzian_cells import (
        MarkedCellOrbitAlgorithm, perfect_complex,
    )

    ambient = Lattices(ZZ)([[0, 1, 0], [1, 0, 0], [0, 0, 2]])
    labels = tuple(ambient.module_generating_set())
    results = []
    for retain in (True, False):
        complex_ = perfect_complex(ambient, NONNEGATIVE_SQUARES, retain_facet_orbits=retain)
        orbits = MarkedCellOrbitAlgorithm(complex_, 0).global_orbits()
        assert len(complex_.quotient_cells()) == 1
        assert len(complex_.adjacencies()) == 1
        assert bool(complex_._traversal_facet_orbits) is retain
        results.append(tuple(tuple(
            tuple(int(mark.to_vector()(label)) for label in labels)
            for mark in orbit
        ) for orbit in orbits))
    assert results[0] == results[1]
    assert len(results[0]) == 1


def test_rank_seventeen_primal_extreme_rays_do_not_construct_exponential_dual(monkeypatch):
    """The 33-ray crosspolytope has 32 primitive extreme rays and 2^16 facets."""
    ambient = Lattices(ZZ)([
        [int(i == j) * (1 if i == 0 else -1) for j in range(17)]
        for i in range(17)
    ])
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(ambient, NONNEGATIVE_SQUARES)
    assert len(cell.rows) == 33
    monkeypatch.setattr(
        backend, "cell_cone",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("primitive null rays must not build the dual PPL cone")
        ),
    )
    rays = tuple(backend.cell_primitive_rays(cell, square=0))
    assert len(rays) == 32
    assert len(set(rays)) == 32
    assert all(int(ray.q()) == 0 and int(ray.content()) == 1 for ray in rays)
