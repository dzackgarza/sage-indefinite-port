"""Fundamental domains of reflective Lorentzian lattices (Allcock's edgewalk).

Simple roots are compared up to isometry: two root systems agree when their Gram matrices
agree up to a simultaneous permutation of the roots. Recorded data: polyhedral_common's
8,821 reflective forms (root counts) and 103 explicit root systems, Allcock's 8,595
rank-3 lattices with explicit roots, OSCAR's Vinberg tests, and Dutour Sikirić--Hulek's
ten simple roots of U + E8(-1).
"""

import pytest
from collections import Counter
from sage_extensions.all import ZZ as PreambleZZ
from sage_extensions.all import HyperbolicLattices
from sage.all import QQ, ZZ, Graph, matrix, vector

from tests.acceptance.consumer import EDGEWALK, coordinates, lattice, require
from tests.fixtures.oracle_fixtures import (
    load_allcock_rank3_reflective,
    load_e10_fundamental_domain,
    load_oscar_lattice_oracles,
    load_reflective_forms,
    load_root_systems,
)


def _root_gram(gram, roots):
    form = matrix(ZZ, gram)
    vectors = [vector(ZZ, root) for root in roots]
    return [[int(left * form * right) for right in vectors] for left in vectors]


def _canonical_root_graph(root_gram):
    """The labelled root graph up to simultaneous permutation of the roots."""
    size = len(root_gram)
    graph = Graph([list(range(size)), [(i, j, root_gram[i][j]) for i in range(size) for j in range(i + 1, size) if root_gram[i][j]]], format="vertices_and_edges")
    norms = sorted({root_gram[i][i] for i in range(size)})
    partition = [[i for i in range(size) if root_gram[i][i] == norm] for norm in norms]
    return norms, [len(cell) for cell in partition], graph.canonical_label(partition=partition, edge_labels=True, algorithm="sage").copy(immutable=True)


def _ambient_root_transport(gram, source_roots, target_roots):
    r"""Return an actual integral lattice automorphism carrying one root set to the other.

    Source and target roots come from distinct implementations. Equality of
    their Gram graphs is not sufficient for an ambient integral isometry:
    the roots may span a proper finite-index sublattice, so the overlattice
    gluing must also be preserved. For a full-rank root configuration, any
    ambient isometry is uniquely determined by its action on *any* rational
    basis chosen from the source roots. Enumerate the finitely many candidate
    images of that basis, preserving all signed pairings and each root's
    full pairing multiset; solve for the unique rational matrix and require
    that it be integral, unimodular, preserve the original Gram, and map the
    ENTIRE unordered source root configuration onto the target. The latter
    condition proves this is an actual isometry of the original lattice,
    not just a rational isometry of two root spans with the same Coxeter
    graph. The search is exhaustive for full-rank configurations because
    every ambient root-set isometry restricts to one of these basis maps.
    """
    original_gram = matrix(ZZ, gram)
    source = matrix(ZZ, source_roots)
    target = matrix(ZZ, target_roots)
    rank = original_gram.nrows()
    if source.nrows() != target.nrows() or source.ncols() != rank or source.rank() != rank or target.rank() != rank:
        return None

    pairing_source = source * original_gram * source.transpose()
    pairing_target = target * original_gram * target.transpose()
    source_colors = [
        (pairing_source[i, i], tuple(sorted(pairing_source.row(i))))
        for i in range(source.nrows())
    ]
    target_colors = [
        (pairing_target[j, j], tuple(sorted(pairing_target.row(j))))
        for j in range(target.nrows())
    ]
    if Counter(source_colors) != Counter(target_colors):
        return None

    independent = []
    for i in range(source.nrows()):
        if source.matrix_from_rows(independent + [i]).rank() > len(independent):
            independent.append(i)
            if len(independent) == rank:
                break
    basis = source.matrix_from_rows(independent).change_ring(QQ)
    candidate_positions = [
        tuple(j for j, color in enumerate(target_colors) if color == source_colors[i])
        for i in independent
    ]
    expected_vectors = Counter(tuple(row) for row in target.rows())

    def extend(images):
        depth = len(images)
        if depth == rank:
            candidate_basis = target.matrix_from_rows(images).change_ring(QQ)
            transport = basis.inverse() * candidate_basis
            if (all(entry in ZZ for entry in transport.list())
                    and abs(transport.det()) == 1
                    and transport * original_gram * transport.transpose() == original_gram
                    and Counter(tuple(row) for row in source * transport) == expected_vectors):
                return transport.change_ring(ZZ)
            return None
        i = independent[depth]
        for j in candidate_positions[depth]:
            if j in images:
                continue
            if any(pairing_source[i, independent[k]] != pairing_target[j, images[k]]
                   for k in range(depth)):
                continue
            result = extend(images + [j])
            if result is not None:
                return result
        return None

    return extend([])


def _hyperbolic(gram):
    """The lattice in HyperbolicLattices' signature (n, 1); a (1, n) Gram matrix is negated, which keeps coordinates."""
    negative = sum(1 for value in matrix(QQ, gram).eigenvalues() if value < 0)
    form = matrix(ZZ, gram)
    assert negative == 1 or negative == form.nrows() - 1, f"{gram} is not hyperbolic"
    oriented = form if negative == 1 else -form
    return HyperbolicLattices(PreambleZZ)(lattice([list(row) for row in oriented.rows()]))


def _computed_roots(gram):
    return [coordinates(root) for root in _hyperbolic(gram).allcock_edgewalk().simple_roots()]


def _assert_same_root_system(gram, recorded_roots) -> None:
    computed = _computed_roots(gram)
    assert len(computed) == len(recorded_roots)
    assert _canonical_root_graph(_root_gram(gram, computed)) == _canonical_root_graph(_root_gram(gram, recorded_roots))
    # Source data records roots in the ORIGINAL integral lattice. Its
    # acceptance contract is ambient integral equivalence, not equality of
    # abstract Coxeter diagrams (which loses overlattice gluing). The
    # certified matrix proves the stronger claim for full-rank root sets.
    transporter = _ambient_root_transport(gram, recorded_roots, computed)
    assert transporter is not None, "matching root Gram graphs need not extend to an integral ambient isometry"


@pytest.mark.parametrize("case", load_reflective_forms(), ids=lambda case: case["id"])
def test_simple_root_counts_match_polyhedral_common(case) -> None:
    require(EDGEWALK)

    report = _hyperbolic(case["gram"]).allcock_edgewalk()

    assert report.is_reflective()
    assert len(report.simple_roots()) == case["num_simple_roots"]


@pytest.mark.parametrize("case", load_root_systems(), ids=lambda case: case["id"])
def test_root_systems_match_polyhedral_common(case) -> None:
    require(EDGEWALK)

    _assert_same_root_system(case["gram"], case["roots"])


def _allcock_case_id(case) -> str:
    r"""The retained source identity of a rank-3 Allcock lattice.

    ``case["id"]`` is built from the lattice id alone
    (``references/extract/allcock_rank3.py:56`` uses only capture group 5), so it names up
    to 374 distinct lattices: 8,595 entries share 176 ids. ``(weyl_group_id, lattice_id)``
    is Allcock's own identity and is unique across every entry, so composing it keeps the
    node id both readable and unambiguous. Nothing but the reported id changes -- the
    parametrized data, its order and every assertion are untouched.
    """
    return f"{case['id']}_w{case['weyl_group_id']}_l{case['lattice_id']}"


@pytest.mark.parametrize("case", load_allcock_rank3_reflective(), ids=_allcock_case_id)
def test_root_systems_match_allcock(case) -> None:
    require(EDGEWALK)

    _assert_same_root_system(case["gram"], case["simple_roots"])


_VINBERG = {case["id"]: case for case in load_oscar_lattice_oracles()["vinberg"]}


@pytest.mark.parametrize("case", list(_VINBERG.values()), ids=lambda case: case["id"])
def test_simple_root_counts_match_oscar_vinberg(case) -> None:
    require(EDGEWALK)

    assert len(_computed_roots(case["gram"])) == case["num_simple_roots"]


def test_vinberg_coxeter_diagram_has_oscars_edges() -> None:
    require(EDGEWALK)

    case = _VINBERG["vinberg_scharlau_U_A2"]
    roots = _computed_roots(case["gram"])
    root_gram = _root_gram(case["gram"], roots)

    # Oscar.jl/test/NumberTheory/vinberg.jl:43-47 counts intersecting walls,
    # including orthogonal ones, by k^2 / (q_i q_j) <= 1.
    edges = sum(
        1
        for i in range(len(roots))
        for j in range(i + 1, len(roots))
        if root_gram[i][j] ** 2 <= root_gram[i][i] * root_gram[j][j]
    )

    assert edges == case["coxeter_diagram_edges"]


def test_vinberg_roots_match_oscars_chamber_up_to_isometry() -> None:
    require(EDGEWALK)

    case = _VINBERG["vinberg_diag_1_-1_-3"]

    _assert_same_root_system(case["gram"], case["simple_roots_of_test_chamber"])


def test_e10_simple_roots_match_dutour_sikiric_hulek() -> None:
    require(EDGEWALK)

    domain = load_e10_fundamental_domain()

    _assert_same_root_system(domain["gram"], [root["vector"] for root in domain["simple_roots"]])
