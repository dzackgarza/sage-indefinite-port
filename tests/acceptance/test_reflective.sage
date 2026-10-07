"""Fundamental domains of reflective Lorentzian lattices (Allcock's edgewalk).

Simple roots are compared up to isometry: two root systems agree when their Gram matrices
agree up to a simultaneous permutation of the roots. Recorded data: polyhedral_common's
8,821 reflective forms (root counts) and 103 explicit root systems, Allcock's 8,595
rank-3 lattices with explicit roots, OSCAR's Vinberg tests, and Dutour Sikirić--Hulek's
ten simple roots of U + E8(-1).
"""

import pytest

from tests.acceptance.consumer import coordinates, lattice
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
    return norms, [len(cell) for cell in partition], graph.canonical_label(partition=partition, edge_labels=True).copy(immutable=True)


def _computed_roots(gram):
    return [coordinates(root) for root in lattice(gram).allcock_edgewalk().simple_roots()]


def _assert_same_root_system(gram, recorded_roots) -> None:
    computed = _computed_roots(gram)
    assert len(computed) == len(recorded_roots)
    assert _canonical_root_graph(_root_gram(gram, computed)) == _canonical_root_graph(_root_gram(gram, recorded_roots))


@pytest.mark.parametrize("case", load_reflective_forms(), ids=lambda case: case["id"])
def test_simple_root_counts_match_polyhedral_common(case) -> None:
    report = lattice(case["gram"]).allcock_edgewalk()

    assert report.is_reflective()
    assert len(report.simple_roots()) == case["num_simple_roots"]


@pytest.mark.parametrize("case", load_root_systems(), ids=lambda case: case["id"])
def test_root_systems_match_polyhedral_common(case) -> None:
    _assert_same_root_system(case["gram"], case["roots"])


@pytest.mark.parametrize("case", load_allcock_rank3_reflective(), ids=lambda case: case["id"])
def test_root_systems_match_allcock(case) -> None:
    _assert_same_root_system(case["gram"], case["simple_roots"])


_VINBERG = {case["id"]: case for case in load_oscar_lattice_oracles()["vinberg"]}


@pytest.mark.parametrize("case", list(_VINBERG.values()), ids=lambda case: case["id"])
def test_simple_root_counts_match_oscar_vinberg(case) -> None:
    assert len(_computed_roots(case["gram"])) == case["num_simple_roots"]


def test_vinberg_coxeter_diagram_has_oscars_edges() -> None:
    case = _VINBERG["vinberg_scharlau_U_A2"]
    roots = _computed_roots(case["gram"])
    root_gram = _root_gram(case["gram"], roots)

    edges = sum(1 for i in range(len(roots)) for j in range(i + 1, len(roots)) if root_gram[i][j])

    assert edges == case["coxeter_diagram_edges"]


def test_vinberg_roots_match_oscars_chamber_up_to_isometry() -> None:
    case = _VINBERG["vinberg_diag_1_-1_-3"]

    _assert_same_root_system(case["gram"], case["simple_roots_of_test_chamber"])


def test_e10_simple_roots_match_dutour_sikiric_hulek() -> None:
    domain = load_e10_fundamental_domain()

    _assert_same_root_system(domain["gram"], [root["vector"] for root in domain["simple_roots"]])
