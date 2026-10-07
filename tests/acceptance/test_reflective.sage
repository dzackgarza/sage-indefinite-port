"""Fundamental domains of reflective Lorentzian lattices (Allcock's edgewalk).

Simple roots are compared up to isometry: two root systems agree when their Gram matrices
agree up to a simultaneous permutation of the roots. Recorded data: polyhedral_common's
8,821 reflective forms (root counts) and 103 explicit root systems, Allcock's 8,595
rank-3 lattices with explicit roots, OSCAR's Vinberg tests, and Dutour Sikirić--Hulek's
ten simple roots of U + E8(-1).
"""

import pytest
from dzack_research.preamble.all import ZZ as PreambleZZ
from dzack_research.preamble.all import HyperbolicLattices

from sage_indefinite_port.readiness import UnfinishedCapability
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
    return norms, [len(cell) for cell in partition], graph.canonical_label(partition=partition, edge_labels=True).copy(immutable=True)


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


@pytest.mark.parametrize("case", load_reflective_forms(), ids=lambda case: case["id"])
@pytest.mark.xfail(reason="Allcock edgewalk: owned by #31", raises=UnfinishedCapability, strict=True)
def test_simple_root_counts_match_polyhedral_common(case) -> None:
    require(EDGEWALK)

    report = _hyperbolic(case["gram"]).allcock_edgewalk()

    assert report.is_reflective()
    assert len(report.simple_roots()) == case["num_simple_roots"]


@pytest.mark.parametrize("case", load_root_systems(), ids=lambda case: case["id"])
@pytest.mark.xfail(reason="Allcock edgewalk: owned by #31", raises=UnfinishedCapability, strict=True)
def test_root_systems_match_polyhedral_common(case) -> None:
    require(EDGEWALK)

    _assert_same_root_system(case["gram"], case["roots"])


@pytest.mark.parametrize("case", load_allcock_rank3_reflective(), ids=lambda case: case["id"])
@pytest.mark.xfail(reason="Allcock edgewalk: owned by #31", raises=UnfinishedCapability, strict=True)
def test_root_systems_match_allcock(case) -> None:
    require(EDGEWALK)

    _assert_same_root_system(case["gram"], case["simple_roots"])


_VINBERG = {case["id"]: case for case in load_oscar_lattice_oracles()["vinberg"]}


@pytest.mark.parametrize("case", list(_VINBERG.values()), ids=lambda case: case["id"])
@pytest.mark.xfail(reason="Allcock edgewalk: owned by #31", raises=UnfinishedCapability, strict=True)
def test_simple_root_counts_match_oscar_vinberg(case) -> None:
    require(EDGEWALK)

    assert len(_computed_roots(case["gram"])) == case["num_simple_roots"]


@pytest.mark.xfail(reason="Allcock edgewalk: owned by #31", raises=UnfinishedCapability, strict=True)
def test_vinberg_coxeter_diagram_has_oscars_edges() -> None:
    require(EDGEWALK)

    case = _VINBERG["vinberg_scharlau_U_A2"]
    roots = _computed_roots(case["gram"])
    root_gram = _root_gram(case["gram"], roots)

    edges = sum(1 for i in range(len(roots)) for j in range(i + 1, len(roots)) if root_gram[i][j])

    assert edges == case["coxeter_diagram_edges"]


@pytest.mark.xfail(reason="Allcock edgewalk: owned by #31", raises=UnfinishedCapability, strict=True)
def test_vinberg_roots_match_oscars_chamber_up_to_isometry() -> None:
    require(EDGEWALK)

    case = _VINBERG["vinberg_diag_1_-1_-3"]

    _assert_same_root_system(case["gram"], case["simple_roots_of_test_chamber"])


@pytest.mark.xfail(reason="Allcock edgewalk: owned by #31", raises=UnfinishedCapability, strict=True)
def test_e10_simple_roots_match_dutour_sikiric_hulek() -> None:
    require(EDGEWALK)

    domain = load_e10_fundamental_domain()

    _assert_same_root_system(domain["gram"], [root["vector"] for root in domain["simple_roots"]])
