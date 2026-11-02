"""The facet-route certificates, checked against the engine's own facet enumeration.

Ground truth is independent of the module under test in both directions used here:

* for the port's own cells, the facets come from `cell_cone`, i.e. from Normaliz's exact H-description, and
  the certificate layer has to recover each facet's primitive normal **from its incidence set alone** and then
  certify it;
* for the cube cones the facet family is closed form (`x_0 +- x_j >= 0`, `2(n-1)` facets) and adjacency is
  analytic: all pairs are adjacent except the `n-1` opposite pairs, so `C(2(n-1), 2) - (n-1)` pairs.

The structural fact both rank tests lean on is checked here too: `LorentzianPerfectCell.__post_init__` enforces
`rank[1|A] = n`, which puts the rows in one affine hyperplane not through the origin, so the cone is a cone over
an `(n-1)`-polytope and is pointed and full-dimensional a priori.
"""
import pytest
from sage.all import GF, QQ, ZZ, matrix

from sage_indefinite_port.indefinite.facet_certificates import (
    adjacent_incidences,
    facet_normal_from_incidence,
    incidence_rank_certificate,
    is_facet_certificate,
    orient_facet_normal,
    pointedness_functional,
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


def _normal_coordinates(covector, lattice_):
    return tuple(
        int(QQ(int(covector(lattice_.module_generator(label)).numerator()))
            / QQ(int(covector(lattice_.module_generator(label)).denominator())))
        for label in lattice_.module_generating_set()
    )


def test_rank_certificate_is_a_nonzero_minor_and_never_over_claims() -> None:
    """The certificate is a determinant Sage evaluates, and it must agree with exact `Q`-rank.

    This is a check of the primitive on chosen cases. It is emphatically not evidence that a facet enumeration
    is complete: a missed rank would silently drop a genuine facet or ridge, and no sample of examples can
    certify an enumeration. What can is the closure argument over the facet graph, which the walk owns.
    """
    cases = [
        [[1, 2, 3], [2, 4, 6], [0, 1, 1]],          # dependent row included, rank 2
        [[1, 0], [0, 1], [1, 1]],                    # rank 2 in dimension 2
        [[2, 4], [3, 6]],                            # rank 1
        [[0, 0], [0, 0]],                            # rank 0
    ]
    for rows in cases:
        exact = matrix(QQ, rows).rank()
        certificate = incidence_rank_certificate(rows, exact) if exact else None
        if exact:
            assert certificate is not None, f"rank {exact} must be certifiable for {rows}"
            pivot_rows, pivot_columns = certificate
            minor = matrix(ZZ, rows).matrix_from_rows_and_columns(pivot_rows, pivot_columns)
            assert minor.det() != 0, "the returned certificate must be a nonzero minor"
            assert len(pivot_rows) == exact
        assert incidence_rank_certificate(rows, exact + 1) is None, "the certificate must not over-claim"
    assert incidence_rank_certificate([[0, 0], [0, 0]], 1) is None


@pytest.mark.parametrize("name", sorted(SPECIMENS))
def test_every_engine_facet_normal_is_recovered_from_its_incidence(name) -> None:
    """Two independent routes to the same normal: Normaliz's H-description, and the incidence's cofactors."""
    backend, cell = _cell(SPECIMENS[name])
    dimension = len(cell.rows[0])
    rays = [tuple(int(entry) for entry in row) for row in cell.rows]

    augmented = matrix(QQ, [[1, *row] for row in rays])
    assert augmented.rank() == dimension, "the cell's rows must lie on one affine hyperplane of rank n-1"

    covectors = backend.cell_cone(cell).halfspace_covectors()
    assert covectors
    for covector in covectors:
        coordinates = _normal_coordinates(covector, cell.lattice)
        incidence = tuple(
            index
            for index, row in enumerate(rays)
            if sum(coordinates[i] * row[i] for i in range(dimension)) == 0
        )
        raw = facet_normal_from_incidence([rays[index] for index in incidence], dimension)
        # The kernel basis sign is a coordinate fact; the cone's orientation is the one non-negative on all rays.
        recovered, _values = orient_facet_normal(raw, rays)
        assert is_facet_certificate(recovered, incidence, rays), (
            f"the recovered normal {recovered} failed its own certificate on incidence {incidence}"
        )
        assert recovered == tuple(coordinates) or recovered == tuple(-entry for entry in coordinates), (
            "the incidence alone must determine the engine's facet normal up to sign"
        )


@pytest.mark.parametrize("name", sorted(SPECIMENS))
def test_adjacency_certificate_matches_the_ridge_rank_criterion(name) -> None:
    """Adjacency is decided from incidences; cross-check it against the exact `QQ` rank it certifies."""
    backend, cell = _cell(SPECIMENS[name])
    dimension = len(cell.rows[0])
    rays = [tuple(int(entry) for entry in row) for row in cell.rows]
    facets = []
    for covector in backend.cell_cone(cell).halfspace_covectors():
        coordinates = _normal_coordinates(covector, cell.lattice)
        facets.append(
            frozenset(
                index
                for index, row in enumerate(rays)
                if sum(coordinates[i] * row[i] for i in range(dimension)) == 0
            )
        )
    for first in range(len(facets)):
        for second in range(first + 1, len(facets)):
            joint = sorted(facets[first] & facets[second])
            exact = matrix(QQ, [rays[index] for index in joint]).rank() if joint else 0
            certified = adjacent_incidences(facets[first], facets[second], rays, dimension)
            assert certified == (exact == dimension - 2), (
                f"adjacency disagreement on facets {first}, {second}: certificate {certified}, rank {exact}"
            )


def test_cube_cone_adjacency_matches_the_closed_form() -> None:
    """A closed-form specimen: `2(n-1)` facets, every pair adjacent except the `n-1` opposite pairs."""
    for n in (4, 6):
        rays = [tuple(int(entry) for entry in row) for row in
                __import__("itertools").product((-1, 1), repeat=n - 1)]
        rays = [(*row, 1) for row in rays]
        facets = []
        for coordinate in range(n - 1):
            for sign in (-1, 1):
                # facet x_coordinate = sign  <=>  1 - sign * x_coordinate = 0, with the constant last
                normal = [0] * n
                normal[coordinate] = -sign
                normal[n - 1] = 1
                incidence = frozenset(
                    index
                    for index, row in enumerate(rays)
                    if sum(normal[i] * row[i] for i in range(n)) == 0
                )
                facets.append((tuple(normal), incidence))
        assert len(facets) == 2 * (n - 1)
        for normal, incidence in facets:
            assert is_facet_certificate(normal, incidence, rays)
        adjacent = sum(
            1
            for first in range(len(facets))
            for second in range(first + 1, len(facets))
            if adjacent_incidences(facets[first][1], facets[second][1], rays, n)
        )
        expected = (2 * (n - 1)) * (2 * (n - 1) - 1) // 2 - (n - 1)
        assert adjacent == expected, f"n={n}: {adjacent} adjacent pairs, closed form says {expected}"


@pytest.mark.parametrize("name", sorted(SPECIMENS))
def test_pointedness_is_certified_by_the_constructor_hypothesis(name) -> None:
    r"""Pointedness from `rank[1|A] = n`, checked row by row on the actual cell.

    The cell's rows lie in one affine hyperplane through no origin; the functional `w'` from the kernel of
    `[1|A]` therefore takes the *same nonzero value* on every ray, so the cone is the cone over an
    `(n-1)`-polytope and is pointed and full-dimensional. This is not an inference from `rank A = n` (which is
    only full-dimensionality) and not an extrapolation from small samples: it is a proof for every cell the port
    constructs, evaluated on this cell's own rows.
    """
    backend, cell = _cell(SPECIMENS[name])
    dimension = len(cell.rows[0])
    functional, pairings = pointedness_functional(cell)

    assert all(value > 0 for value in pairings), "the affine functional must be strictly positive on every ray"
    assert len(set(pairings)) == 1, (
        "the rows lie in an affine hyperplane, so the functional takes one constant value on them, "
        f"measured {sorted(set(pairings))[:4]}"
    )
    assert matrix(QQ, [[1, *row] for row in cell.rows]).rank() == dimension


def test_orient_refuses_a_normal_that_annihilates_every_ray() -> None:
    """The earlier version returned all-zero pairings as if they were a valid orientation; it must raise."""
    rays = [(1, 0), (0, 1)]
    with pytest.raises(ArithmeticError):
        orient_facet_normal((0, 0), rays)
    with pytest.raises(ArithmeticError):
        orient_facet_normal((1, -1), rays)


def test_a_nonzero_minor_modulo_a_prime_certifies_a_rational_lower_bound() -> None:
    r"""Correcting an earlier claim of mine: non-vanishing mod `p` *does* prove `rank >= k` over `Q`.

    `det != 0` in `Z/pZ` implies `det != 0` in `Z`, hence in `Q`. What modular arithmetic cannot certify is an
    upper bound, which is why the upper bounds in this module are structural rather than computed.
    """
    prime = 101
    rows = [[1, 2], [2, 4], [3, 7]]             # rank 2 over Q; the 2x2 minor rows 0,2 is 1*7-2*3 = 1
    minor = matrix(GF(prime), [[rows[0][0], rows[0][1]], [rows[2][0], rows[2][1]]])
    assert minor.det() != 0, "a nonzero modular determinant is a rational rank certificate"
    assert matrix(QQ, rows).rank() >= 2
