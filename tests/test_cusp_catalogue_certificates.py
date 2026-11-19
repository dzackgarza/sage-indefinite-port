"""Independent arithmetic certificates for the frozen rank-13 and rank-17 cusp families.

The tests verify actual representatives, pairwise inequivalence invariants,
all possible local definite genera, and the exact Smith--Minkowski--Siegel
mass. They do not use a failed/timeout search as a classification verdict.
"""
from sage.libs.pari import pari
from sage.matrix.constructor import matrix
from sage.matrix.special import diagonal_matrix
from sage.quadratic_forms.genera.genus import Genus, genera
from sage.quadratic_forms.quadratic_form import QuadraticForm
from sage.rings.integer_ring import ZZ as SageZZ
from sage.rings.rational_field import QQ as SageQQ
from sage_extensions.all import Lattices, ZZ

from sage_indefinite_port.indefinite.lorentzian_cells import _gram_matrix


def _mark(lattice, coordinates):
    labels = tuple(lattice.module_generating_set())
    result = lattice.linear_combination({
        labels[index]: value for index, value in coordinates.items()
    })
    assert result.q() == 0 and result.content() == 1
    return result


def _positive_quotient(mark):
    quotient = mark.isotropic_reduction().quotient_lattice()
    signature = tuple(int(value) for value in quotient.signature_pair())
    assert signature == (0, int(quotient.rank()))
    return -_gram_matrix(quotient)


def test_five_distinct_rank_seventeen_unimodular_cusps():
    r"""Conway--Sloane (1982): all positive-definite unimodular rank-15
    lattices form five isometry classes. Every primitive isotropic vector
    in the odd indefinite unimodular rank-17 lattice splits a U-plane
    and hence determines one such rank-15 lattice. These five vectors
    have pairwise different numbers of norm-one vectors, so they are
    inequivalent and exhaust the possible orbit classes.
    """
    gram = diagonal_matrix(SageZZ, [1] + [-1] * 16)
    lattice = Lattices(ZZ)([list(row) for row in gram.rows()])
    candidates = [
        {0: 1, 16: -1},
        {0: 4, **{i: -1 for i in range(1, 17)}},
        {0: 4, **{i: -1 for i in range(4, 16)}, 16: -2},
        {0: 3, **{i: -1 for i in (6, 7, 8, 9, 10, 11, 13, 14, 15)}},
        {
            0: 6, 1: 1,
            **{i: -1 for i in (3, 4, 5, 6, 7, 8, 11)},
            **{i: -2 for i in (9, 10, 12, 13, 14, 15, 16)},
        },
    ]
    invariants = []
    for coordinates in candidates:
        positive = _positive_quotient(_mark(lattice, coordinates))
        assert positive.nrows() == 15 and positive.det() == 1
        invariants.append(int(pari(positive).qfminim(1)[0]))
    assert sorted(invariants) == [0, 2, 6, 14, 30]


def test_four_rank_thirteen_cusps_exhaust_the_allowed_definite_genus_mass():
    r"""Rank-11/determinant-3 definite forms have exactly two genera.
    Only the genus of Z^9+A2 can be an orthogonal complement to a
    primitive U in the genus of diag(3,-1^12). These four primitive
    null vectors lie in four distinct quotient classes, and the sum
    of their reciprocal automorphism orders equals the entire mass
    of that allowed genus. Hence no additional null-vector orbit exists.
    """
    indefinite = diagonal_matrix(SageZZ, [3] + [-1] * 12)
    lattice = Lattices(ZZ)([list(row) for row in indefinite.rows()])
    a2 = diagonal_matrix(SageZZ, [1] * 9 + [2, 2])
    a2[9, 10] = a2[10, 9] = 1
    other = diagonal_matrix(SageZZ, [1] * 10 + [3])
    possible = tuple(genera([11, 0], 3))
    assert len(possible) == 2
    assert Genus(a2) != Genus(other)
    assert all(any(Genus(reference) == genus for genus in possible)
               for reference in (a2, other))

    hyperbolic = matrix(SageZZ, [[0, 1], [1, 0]])

    def add_negative(positive):
        return matrix(
            SageZZ, 13, 13,
            lambda i, j: (
                hyperbolic[i, j] if i < 2 and j < 2
                else -positive[i - 2, j - 2] if i >= 2 and j >= 2
                else 0
            ),
        )

    assert Genus(add_negative(a2)) == Genus(indefinite)
    assert Genus(add_negative(other)) != Genus(indefinite)

    candidates = [
        {0: 4, 1: -6, 2: -3, 10: -1, 11: -1, 12: -1},
        {0: 5, 1: -7, 2: -4, **{i: -1 for i in range(3, 13)}},
        {0: 4, 1: -6, 2: -2, **{i: -1 for i in (4, 5, 7, 8, 9, 10, 11, 12)}},
        {0: 3, 1: -3, 2: -3, **{i: -1 for i in range(4, 13)}},
    ]
    expected_orders = [2229534720, 958003200, 278691840, 16721510400]
    expected_theta = [
        (18, 168, 948), (0, 132, 1056), (6, 144, 1020), (2, 248, 740),
    ]
    reciprocals = SageQQ(0)
    for coordinates, expected_order, expected_counts in zip(
        candidates, expected_orders, expected_theta, strict=True,
    ):
        positive = _positive_quotient(_mark(lattice, coordinates))
        assert positive.nrows() == 11 and positive.det() == 3
        assert Genus(positive) == Genus(a2)
        counts = tuple(int(pari(positive).qfminim(bound)[0]) for bound in (1, 2, 3))
        assert counts == expected_counts
        order = int(pari(positive).qfauto()[0])
        assert order == expected_order
        reciprocals += SageQQ(1) / order
    mass = QuadraticForm(SageZZ, 2 * a2).conway_mass()
    assert mass == SageQQ(1891) / 367873228800
    assert reciprocals == mass
