"""`L = U + U`, `Gamma = Stab_O(L)(e1)` pointwise, and the family `v_m = e2 + m f1` (issues #53 and #62).

The obligation is narrow and easy to get wrong in one specific way: `Gamma` here is the **pointwise** stabilizer
of `e1`, not `O(L)` and not a finite-index surrogate. The `v_m` family exists to show the quotient `L[0]/Gamma` has
infinitely many classes -- not to replace the quotient by a parametrised slice.

Measured facts this file establishes, all through the retained APIs:

* `b(e1, v_m) = m` is a `Gamma`-invariant, because `gamma(e1) = e1` and `gamma` preserves the pairing. Distinct
  `m` therefore give distinct `Gamma`-orbits -- a proved negative, not a failed search.
* `Aut((L, Z e1))` is strictly larger than the pointwise stabilizer: the line group contains an element with
  `g(e1) = -e1`, the pointwise group contains none. Based vectors and unoriented lines are different quotients.
* a **non-identity** witnessed `Gamma`-transport exists for the family: the isometry swapping `e2` and `f2` (and
  fixing `e1`, `f1`) lies in `Gamma` and carries `v_m` to `f2 + m f1`, a different vector with the same invariant.

What is *not* established, and must not be inferred: that vectors sharing an invariant are equivalent. For fixed
`m` the orbit is not decided here beyond the witnesses actually produced, and having generators of `Gamma` says
nothing about deciding membership in it -- that is subgroup membership in infinite index, which this port does not
decide. `L[0]/Gamma` as a whole is therefore **not** supplied as a decidable quotient; what is supplied is the
exact action input (generators fixing `e1`), a proved separating invariant, and a certificate checker for
transport witnesses.
"""
import pytest
from sage.all import ZZ, matrix

from sage_indefinite_port.indefinite.recursive import EichlerRecursion, _isometry_from_row_matrix
from tests.acceptance.consumer import element, lattice

_GRAM = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]  # e1, f1, e2, f2 with b(e_i, f_i) = 1
_RECURSION = EichlerRecursion()
_LATTICE = None


def _uu():
    global _LATTICE
    if _LATTICE is None:
        _LATTICE = lattice(_GRAM)
    return _LATTICE


def _e1():
    return element(_uu(), (1, 0, 0, 0))


def _f1():
    return element(_uu(), (0, 1, 0, 0))


def _e2():
    return element(_uu(), (0, 0, 1, 0))


def _f2():
    return element(_uu(), (0, 0, 0, 1))


def _v(m):
    """`v_m = e2 + m f1`, primitive isotropic with `b(e1, v_m) = m`."""
    return _e2() + _uu().scalar_multiple(m, _f1())


def _coordinates(vector):
    values = vector.to_vector()
    return tuple(int(values(label)) for label in _uu().module_generating_set())


def _swap_e2_f2():
    """The isometry swapping `e2` and `f2`, fixing `e1` and `f1`: an element of `Gamma`."""
    return _isometry_from_row_matrix(_uu(), _uu(), matrix(ZZ, [[1, 0, 0, 0], [0, 1, 0, 0],
                                                               [0, 0, 0, 1], [0, 0, 1, 0]]))


def test_the_family_has_the_stated_invariant_and_is_primitive_isotropic() -> None:
    r"""`q(v_m) = 0`, `v_m` primitive, `b(e1, v_m) = m` -- checked on the specimen, not assumed."""
    L = _uu()
    for m in range(5):
        vector = _v(m)
        assert int(vector.q()) == 0, f"v_{m} must be isotropic"
        assert int(vector.content()) == 1, f"v_{m} must be primitive"
        assert int(L.b(_e1(), vector)) == m, f"b(e1, v_{m}) must be {m}"


def test_pointwise_stabilizer_generators_really_fix_e1() -> None:
    r"""The retained generators of `Aut((L, e1))` are checked one by one against `e1`."""
    gamma = _RECURSION.pointwise_stabilizer_generators(_e1())
    assert gamma, "the pointwise stabilizer of an isotropic vector must not be trivial"
    for generator in gamma:
        assert generator(_e1()) == _e1(), "every generator must fix e1 pointwise"


def test_line_stabilizer_is_strictly_larger_than_the_pointwise_one() -> None:
    r"""`Aut((L, Z e1))` contains `e1 |-> -e1`; the pointwise stabilizer contains no such element.

    This is the based-vector versus unoriented-line distinction: the two quotients `L[0]/Gamma` for
    `Gamma = Stab(e1)` and `Gamma = Stab(Z e1)` are different problems, and the port must not conflate them.
    """
    pointwise = _RECURSION.pointwise_stabilizer_generators(_e1())
    line = _RECURSION.line_stabilizer_generators(_e1())

    assert any(generator(_e1()) == -_e1() for generator in line), (
        "the line stabilizer must contain an element acting as -1 on the line"
    )
    assert not any(generator(_e1()) == -_e1() for generator in pointwise), (
        "the pointwise stabilizer must fix e1, so no generator may send it to -e1"
    )


def test_the_invariant_proves_distinct_gamma_orbits_for_distinct_m() -> None:
    r"""`b(e1, v_m) = m` separates the whole family: infinitely many `Gamma`-orbits in one shell.

    This is the executable form of the #53 statement. `pointwise_stabilizer_separates` returning `True` is a
    decision -- the two vectors cannot be in one orbit of the pointwise stabilizer -- not a search that gave up.
    """
    marked = (_e1(),)
    for m in range(4):
        for n in range(m + 1, 4):
            assert _RECURSION.pointwise_stabilizer_separates(marked, _v(m), _v(n)), (
                f"v_{m} and v_{n} must be separated by b(e1, .)"
            )
            # and the invariant the separation is read off is exactly the pairing
            assert _RECURSION.pointwise_stabilizer_invariant(marked, _v(m))[2] == m
            assert _RECURSION.pointwise_stabilizer_invariant(marked, _v(n))[2] == n


def test_a_nonidentity_gamma_transport_is_witnessed_and_certified() -> None:
    r"""The `e2 <-> f2` swap lies in `Gamma` and carries `v_m` to a different vector of the same invariant.

    This is the positive direction: the certificate checker accepts a real Gamma-element with the required action,
    and the image is checked to be neither the source nor separated from it by the invariant.
    """
    marked = (_e1(),)
    swap = _swap_e2_f2()

    for m in (0, 1, 3):
        source = _v(m)
        target = swap(source)
        assert target != source, "the swap must move the marking, otherwise the witness is trivial"
        assert int(_uu().b(_e1(), target)) == m, "the image must keep the invariant"
        assert not _RECURSION.pointwise_stabilizer_separates(marked, source, target), (
            "the invariant must not separate a vector from its image under a Gamma-element"
        )
        _RECURSION.verify_pointwise_stabilizer_transport(marked, source, target, swap)


def test_the_certificate_checker_refuses_a_non_gamma_element() -> None:
    r"""An isometry that does not fix `e1` is refused, so the checker is not vacuous.

    The reflection `e2 <-> -e2` does not fix `e1`? It does fix `e1`; what fails is the *transport* condition for
    the wrong target. Both failure modes are exercised: a candidate outside `Gamma`, and a candidate inside it that
    does not reach the claimed target.
    """
    marked = (_e1(),)
    swap = _swap_e2_f2()

    # a genuine Gamma-element, but asked to transport onto something it does not reach
    with pytest.raises(ArithmeticError):
        _RECURSION.verify_pointwise_stabilizer_transport(marked, _v(1), _v(1) + _e2(), swap)

    # an isometry that does not fix e1: the U(e1,f1) swap sends e1 to f1
    not_gamma = _isometry_from_row_matrix(_uu(), _uu(), matrix(ZZ, [[0, 1, 0, 0], [1, 0, 0, 0],
                                                                   [0, 0, 1, 0], [0, 0, 0, 1]]))
    assert not_gamma(_e1()) == _f1()
    with pytest.raises(ArithmeticError):
        _RECURSION.verify_pointwise_stabilizer_transport(marked, _v(0), not_gamma(_v(0)), not_gamma)


def test_equal_invariant_is_not_read_as_equivalence() -> None:
    r"""Sharing the invariant produces silence, never an equivalence claim.

    `v_2 = e2 + 2 f1` and `2 f1 + f2` have the same `(q, content, b(e1, .)) = (0, 1, 2)`, so
    `pointwise_stabilizer_separates` returns `False` -- "not separated", not "equivalent". (Those two happen to be
    related by the `e2 <-> f2` swap, which is consistent: silence is compatible with either answer.) The point is
    the API surface: there is no `pointwise_stabilizer_equivalent`, so the silence cannot be misread as a decision.
    """
    marked = (_e1(),)
    source = _v(2)
    same_invariant = _uu().scalar_multiple(2, _f1()) + _f2()   # coordinates (0, 2, 0, 1)

    assert _RECURSION.pointwise_stabilizer_invariant(marked, source) == (0, 1, 2)
    assert _RECURSION.pointwise_stabilizer_invariant(marked, same_invariant) == (0, 1, 2)
    assert not _RECURSION.pointwise_stabilizer_separates(marked, source, same_invariant)

    assert not hasattr(_RECURSION, "pointwise_stabilizer_equivalent"), (
        "no equivalence decision is provided for the pointwise stabilizer: adding one that reads a non-separation "
        "as equivalence would report unproved orbits"
    )
