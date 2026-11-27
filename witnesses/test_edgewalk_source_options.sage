r"""Source options and binary guards of Allcock's edgewalk (#31).

Pinned semantics: ``edgewalk.h:2100-2143`` normalises the generator for both
FileVertex and FileVertexRoots, while ``:1764-1774`` runs its group finiteness
tester only when EarlyTerminationIfNotReflective is true. The rank-two cycle
is ``two_dim_lorentzian.h:70-326``.

These are implementer diagnostics, not replacements for the independent
unchanged cases in tests/acceptance/test_reflective.sage. The exhaustive-mode
check does compare its whole integral chamber with Allcock's RK3_all line 13.
"""

import pytest
from sage.all import ZZ, matrix, vector
from sage_extensions.all import ZZ as PreambleZZ, HyperbolicLattices

from tests.acceptance.consumer import coordinates, element, lattice
from tests.fixtures.oracle_fixtures import load_allcock_rank3_reflective
from sage_indefinite_port.indefinite.edgewalk import _Edgewalk
from sage_indefinite_port.indefinite.edgewalk_rank2 import _FormBinaryCycle, anisotropic_cycle


def hyperbolic(rows):
    return HyperbolicLattices(PreambleZZ)(lattice(rows))


def test_unreduced_vinberg_start_retains_the_original_corner() -> None:
    r"""ApplyReduction=false is the identity presentation (Indefinite_LLL.h:934).

    GetV0_vector still calls GetShortIntegralVector, which reduces -G
    internally (vinberg.h:307-311, Positivity.h:410-482,516-534).
    With the port's fixed reduction seed the internal basis is
    ((0,0,-1),(1,-2,-2),(0,-1,0)); the first positive diagonal direction
    of the reduced -G is (1,0,1). Its normalized approximation rounds
    both nonzero coordinates to 1 at the first step, yielding
    (0,-1,-1) in the unreduced input. Verify the entire chamber by
    an integral ambient transporter, not just the count.
    """
    from itertools import permutations

    case = next(case for case in load_allcock_rank3_reflective()
                if case["source"]["line"] == 13)
    gram = matrix(ZZ, case["gram"])
    basis_change = matrix(ZZ, [[1, 2, 2], [0, 1, 0], [0, 0, 1]])
    assert basis_change.det() == 1
    changed = basis_change * gram * basis_change.transpose()
    inner = matrix(ZZ, [[0, 0, -1], [1, -2, -2], [0, -1, 0]])
    assert abs(inner.det()) == 1
    assert inner * (-changed) * inner.transpose() == matrix(ZZ, [[0, 0, 1], [0, -2, 0], [1, 0, 0]])
    diagonalizer = matrix(ZZ, [[0, 1, 0], [1, 0, 1], [-1, 0, 1]])
    assert diagonalizer * inner * (-changed) * inner.transpose() * diagonalizer.transpose() == matrix(ZZ, [[-2, 0, 0], [0, 2, 0], [0, 0, -2]])
    assert tuple((diagonalizer * inner).row(1)) == (0, -1, -1)
    recorded = matrix(ZZ, case["simple_roots"]) * basis_change.inverse()
    for sign in (1, -1):
        L = hyperbolic((sign * changed).rows())
        report = L.allcock_edgewalk(option_initial_vertex="vinberg", apply_reduction=False)
        assert report.is_reflective()
        assert report.has_computed_simple_roots()
        first = tuple(report.vertices()[0].primitive_rays())[0]
        assert coordinates(first) == (0, -1, -1)
        assert int(first.q()) == -2 * sign
        actual = matrix(ZZ, [coordinates(root) for root in report.simple_roots()])
        assert actual.nrows() == recorded.nrows() == 3
        assert any(
            all(x in ZZ for x in transport.list())
            and abs(transport.det()) == 1
            and transport * changed * transport.transpose() == changed
            for transport in (
                recorded.inverse() * actual.matrix_from_rows(order)
                for order in permutations(range(3))
            )
        )
        reduced_start = L.allcock_edgewalk(option_initial_vertex="vinberg")
        assert reduced_start.is_reflective()
        start = tuple(reduced_start.vertices()[0].primitive_rays())[0]
        resumed = L.allcock_edgewalk(initial_vertex=start)
        assert tuple(coordinates(r) for r in reduced_start.simple_roots()) == tuple(
            coordinates(r) for r in resumed.simple_roots()
        )
        assert tuple(
            (tuple(coordinates(r) for r in v.primitive_rays()),
             tuple(coordinates(r) for r in v.wall_roots()))
            for v in reduced_start.vertices()
        ) == tuple(
            (tuple(coordinates(r) for r in v.primitive_rays()),
             tuple(coordinates(r) for r in v.wall_roots()))
            for v in resumed.vertices()
        )


def test_anisotropic_binary_selector_uses_the_same_vinberg_branch() -> None:
    r"""For q=2*x^2-y^2, the chamber in y>0 lies between x=0 and x=-y/2.

    The walls (1,0), (-1,1) have squares 2 and 1 and are reflective.
    An integral primitive reflective root (a,b) has positive square
    dividing 4 (twice the exponent of the discriminant group). Square
    4 forces both coordinates even, so the possible squares are 1,2.
    A wall strictly between these endpoints would have 0<|b|<|a|,
    hence square >a^2, impossible for |a|>=2, and no integral b exists
    for |a|=1. Thus these are adjacent walls of the chamber.

    The isotropic_vinberg selector must fall back to precisely the same
    Vinberg branch when the form is anisotropic (edgewalk.h:2119-2126).
    The shear makes bypassing pre-reduction observable.
    """
    from itertools import permutations

    gram = matrix(ZZ, [[2, 0], [0, -1]])
    change = matrix(ZZ, [[2, 3], [1, 2]])
    assert change.det() == 1
    changed = change * gram * change.transpose()
    expected = matrix(ZZ, [[1, 0], [-1, 1]]) * change.inverse()
    L = hyperbolic(changed.rows())
    assert not L.is_isotropic()
    for reduce in (True, False):
        explicit = L.allcock_edgewalk(option_initial_vertex="vinberg", apply_reduction=reduce)
        fallback = L.allcock_edgewalk(option_initial_vertex="isotropic_vinberg", apply_reduction=reduce)
        for report in (explicit, fallback):
            assert report.is_reflective()
            actual = matrix(ZZ, [coordinates(root) for root in report.simple_roots()])
            assert actual.nrows() == 2
            assert any(
                all(x in ZZ for x in transport.list())
                and abs(transport.det()) == 1
                and transport * changed * transport.transpose() == changed
                for transport in (
                    expected.inverse() * actual.matrix_from_rows(order)
                    for order in permutations(range(2))
                )
            )
        assert tuple(coordinates(r) for r in explicit.simple_roots()) == tuple(
            coordinates(r) for r in fallback.simple_roots()
        )
        assert tuple(tuple(coordinates(r) for r in v.primitive_rays()) for v in explicit.vertices()) == tuple(
            tuple(coordinates(r) for r in v.primitive_rays()) for v in fallback.vertices()
        )
        if not reduce:
            assert coordinates(tuple(explicit.vertices()[0].primitive_rays())[0]) == (-1, 0)


def test_binary_file_vertex_and_file_vertex_roots() -> None:
    L = hyperbolic([[0, 1], [1, 0]])
    k, r = element(L, (1, -1)), element(L, (1, 1))
    assert tuple(int(value) for value in (k.q(), r.q(), L.b(k, r))) == (-2, 2, 0)

    assert r - L.b(r, r) * r == -r
    assert k - L.b(k, r) * r == k

    default = L.allcock_edgewalk()
    assert not default.is_reflective()
    assert len(default.vertices()) == 0
    assert default.reason_non_reflective() == "list of roots is empty"

    for prescribed in (k, 2 * k, -3 * k, (k, (r,)), (2 * k, (r,)), (-3 * k, (r,))):
        report = L.allcock_edgewalk(initial_vertex=prescribed)
        assert not report.is_reflective()
        assert len(report.simple_roots()) == 0
        assert "finite integral orthogonal group" in report.reason_non_reflective()
        assert len(report.vertices()) == 2
        rays = [ray for vertex in report.vertices() for ray in vertex.primitive_rays()]
        assert all(int(ray.content()) == 1 for ray in rays)
        assert any(int(ray.q()) == -2 for ray in rays)

    for null in (element(L, (0, 1)), element(L, (0, 2)), element(L, (0, -3))):
        report = L.allcock_edgewalk(initial_vertex=null)
        assert not report.is_reflective()
        assert len(report.vertices()) == 0
        assert report.reason_non_reflective() == "list of roots is empty"

    with pytest.raises(ValueError, match="invalid wall roots"):
        L.allcock_edgewalk(initial_vertex=(k, (element(L, (1, 0)),)))


def test_exhaustive_mode_skips_early_finiteness_and_preserves_chamber(monkeypatch) -> None:
    case = next(
        case for case in load_allcock_rank3_reflective()
        if case["source"]["line"] == 13
    )
    assert case["gram"] == [[2, 0, 0], [0, 0, -1], [0, -1, 0]]
    L = hyperbolic(case["gram"])

    def forbidden(*args):
        raise AssertionError("EarlyTerminationIfNotReflective=false invoked the early-only finiteness test")

    monkeypatch.setattr(_Edgewalk, "common_fixed_space", forbidden)
    report = L.allcock_edgewalk(early_termination=False)
    assert report.is_reflective()
    assert report.reason_non_reflective() is None
    computed = [coordinates(root) for root in report.simple_roots()]
    assert len(computed) == 3
    gram = matrix(ZZ, case["gram"])
    recorded = matrix(ZZ, case["simple_roots"])
    transport = matrix(ZZ, [[1, 0, 0], [0, 0, -1], [0, -1, 0]])
    assert transport * gram * transport.transpose() == gram
    assert transport.det() == -1
    assert sorted(tuple(row) for row in (recorded * transport).rows()) == sorted(tuple(row) for row in matrix(ZZ, computed).rows())


def test_default_mode_retains_the_exact_nonreflectivity_certificate() -> None:
    report = hyperbolic([[0, 1, 0], [1, 0, 0], [0, 0, 46]]).allcock_edgewalk()
    assert not report.is_reflective()
    assert len(report.vertices()) == 3
    assert len(report.simple_roots()) == 0
    assert "no nonzero common fixed timelike vector" in report.reason_non_reflective()


def test_anisotropic_cycle_period_and_input_guards() -> None:
    result = anisotropic_cycle((2, 0, -1), 2, (1, 0))
    assert result is not None
    period, representatives = result
    assert period == matrix(ZZ, [[3, 4], [2, 3]])
    assert representatives == ((1, 1), (3, 4))
    gram = matrix(ZZ, [[2, 0], [0, -1]])
    assert period * gram * period.transpose() == gram
    assert period.det() == 1
    assert [int(vector(ZZ, row) * gram * vector(ZZ, row)) for row in representatives] == [1, 2]

    with pytest.raises(ValueError, match="anisotropic form"):
        anisotropic_cycle((1, 0, -1), 1, (1, 0))
    with pytest.raises(ValueError, match="signature"):
        anisotropic_cycle((1, 0, 1), 1, (1, 0))
    with pytest.raises(ValueError, match="primitive positive"):
        anisotropic_cycle((2, 0, -1), 2, (2, 0))
    with pytest.raises(ValueError, match="primitive positive"):
        anisotropic_cycle((2, 0, -1), 2, (0, 1))
    with pytest.raises(ValueError, match="positive integer"):
        anisotropic_cycle((2, 0, -1), 2, (1, 0), scale=0)

    cycle = _FormBinaryCycle((2, 0, -1))
    with pytest.raises(ValueError, match="oriented unimodular pair"):
        cycle.promised_step(2, (1, 0), (0, 2))
    with pytest.raises(ArithmeticError, match="positive right norm"):
        cycle.canonical_companion(2, (0, 1), (1, 0))
