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


def test_binary_file_vertex_and_file_vertex_roots() -> None:
    L = hyperbolic([[0, 1], [1, 0]])
    k, r = element(L, (1, -1)), element(L, (1, 1))
    assert tuple(int(value) for value in (k.q(), r.q(), L.b(k, r))) == (-2, 2, 0)

    # The displayed root is reflective on the integral lattice: reflection is
    # -1 on r=e0+e1 and +1 on k=e0-e1, so it maps (e0,e1) to (-e1,-e0).
    assert r - L.b(r, r) * r == -r
    assert k - L.b(k, r) * r == k

    default = L.allcock_edgewalk()
    assert not default.is_reflective()
    assert len(default.vertices()) == 0
    assert default.reason_non_reflective() == "list of roots is empty"

    # FileVertex and FileVertexRoots must not be confused with the default
    # isotropic_vinberg start. O(U) is finite (its null-line scalings are
    # integral units), so the exact nonreflective verdict has infinite
    # covolume, not a fictitious empty wall list at the prescribed finite k.
    for prescribed in (k, 2 * k, -3 * k, (k, (r,)), (2 * k, (r,)), (-3 * k, (r,))):
        report = L.allcock_edgewalk(initial_vertex=prescribed)
        assert not report.is_reflective()
        assert len(report.simple_roots()) == 0
        assert "finite integral orthogonal group" in report.reason_non_reflective()
        assert len(report.vertices()) == 2
        rays = [ray for vertex in report.vertices() for ray in vertex.primitive_rays()]
        assert all(int(ray.content()) == 1 for ray in rays)
        assert any(int(ray.q()) == -2 for ray in rays)

    # A bare null FileVertex follows get_simple_cone, including ray
    # normalisation, and returns the source's cusp nonreflectivity report.
    for null in (element(L, (0, 1)), element(L, (0, 2)), element(L, (0, -3))):
        report = L.allcock_edgewalk(initial_vertex=null)
        assert not report.is_reflective()
        assert len(report.vertices()) == 0
        assert report.reason_non_reflective() == "list of roots is empty"

    # The empty-wall isotropic fallback used above must not accept bad roots
    # merely because the ambient binary lattice is isotropic.
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
    # The transporter is explicit, so this checks the ENTIRE integral root
    # multiset and its lattice gluing rather than merely a Coxeter graph.
    # It is a signed swap of the last two coordinates, with determinant -1.
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
    # For q=2*x^2-y^2, P sends (e0,e1) to (3,4),(2,3). Direct evaluation
    # gives q(3,4)=2, q(2,3)=-1 and b((3,4),(2,3))=0; det(P)=1.
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
        anisotropic_cycle((1, 0, -1), 1, (1, 0))  # split: discriminant one
    with pytest.raises(ValueError, match="signature"):
        anisotropic_cycle((1, 0, 1), 1, (1, 0))  # not Lorentzian
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
