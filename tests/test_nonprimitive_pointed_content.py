"""Exact content reduction for nonprimitive pointed lattice isometries (#48)."""
from tests.acceptance.consumer import lattice
from sage_indefinite_port.indefinite.recursive import EichlerRecursion


def test_imprimitive_isotropic_transport_cancels_content():
    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    recursion = EichlerRecursion()
    witness = recursion.pointed_isometry(2 * e, 2 * f)
    assert witness is not None
    assert witness(2 * e) == 2 * f
    assert recursion.pointed_isometry(2 * e, 4 * f) is None


def test_imprimitive_nonzero_norm_transport_cancels_content():
    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    recursion = EichlerRecursion()
    witness = recursion.pointed_isometry(2 * (e + 2 * f), 2 * (2 * e + f))
    assert witness is not None
    assert witness(2 * (e + 2 * f)) == 2 * (2 * e + f)
    assert recursion.pointed_isometry(2 * (e + f), e + 4 * f) is None


def test_imprimitive_isotropic_stabilizer_fixes_original_mark():
    plane = lattice([[0, 1], [1, 0]])
    e, f = tuple(plane.module_generators())
    family = EichlerRecursion().pointed_automorphism_generators(3 * e)
    assert all(g(3 * e) == 3 * e and g(e) == e for g in family)
    assert all(plane.b(g(x), g(y)) == plane.b(x, y) for g in family for x in (e, f) for y in (e, f))
