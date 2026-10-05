from dzack_research.preamble.all import Lattices, ZZ

import pytest

from sage_indefinite_port.indefinite.eichler import (
    EichlerOrbitCover,
    InfiniteLocusError,
    eichler_transvection,
    find_hyperbolic_pair,
    square_divisors,
)


def test_eichler_transvection_delegates_to_preamble_isometry() -> None:
    lattice = Lattices(ZZ)("U") + Lattices(ZZ)("U") + Lattices(ZZ)("E8")
    e, _f, _e_prime, _f_prime, *complement = lattice.module_generators()
    x = complement[0]
    transvection = eichler_transvection(e, x)

    assert transvection.parent() is lattice.Aut()
    assert transvection(e) == e
    assert lattice.b(e, x) == ZZ.zero()


def test_square_divisors_match_nonprimitive_norm_decomposition() -> None:
    assert square_divisors(72) == (1, 2, 3, 6)
    assert square_divisors(-72) == (1, 2, 3, 6)


def test_two_u_cover_preserves_primitive_and_nonprimitive_semantics() -> None:
    complement = Lattices(ZZ)("A2")
    model = EichlerOrbitCover(complement.two_u_eichler_model())

    primitive = model.covering_representatives(2, primitive=True)
    nonprimitive = model.covering_representatives(8, primitive=False)
    primitive_norm = primitive.representatives[0].q()
    nonprimitive_norm = nonprimitive.representatives[0].q()

    assert all(vector.is_primitive() and vector.q() == primitive_norm for vector in primitive)
    assert all(vector.q() == nonprimitive_norm for vector in nonprimitive)
    with pytest.raises(InfiniteLocusError):
        model.covering_representatives(0, primitive=False)


def test_find_hyperbolic_pair_matches_isotropic_fixture_cases() -> None:
    isotropic = Lattices(ZZ)([[4, 0, 0], [0, 0, -1], [0, -1, -2]])
    v, w = find_hyperbolic_pair(isotropic)
    assert v.is_primitive() and v.is_isotropic()
    assert w.is_isotropic()
    assert isotropic.b(v, w) > ZZ.zero()

    anisotropic = Lattices(ZZ)([[-516, 36, 72], [36, -2, -5], [72, -5, -10]])
    with pytest.raises(ValueError):
        find_hyperbolic_pair(anisotropic)


def test_choose_splitting_vector_returns_positive_vector() -> None:
    model = EichlerOrbitCover(Lattices(ZZ)("A2").two_u_eichler_model())
    vector = model.choose_splitting_vector()
    assert vector.q() > model.lattice().base_ring().zero()
