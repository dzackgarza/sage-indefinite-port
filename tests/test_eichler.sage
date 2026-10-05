from dzack_research.preamble.all import Lattices, ZZ

from sage_indefinite_port.indefinite.eichler import eichler_transvection, square_divisors


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
