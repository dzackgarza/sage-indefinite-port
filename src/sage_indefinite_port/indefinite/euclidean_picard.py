r"""Generators of the Euclidean Picard groups ``PU(2,1; O_d)``, ``d = 2, 7, 11`` (Zhao, arXiv:1006.3331).

Theorems 5.2-5.4 of the source give five matrices generating the projective
group. The kernel of ``U(2,1; O_d) -> PU(2,1; O_d)`` is the scalar units
``+-1``, so the five lifts and ``-I`` generate ``U(2,1; O_d)``. The form is
``z_1 conj(w_3) + z_2 conj(w_2) + z_3 conj(w_1)``.
"""

from sage.matrix.constructor import matrix
from sage.rings.number_field.number_field import QuadraticField
from sage_extensions.categories.hermitian_lattices import HermitianLattice


def field(d: int):
    r"""Return ``Q(sqrt(-d))`` for a source value of ``d``."""
    if d not in (2, 7, 11):
        raise ValueError(f"the source covers d = 2, 7, 11, not {d}")
    return QuadraticField(-d, "a")


def form(d: int):
    r"""Return the source Hermitian Gram matrix."""
    E = field(d)
    return matrix(E, 3, 3, [0, 0, 1, 0, 1, 0, 1, 0, 0])


def order_basis(d: int) -> tuple:
    r"""Return a ``ZZ``-basis of the maximal order."""
    E = field(d)
    return (E.one(), E.gen() if d == 2 else (1 + E.gen()) / 2)


def generators(d: int) -> dict:
    r"""Return the five source matrices of Theorems 5.2-5.4."""
    E = field(d)
    a = E.gen()
    omega = (1 + a) / 2
    bar = omega.conjugate()

    def entries(*rows):
        return matrix(E, 3, 3, [entry for row in rows for entry in row])

    I0 = entries((0, 0, 1), (0, -1, 0), (1, 0, 0))
    R1 = entries((1, 0, 0), (0, -1, 0), (0, 0, 1))
    T = entries((1, 0, a), (0, 1, 0), (0, 0, 1))
    match d:
        case 2:
            R2 = entries((1, 2, -2), (0, -1, 2), (0, 0, 1))
            R3 = entries((1, -a, -1), (0, -1, a), (0, 0, 1))
        case 7:
            R2 = entries((1, 1, -bar), (0, -1, 1), (0, 0, 1))
            R3 = entries((1, bar, -1), (0, -1, omega), (0, 0, 1))
        case 11:
            R2 = entries((1, 1, -bar), (0, -1, 1), (0, 0, 1))
            R3 = entries((1, bar, -1 - bar), (0, -1, omega), (0, 0, 1))
    return {"I0": I0, "R1": R1, "R2": R2, "R3": R3, "T": T}


def is_element(d: int, g) -> bool:
    r"""Return whether ``g`` is an integral unitary matrix of the source form."""
    order = field(d).ring_of_integers()
    return all(entry in order for entry in g.list()) and g.conjugate_transpose() * form(d) * g == form(d)


def full_unitary_generators(d: int) -> tuple:
    r"""Return the five lifts with the scalar kernel generator ``-I``."""
    E = field(d)
    return (*generators(d).values(), -matrix.identity(E, 3))


def lattice(d: int) -> HermitianLattice:
    r"""Return ``O_d^3`` with the source form."""
    return HermitianLattice(order_basis(d), form(d).transpose())
