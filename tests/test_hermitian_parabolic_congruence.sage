r"""Level structure on the integral parabolic of the Gaussian line in Z[i]^3, J = antidiag(1, 1, 1).

The Schreier generators of Gamma(2) in the setwise stabilizer reduce to the identity and preserve
the lattice and the form, and the index is the order of the finite image. The pointwise
stabilizer fixes the marked vector, and the setwise stabilizer contains a unit action on it.
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.matrix.special import identity_matrix as sage_identity
from sage.rings.number_field.number_field import CyclotomicField as SageCyclotomicField
from sage_extensions.categories.hermitian_lattices import HermitianLattice, HermitianWittFrame
from sage_indefinite_port.indefinite.picard_congruence import congruence_generators, finite_image, reduction


def _frame():
    field = SageCyclotomicField(int(4))
    zero, one = field.zero(), field.one()
    gram = sage_matrix(field, int(3), int(3), [zero, zero, one, zero, one, zero, one, zero, zero])
    lattice = HermitianLattice((one, field.gen()), gram)
    return field, HermitianWittFrame(lattice, sage_identity(field, int(3)), sage_matrix(field, int(1), int(1), [one]), int(1))


def test_level_two_subgroup_of_the_setwise_stabilizer():
    field, frame = _frame()
    generators = tuple(frame.stabilizer_generators(False))
    one = sage_identity(field, int(3))
    reduce = reduction(int(2))
    image = finite_image(generators, reduce, one)
    index, schreier = congruence_generators(generators, reduce, one, [one])
    assert int(index) == len(image) and len(image) > int(1)
    assert all(reduce(g) == reduce(one) and frame.preserves_lattice(g) for g in schreier)
    assert all(g.conjugate_transpose() * frame.j() * g == frame.j() for g in schreier)


def test_pointwise_stabilizer_fixes_the_marked_vector():
    field, frame = _frame()
    e = sage_matrix(field, int(3), int(1), [field.one(), field.zero(), field.zero()])
    assert all(g * e == e for g in frame.stabilizer_generators(True))
    units = [g[int(0), int(0)] for g in frame.stabilizer_generators(False)]
    assert field.gen() in units or -field.gen() in units
