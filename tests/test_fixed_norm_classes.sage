r"""O(L)-classes of fixed-norm vectors of L = diag(2,2,-6,-6) = (ZZ<1,i,j,k>, 2 Nrd) in (-1,3).

Expected answers by hand (b(v, v) = 2 Nrd z):

- b = 2, Nrd z = 1: left multiplication by z lies in O(L) and sends 1 to z, so one class.
- b = -2, Nrd z = -1: for a fixed e of norm -1, z = (z e^{-1}) e with Nrd(z e^{-1}) = 1, so one class.
- b = 3: b(v, v) is even on L, so there is no vector.
- b = 8, Nrd z = 4: a nonprimitive z = d u has d^2 Nrd u = 4, so d = 2 and Nrd u = 1, and these form
  the single class of 2. The primitive 1+3i+j+k has Nrd 1+9-3-3 = 4, so a primitive class exists.
- (2+j)^2 = 7+4j has Nrd 49-48 = 1, so 1 and 7+4j are equivalent; 1+i has Nrd 2, so it is not.
"""

from sage_extensions.all import *
from sage.matrix.constructor import matrix as sage_matrix
from sage.rings.rational_field import QQ as SageQQ
from sage_indefinite_port.indefinite.quaternion_lattices import QuaternionNormLattice
from sage_indefinite_port.indefinite.quaternion_voronoi import QuaternionVoronoi
from sage_indefinite_port.indefinite.fixed_norm_classes import FixedNormClasses, fixed_norm_vector_classes


def _setup():
    model = QuaternionNormLattice(int(-1), int(3), int(2))
    return model, QuaternionVoronoi(int(-1), int(3)), model.normalizer_classes_in_box(int(3))


def test_norm_two_and_minus_two_each_form_one_class():
    model, voronoi, normalizers = _setup()
    assert len(fixed_norm_vector_classes(model, voronoi, normalizers, int(2), False)) == int(1)
    assert len(fixed_norm_vector_classes(model, voronoi, normalizers, int(-2), False)) == int(1)


def test_odd_norm_has_no_vectors():
    model, voronoi, normalizers = _setup()
    assert fixed_norm_vector_classes(model, voronoi, normalizers, int(3), False) == ()


def test_norm_eight_has_one_nonprimitive_class_and_a_primitive_class():
    model, voronoi, normalizers = _setup()
    classes = FixedNormClasses(model, voronoi, normalizers)
    every = classes.classes(int(4), False)
    primitive = classes.classes(int(4), True)
    assert len(every) == len(primitive) + int(1)
    witness = model.element(int(1), int(3), int(1), int(1))
    assert any(classes.transporter(p, witness) is not None for p in primitive)


def test_transporter_carries_one_to_seven_plus_four_j_and_preserves_the_form():
    model, voronoi, normalizers = _setup()
    classes = FixedNormClasses(model, voronoi, normalizers)
    one = model.element(int(1), int(0), int(0), int(0))
    target = model.element(int(7), int(0), int(4), int(0))
    columns = classes.transporter(one, target)
    gram = sage_matrix(SageQQ, [[int(model.lattice().gram_matrix()[i, j]) for j in range(int(4))] for i in range(int(4))])
    assert columns.transpose() * gram * columns == gram
    assert abs(columns.det()) == int(1) and columns.denominator() == int(1)
    assert list(columns.column(int(0))) == list(target)
    assert classes.transporter(one, model.element(int(1), int(1), int(0), int(0))) is None
