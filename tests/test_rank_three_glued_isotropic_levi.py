"""Source-faithful isotropic induction with a proper integral rank-three Levi image.

Let L = U^2 + U(3), with I = <e_1,e_2,e_3>.  On the ordered
coordinates of I, the kernel of its mod-three pairing with L is the
line <e_3>.  The integral Levi image Stab_O(L)(I) -> GL(I) is exactly the
preimage of the stabilizer of that line in GL_3(F_3): preservation is
necessary, and every such A extends by the integral contragredient
D A^(-T) D^(-1) on <f_1,f_2,f_3>, D=diag(1,1,3).
Consequently its index is |P^2(F_3)| = (3^3-1)/(3-1) = 13.

The projective line of the pulled-back pairing kernel is a complete,
independently computed invariant for the based frames of this I.
"""

from sage.all import GF, matrix
from sage_extensions.all import Cat, Lattices, ZZ
from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory
from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
from sage_indefinite_port.indefinite.marked_lattices import (
    _flag_action,
    primitive_marking_classes,
)


def test_glued_rank_three_based_isotropic_classes_detect_proper_levi_image():
    lattices = Lattices(ZZ)
    ambient = lattices("U") + lattices("U") + lattices([[0, 3], [3, 0]])
    source = lattices([[0] * 3 for _ in range(3)])
    category = CosliceCategory(lattices, source)
    classes = tuple(primitive_marking_classes(category, source, ambient))
    field = GF(3)
    projective_lines = set()
    for marked in classes:
        arrow = marked.arrow()
        assert arrow.is_primitive() and arrow.codomain() is ambient
        basis = tuple(arrow(point) for point in source.module_generators())
        pairing = matrix(field, [
            [int(ambient.b(v, point)) % 3 for point in ambient.module_generators()]
            for v in basis
        ])
        kernel = pairing.left_kernel()
        assert kernel.dimension() == 1
        coordinates = tuple(kernel.basis()[0])
        pivot = next(entry for entry in coordinates if entry)
        projective_lines.add(tuple(int(entry / pivot) for entry in coordinates))
    assert len(classes) == len(projective_lines) == 13


def test_glued_rank_three_partial_flags_split_by_first_line_divisibility():
    r"""For I in U^2+U(3), the Levi fixes <e_3> mod 3 and acts
    transitively on the remaining projective points. Its two orbits on
    primitive first lines have divisibility 3 and 1 respectively.
    """
    lattices = Lattices(ZZ)
    ambient = lattices("U") + lattices("U") + lattices([[0, 3], [3, 0]])
    source = lattices([[0] * 3 for _ in range(3)])
    first = source.primitive_isotropic_subobject(next(iter(source.module_generators())))
    category = Cat().Mor(FiniteOrdinalCategory(3), lattices)
    classes = tuple(primitive_marking_classes(category, first.inclusion(), ambient))
    divisibilities = set()
    for marked in classes:
        dimensions, basis, _images = _flag_action(marked)
        assert dimensions == (1, 3)
        divisibilities.add(int(basis[0].divisor()))
    assert len(classes) == 2 and divisibilities == {1, 3}
