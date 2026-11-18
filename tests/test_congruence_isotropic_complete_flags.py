r"""Complete isotropic flags in 2U modulo the level-two principal subgroup.

The quadratic space L/2L has O^+(4,2) of order 72. Its maximal singular
planes form two rulings of three planes each; each plane has three lines.
The integral group O(2U) is transitive on primitive complete null flags
by hyperbolic completion. Its mod-two image is all O^+(4,2), certified
by the public representation order, so its level-two kernel has 18
complete-flag orbits, with one class per finite flag.
"""

from sage.all import GF, matrix
from sage_extensions.all import Cat, Lattices, ZZ
from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory
from sage_indefinite_port.indefinite.marked_lattices import (
    _flag_action,
    primitive_marking_classes,
)


def test_two_u_complete_null_flag_classes_under_gamma_two():
    lattices = Lattices(ZZ)
    ambient = lattices("U") + lattices("U")
    group = ambient.O()
    assert int(group.congruence_representation(2).codomain().cardinality()) == 72
    gamma = group.congruence_preimage(2)
    source = lattices([[0, 0], [0, 0]])
    first = source.primitive_isotropic_subobject(next(iter(source.module_generators())))
    category = Cat().Mor(FiniteOrdinalCategory(3), lattices)
    full = tuple(primitive_marking_classes(category, first.inclusion(), ambient))
    assert len(full) == 1
    classes = tuple(primitive_marking_classes(category, first.inclusion(), ambient, acting_group=gamma))
    assert len(classes) == 18
    labels = tuple(ambient.module_generating_set())
    residues = set()
    for marked in classes:
        ranks, basis, _actions = _flag_action(marked)
        assert ranks == (1, 2)
        assert all(ambient.b(x, y) == 0 for x in basis for y in basis)
        residue = tuple(
            tuple(tuple(int(c) for c in row) for row in matrix(
                GF(2), [[int(point.to_vector()(label)) % 2 for label in labels] for point in basis[:rank]]
            ).echelon_form().rows())
            for rank in ranks
        )
        residues.add(residue)
    assert len(residues) == 18


def test_two_u_based_maximal_null_planes_under_gamma_two():
    r"""Six finite null planes each have six ordered mod-two bases, so the
    Γ(2)-orbits of primitive *based* maximal null planes number 36.
    This is twice the number of complete-flag classes, as required by
    the full source GL2 action versus its parabolic flag subgroup.
    """
    from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory

    lattices = Lattices(ZZ)
    ambient = lattices("U") + lattices("U")
    source = lattices([[0, 0], [0, 0]])
    category = CosliceCategory(lattices, source)
    gamma = ambient.O().congruence_preimage(2)
    ordinary = tuple(primitive_marking_classes(category, source, ambient))
    assert len(ordinary) == 1
    classes = tuple(primitive_marking_classes(category, source, ambient, acting_group=gamma))
    assert len(classes) == 36
    labels = tuple(ambient.module_generating_set())
    rows = {
        tuple(tuple(int(arrow(point).to_vector()(label)) % 2 for label in labels)
              for point in source.module_generators())
        for marked in classes
        for arrow in (marked.arrow(),)
    }
    assert len(rows) == 36
    assert all(marked.arrow().codomain() is ambient and marked.arrow().is_primitive()
               for marked in classes)


def test_mixed_degeneration_rank_two_moving_transporter_skips_nondegenerate_invariant():
    r"""For L=U⊕0, the full isotropic plane <e,r> is carried to
    <f,r> by the integral swap e↔f. Since the characteristic radical r
    lies in each plane, the pairing L→I* has rank one, not two: the
    nondegenerate source invariant is undefined, but the degenerate
    primitive-subobject transporter remains exact.
    """
    from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry

    ambient = Lattices(ZZ)([[0, 1, 0], [1, 0, 0], [0, 0, 0]])
    e, f, radical = tuple(ambient.module_generators())
    source = ARROWS(ambient.primitive_isotropic_subobject(e, radical).inclusion())
    target = ARROWS(ambient.primitive_isotropic_subobject(f, radical).inclusion())
    witness = isometry(ARROWS, source, target)
    assert witness is not None
    assert all(target.arrow().is_in_image(witness(point)) for point in (e, radical))
    assert all(source.arrow().is_in_image((~witness)(point)) for point in (f, radical))
    assert all(ambient.b(witness(x), witness(y)) == ambient.b(x, y)
               for x in (e, f, radical) for y in (e, f, radical))
