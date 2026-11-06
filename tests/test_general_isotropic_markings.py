"""Maximal isotropic sublattices of 3U admit ambient transporters."""
from sage_indefinite_port.indefinite.marked_lattices import ARROWS, isometry
from tests.acceptance.consumer import lattice


def test_three_u_maximal_isotropic_sublattices_have_a_transporter():
    gram = [[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)]
    ambient = lattice(gram)
    basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    source = ambient.primitive_isotropic_subobject(*basis[::2])
    target = ambient.primitive_isotropic_subobject(*basis[1::2])
    witness = isometry(ARROWS, ARROWS.object(source.inclusion()), ARROWS.object(target.inclusion()))
    assert witness is not None
    for point in basis[::2]:
        assert target.inclusion().is_in_image(witness(point))
    for point in basis[1::2]:
        assert source.inclusion().is_in_image((~witness)(point))
    for left in basis:
        for right in basis:
            assert ambient.b(witness(left), witness(right)) == ambient.b(left, right)


def test_three_u_partial_flags_have_a_transporter():
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS

    gram = [[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)]
    ambient = lattice(gram)
    basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    source = ambient.primitive_isotropic_subobject(*basis[::2])
    target = ambient.primitive_isotropic_subobject(*basis[1::2])
    source_line = basis[0] + basis[2]
    target_line = basis[1] + basis[3]
    source_term = source.inclusion().domain().primitive_isotropic_subobject(source.inclusion().lift(source_line))
    target_term = target.inclusion().domain().primitive_isotropic_subobject(target.inclusion().lift(target_line))
    source_marking = FLAGS((source_term.inclusion(), source.inclusion()))
    target_marking = FLAGS((target_term.inclusion(), target.inclusion()))
    witness = isometry(FLAGS, source_marking, target_marking)
    assert witness is not None
    target_ambient_line = ambient.primitive_isotropic_subobject(target_line)
    assert target_ambient_line.inclusion().is_in_image(witness(source_line))
    for point in basis[::2]:
        assert target.inclusion().is_in_image(witness(point))
    for left in basis:
        for right in basis:
            assert ambient.b(witness(left), witness(right)) == ambient.b(left, right)


def test_three_u_complete_flags_have_a_transporter():
    from sage_extensions.all import Cat, Lattices, ZZ
    from sage_extensions.categories.abstract_categories.products import FiniteOrdinalCategory

    gram = [[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)]
    ambient = lattice(gram)
    basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())
    category = Cat().Mor(FiniteOrdinalCategory(4), Lattices(ZZ))

    def marking(vectors):
        space = ambient.primitive_isotropic_subobject(*vectors)
        domain = space.inclusion().domain()
        plane = domain.primitive_isotropic_subobject(*(space.inclusion().lift(point) for point in vectors[:2]))
        line = plane.inclusion().domain().primitive_isotropic_subobject(plane.inclusion().lift(space.inclusion().lift(vectors[0] + vectors[1])))
        return category((line.inclusion(), plane.inclusion(), space.inclusion()))

    witness = isometry(category, marking(basis[::2]), marking(basis[1::2]))
    assert witness is not None
    for size in (1, 2, 3):
        source_vectors = (basis[0] + basis[2],) if size == 1 else basis[::2][:size]
        target_vectors = (basis[1] + basis[3],) if size == 1 else basis[1::2][:size]
        target_term = ambient.primitive_isotropic_subobject(*target_vectors)
        assert all(target_term.inclusion().is_in_image(witness(point)) for point in source_vectors)
    assert all(ambient.b(witness(left), witness(right)) == ambient.b(left, right) for left in basis for right in basis)


def test_isotropic_rank_above_witt_index_has_no_markings():
    from sage_indefinite_port.indefinite.marked_lattices import primitive_marking_classes

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    shape = lattice([[0 for _ in range(4)] for _ in range(4)])
    assert tuple(primitive_marking_classes(ARROWS, shape, ambient)) == ()


def test_three_u_general_induction_classifies_maximal_isotropic_spaces():
    from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    representatives = IsotropicInduction(EichlerRecursion()).isotropic_flag_bases(ambient, (3,))
    assert len(representatives) == 1
    for basis in representatives:
        assert len(basis) == 3
        assert all(ambient.b(left, right) == 0 for left in basis for right in basis)
        assert ambient.primitive_isotropic_subobject(*basis).inclusion().domain().rank() == 3


def test_complete_flag_marking_retains_each_primitive_term():
    from sage_indefinite_port.indefinite.marked_lattices import _flag_marking_from_basis

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    basis = tuple(ambient.module_generator(label) for label in ambient.module_generating_set())[::2]
    arrows = _flag_marking_from_basis(ambient, basis, (1, 2, 3))
    assert tuple(int(arrow.domain().rank()) for arrow in arrows) == (1, 2, 3)
    assert arrows[-1].codomain() is ambient
    assert all(arrow.is_primitive() for arrow in arrows)
    assert all(left.codomain() is right.domain() for left, right in zip(arrows, arrows[1:]))


def test_partial_flag_above_witt_index_has_no_markings():
    from sage_indefinite_port.indefinite.marked_lattices import FLAGS, primitive_marking_classes

    ambient = lattice([[0, 1], [1, 0]])
    plane = lattice([[0, 0], [0, 0]])
    line = plane.primitive_isotropic_subobject(next(iter(plane.module_generators())))
    assert tuple(primitive_marking_classes(FLAGS, line.inclusion(), ambient)) == ()


def test_three_u_general_induction_classifies_complete_and_partial_flags():
    from sage_indefinite_port.indefinite.isotropic_induction import IsotropicInduction
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion

    ambient = lattice([[int(i // 2 == j // 2 and i != j) for j in range(6)] for i in range(6)])
    induction = IsotropicInduction(EichlerRecursion())
    for dimensions in ((1, 3), (1, 2, 3)):
        representatives = induction.isotropic_flag_bases(ambient, dimensions)
        assert len(representatives) == 1
        for basis in representatives:
            assert len(basis) == 3
            assert all(ambient.b(left, right) == 0 for left in basis for right in basis)
            for rank in dimensions:
                assert ambient.primitive_isotropic_subobject(*basis[:rank]).inclusion().domain().rank() == rank
