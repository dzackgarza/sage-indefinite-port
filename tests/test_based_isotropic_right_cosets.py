"""Full rational/integral right-coset action for a pointwise isotropic mark."""

from itertools import combinations

from sage_extensions.all import Lattices, ZZ
from sage_extensions.categories.abstract_categories.arrow_categories import CosliceCategory

from sage_indefinite_port.indefinite.marked_lattices import (
    ARROWS, RECURSION, isotropic_stabilizer_right_cosets,
)
from tests.acceptance.consumer import lattice


def test_based_glued_isotropic_line_has_four_distinct_integral_right_cosets():
    r"""The rational line-stabilizer index is four for this glued lattice.

    Central negation -id lies in the integral *moving* line stabilizer
    and sends e to -e. Hence the moving rational group is obtained from
    the pointwise rational group by adjoining this integral element;
    its integral subgroup is extended by the same element. Thus their
    right-coset spaces have the same four elements, although the actual
    rational stabilizer groups and based/moving markings differ.
    """
    ambient = lattice([
        [0, 2, 0, 0], [2, 0, 1, 0],
        [0, 1, 0, 2], [0, 0, 2, 0],
    ])
    point = tuple(ambient.module_generators())[0]
    source = Lattices(ZZ)([[0]])
    category = CosliceCategory(Lattices(ZZ), source)
    based = category.object(source.hom([point], codomain=ambient))
    moving = ARROWS.object(ambient.primitive_isotropic_subobject(point).inclusion())

    moving_representatives = tuple(isotropic_stabilizer_right_cosets(ARROWS, moving))
    based_representatives = tuple(isotropic_stabilizer_right_cosets(category, based))
    assert len(moving_representatives) == len(based_representatives) == 4

    based_action = RECURSION.isotropic_stabilizer_action((point,), ())
    inclusion = based_action.lattice_inclusion()
    assert all(representative.domain() is representative.codomain()
               for representative in based_representatives)
    assert all(
        not based_action._carries_lattice((~first) * second, inclusion, inclusion)
        for first, second in combinations(based_representatives, 2)
    )


def test_based_unimodular_maximal_plane_has_one_integral_right_coset():
    r"""For 2U, fixing the two primitive null-plane basis vectors
    leaves exactly the skew-unipotent pointwise parabolic. Its integral
    generator has one skew parameter, so the rational right-coset
    construction introduces no additional integral lattice presentations.
    """
    ambient = Lattices(ZZ)("U") + Lattices(ZZ)("U")
    e, _f, u, _v = tuple(ambient.module_generators())
    source = Lattices(ZZ)([[0, 0], [0, 0]])
    category = CosliceCategory(Lattices(ZZ), source)
    based = category.object(source.hom((e, u), codomain=ambient))
    representatives = tuple(isotropic_stabilizer_right_cosets(category, based))
    assert len(representatives) == 1
    rational = representatives[0].domain()
    assert representatives[0].codomain() is rational
    assert all(rational.b(representatives[0](x), representatives[0](y)) == rational.b(x, y)
               for x in rational.module_generators() for y in rational.module_generators())


def test_based_rank_zero_isotropic_mark_has_unique_right_coset():
    r"""The zero sublattice has perpendicular L itself, so its rational
    extension group already preserves L and G/H is a singleton.
    """
    ambient = Lattices(ZZ)("U")
    source = Lattices(ZZ)(0)
    category = CosliceCategory(Lattices(ZZ), source)
    marked = category.object(source.hom([], codomain=ambient))
    (identity,) = tuple(isotropic_stabilizer_right_cosets(category, marked))
    assert identity.domain() is ambient and identity.codomain() is ambient
    assert all(identity(point) == point for point in ambient.module_generators())
