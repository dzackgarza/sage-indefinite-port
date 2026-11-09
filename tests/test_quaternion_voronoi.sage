r"""Voronoi generators of the units of ZZ<1,i,j,k> in (-1,3).

The generators are units of the order. Every unit found in a coordinate box is written
by descent on phi_S as a word ending in a word of length at most four. This is
evidence of completeness on the box, and it does not prove generation.
"""

from itertools import product

from sage_extensions.all import *
from sage_indefinite_port.indefinite.quaternion_voronoi import QuaternionVoronoi


def test_the_order_has_one_perfect_form_class_and_integral_generators():
    model = QuaternionVoronoi(int(-1), int(3))
    forms = model.perfect_forms()
    generators = model.unit_generators()
    assert len(forms) == int(1)
    assert all(abs(g.reduced_norm()) == int(1) and all(c.denominator() == int(1) for c in g) for g in generators)
    assert any(g.reduced_norm() == int(-1) for g in generators)


def test_units_in_a_box_are_words_in_the_generators():
    model = QuaternionVoronoi(int(-1), int(3))
    algebra = model._algebra
    basis = model._basis
    generators = model.unit_generators()
    moves = [*generators, *(~g for g in generators)]
    terminal = set([algebra.one()])
    frontier = [algebra.one()]
    for _ in range(int(4)):
        frontier = [move * x for x in frontier for move in moves if move * x not in terminal]
        terminal.update(frontier)
    span = range(int(-5), int(6))
    units = [sum(int(c) * q for c, q in zip(point, basis, strict=True)) for point in product(span, repeat=int(4))]
    units = [u for u in units if abs(u.reduced_norm()) == int(1)]
    assert len(units) > int(20)
    assert algebra.one() in units
    for unit in units:
        assert model.reduce_unit(unit, generators) in terminal


def test_the_cell_area_equals_the_mass_formula():
    r"""The maximal order of discriminant 6 has O^1 area (pi/3)(2-1)(3-1) = 2 pi / 3.

    Lambda has index 2 in it, and locally at 2 the unit index is |F_4^x / F_2^x| = 3,
    so Lambda^1 has area 2 pi. A unit of norm -1 reverses orientation, halving it to pi.
    """
    from math import pi

    assert abs(QuaternionVoronoi(int(-1), int(3)).covolume() - pi) < int(1) / float(int(10) ** int(9))


def test_the_norm_one_extraction_applies_to_the_voronoi_generators():
    generators, kernel = QuaternionVoronoi(int(-1), int(3)).norm_one_generators()
    assert all(g.reduced_norm() == int(1) for g in generators)
    assert len(kernel) == int(1)


def test_the_full_orthogonal_group_generators_of_the_anisotropic_lattice():
    r"""The generators are integral isometries and include both determinants."""
    from sage_indefinite_port.indefinite.quaternion_lattices import QuaternionNormLattice
    from sage_indefinite_port.indefinite.quaternion_voronoi import orthogonal_group_generators

    model = QuaternionNormLattice(int(-1), int(3), int(2))
    generators = orthogonal_group_generators(model, QuaternionVoronoi(int(-1), int(3)), model.normalizer_classes_in_box(int(3)))
    points = tuple(model.lattice().module_generators())
    for generator in generators:
        for x in points:
            for y in points:
                assert model.lattice().b(generator(x), generator(y)) == model.lattice().b(x, y)
    assert generators[-1](points[1]) == -points[1]
    assert len(generators) > int(12)

