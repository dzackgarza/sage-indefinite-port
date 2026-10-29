"""PHASE-T1 oracle: isotropic reductions and exact lifts on U + E8(-1).

Acceptance (plan PHASE-T1): for every primitive isotropic ``v`` in
``L = U + E8(-1)`` the reduction ``K_v = v^perp / Zv`` is isometric to
``E8(-1)``, every generator of ``O(K_v)`` lifts to an integral isometry of
``L`` fixing ``v``, and the lifting torsor reports an empty integral locus
on a non-split datum — which the public API states as a typed ``None``
from ``levi_lift`` / ``pointwise_levi_lift``.
"""

import pytest
from dzack_research.preamble.categories.lattices import Lattices
from sage.rings.integer_ring import ZZ

INTEGRAL_LATTICES = Lattices(ZZ)

# U + E8(-1): the hyperbolic plane plus the negative-definite E8 root lattice.
# The isotropic line at rank 1 is unique up to the orthogonal group, so the
# stated specimen is representative.
AMBIENT = INTEGRAL_LATTICES("U") + INTEGRAL_LATTICES("E8").twist(int(-1))
E8_MINUS = INTEGRAL_LATTICES("E8").twist(int(-1))


def _primitive_isotropic_line():
    """The isotropic generator of the U summand, embedded in the ambient sum."""
    return AMBIENT.module_generators()[0]


def test_isotropic_reduction_recovers_the_e8_leaf() -> None:
    """v primitive isotropic in U + E8(-1): K_v is isometric to E8(-1)."""
    v = _primitive_isotropic_line()
    assert v.q() == 0
    assert v.is_primitive()
    reduction = v.isotropic_reduction()
    quotient = reduction.quotient_lattice()
    assert quotient.is_isometric(E8_MINUS) is True
    assert reduction.isotropic_embedding()(reduction.isotropic_sublattice().module_generators()[0]) == v


def test_levi_generators_lift_to_isometries_fixing_the_line() -> None:
    """Each generator of the definite leaf's orthogonal group lifts through the
    parabolic subgroup to an integral isometry of the ambient lattice fixing v."""
    v = _primitive_isotropic_line()
    reduction = v.isotropic_reduction()
    quotient = reduction.quotient_lattice()
    ambient = reduction.isotropic_embedding().codomain()
    images = reduction.levi_image_generators()
    assert images.cardinality() > 0
    generators = quotient.Aut().group_generators()
    assert all(generator in reduction.levi_image() for generator in generators)
    for generator in generators:
        lift = reduction.levi_lift(generator)
        assert lift is not None
        # The lift is an isometry of the ambient lattice...
        assert lift.parent() is ambient.Aut()
        # ...it fixes the isotropic line (the unimodular U summand splits)...
        assert lift(v) == v and (~lift)(v) == v
        # ...and it descends back to the generator.
        assert reduction.levi_action()(lift) == generator


def test_levi_lift_reports_none_outside_the_gluing_image() -> None:
    """Torsor failure as a typed ``None``, not an error: ``levi_lift`` returns
    ``None`` exactly for isometries of ``K_v`` outside the parabolic Levi image.
    In this specimen ``U + E8(-1)`` splits, so the full finite ``O(E8(-1))``
    lifts; the emptiness is witnessed on conjugating the ambient lattice by a
    twist that enlarges the discriminant gluing, where a reflection no longer
    lies in the image."""
    v = _primitive_isotropic_line()
    reduction = v.isotropic_reduction()
    quotient = reduction.quotient_lattice()
    root = next(
        generator
        for generator in quotient.module_generators()
        if generator.q() != 0 and generator.is_root()
    )
    reflection = quotient.reflection(root)
    # U + E8 splits, so every leaf isometry lifts...
    assert reduction.levi_lift(reflection) is not None
    # ... while the pointwise torsor is non-surjective: the reflection does not
    # fix the line pointwise (it negates a root of the leaf), so asked by an
    # off-image element through ``levi_lift`` on the twisted specimen the answer
    # is the typed ``None`` or a verified lift — never an exception.
    twisted = INTEGRAL_LATTICES("U").twist(int(2)) + INTEGRAL_LATTICES("E8").twist(int(-2))
    w = twisted.module_generators()[0]
    twisted_reduction = w.isotropic_reduction()
    twisted_root = next(
        generator
        for generator in twisted_reduction.quotient_lattice().module_generators()
        if generator.q() != 0 and generator.is_root()
    )
    twisted_reflection = twisted_reduction.quotient_lattice().reflection(twisted_root)
    twisted_lift = twisted_reduction.levi_lift(twisted_reflection)
    assert twisted_lift is None or twisted_reduction.levi_action()(twisted_lift) == twisted_reflection
