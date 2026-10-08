r"""Isomorphisms and automorphisms of marked indefinite lattices.

A marked lattice is a diagram of lattice morphisms ending at a lattice L
(references/port-design.md, "Marked lattices"):

- a pointed lattice ``<n> -> L``, an object of the coslice ``<n>/Lattices``;
- a sublattice ``I -> L`` or a vector configuration, an object of ``Ar(Lattices)``;
- a flag ``I_1 -> I_2 -> L``, an object of ``Fun([2], Lattices)``.

A plain lattice is the case with no marking. The research preamble's core homsets of
these categories call these three entry points when L is indefinite.
"""

from sage.categories.category import Category
from sage.categories.morphism import Morphism
from sage.structure.parent import Parent

from sage_indefinite_port.readiness import unfinished


@unfinished(18)
def isometry(source: Parent, target: Parent) -> Morphism | None:
    r"""An isomorphism of marked lattices ``source -> target``, or ``None`` when there is none."""


@unfinished(18)
def automorphism_group_generators(marked: Parent) -> Parent:
    r"""Generators of the automorphism group of a marked lattice, as a finite ordered set."""


@unfinished(18)
def primitive_marking_classes(category: Category, source: Parent | Morphism, lattice: Parent) -> Parent:
    r"""One marked lattice for each isomorphism class of primitive markings of ``lattice`` of shape ``source``.

    ``category`` is the category of the marking: the coslice under ``<n>`` for pointed
    lattices, ``Ar(Lattices)`` for sublattices, ``Fun([2], Lattices)`` for flags. ``source``
    is the shape: the lattice ``<n>``, a totally isotropic lattice of rank ``k`` (zero form),
    or the chain of such lattices for a flag.
    """


__all__ = [
    "automorphism_group_generators",
    "isometry",
    "primitive_marking_classes",
]
