r"""Orbits, equivalence witnesses and stabilizers of primitive isotropic sublattices and flags.

These are the six entry points the research preamble's isotropic-orbit methods import for
an indefinite lattice ``L`` with orthogonal group ``G = L.O()``. Sublattices are preamble
subobjects of ``L``; flags are preamble ``IsotropicFlag`` objects. Ranks are at most 2
(map section 13.5).
"""

from dzack_research.preamble.categories.isotropic_orbits import IsotropicFlag
from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods, LatticeIsometryMor
from dzack_research.preamble.categories.rings.ring_foundation import _OwnedRingElement
from sage.structure.parent import Parent

from sage_indefinite_port.readiness import unfinished


@unfinished(18)
def isotropic_sublattice_orbit_representatives(group: LatticeIsometryMor, rank: _OwnedRingElement) -> Parent:
    r"""One primitive totally isotropic sublattice of the given rank from each ``G``-orbit."""


@unfinished(18)
def isotropic_flag_orbit_representatives(group: LatticeIsometryMor, rank: _OwnedRingElement) -> Parent:
    r"""One isotropic flag of the given top rank from each ``G``-orbit."""


@unfinished(18)
def isotropic_sublattice_equivalence_witness(group: LatticeIsometryMor, left: Parent, right: Parent) -> LatticeIsometryMethods | None:
    r"""An element of ``G`` carrying ``left`` onto ``right``, or ``None`` when they lie in different orbits."""


@unfinished(18)
def isotropic_flag_equivalence_witness(group: LatticeIsometryMor, left: IsotropicFlag, right: IsotropicFlag) -> LatticeIsometryMethods | None:
    r"""An element of ``G`` carrying the flag ``left`` onto ``right`` term by term, or ``None``."""


@unfinished(18)
def isotropic_sublattice_stabilizer_generators(group: LatticeIsometryMor, sublattice: Parent) -> Parent:
    r"""Generators of the setwise stabilizer of ``sublattice`` in ``G``."""


@unfinished(18)
def isotropic_flag_stabilizer_generators(group: LatticeIsometryMor, flag: IsotropicFlag) -> Parent:
    r"""Generators of the stabilizer of ``flag`` in ``G``, preserving every term."""


__all__ = [
    "isotropic_flag_equivalence_witness",
    "isotropic_flag_orbit_representatives",
    "isotropic_flag_stabilizer_generators",
    "isotropic_sublattice_equivalence_witness",
    "isotropic_sublattice_orbit_representatives",
    "isotropic_sublattice_stabilizer_generators",
]
