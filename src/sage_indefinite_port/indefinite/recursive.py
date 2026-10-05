"""Public recursive indefinite-lattice operations on live preamble objects."""

from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.lattices import Lattices

from sage_indefinite_port.indefinite.lorentzian_cells import IndefiniteOrthogonalAlgorithm


def orthogonal_group_generators(homset) -> tuple[LatticeIsometryMethods, ...]:
    lattice = homset.domain()
    group = IndefiniteOrthogonalAlgorithm().orthogonal_group(lattice)
    return tuple(group.generators())


def isometry(source: Lattices.ParentMethods, target: Lattices.ParentMethods):
    return IndefiniteOrthogonalAlgorithm().isometry(source, target)


def vector_equivalence_witness(homset, left, right):
    if left.parent() is not homset.domain() or right.parent() is not homset.domain():
        raise ValueError("vector equivalence requires vectors in the homset lattice")
    return IndefiniteOrthogonalAlgorithm().vector_transporter(left, right)


def vector_stabilizer_generators(homset, element) -> tuple[LatticeIsometryMethods, ...]:
    if element.parent() is not homset.domain():
        raise ValueError("vector stabilizer requires an element of the homset lattice")
    subgroup = IndefiniteOrthogonalAlgorithm().vector_stabilizer(element)
    return tuple(subgroup.generators())


__all__ = [
    "isometry",
    "orthogonal_group_generators",
    "vector_equivalence_witness",
    "vector_stabilizer_generators",
]
