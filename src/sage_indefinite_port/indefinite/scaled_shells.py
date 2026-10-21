"""First-success multiplier searches over exact definite close-vector shells."""

from __future__ import annotations

from dzack_research.preamble.all import NN, Lattices
from dzack_research.preamble.categories.sets.set_categories import (
    CartesianProductsOfSets,
    Sets,
)
from sage.structure.element import RingElement


def first_scaled_shell(
    lattice: Lattices.ParentMethods, target: Lattices.ElementMethods, square_bound: RingElement, max_multiplier: int, *, exact_distance: bool = False
) -> CartesianProductsOfSets.ElementMethods | None:
    """Return the least positive multiplier with a nonempty shell or sphere.

    Fixed-target distance enumeration belongs to the definite lattice. The
    ordered multiplier search and the selection of its first success live here.
    """
    family = lattice.scaled_close_vector_shells(target, square_bound, NN(max_multiplier))
    subsets = lattice.finite_subsets()
    field = lattice.base_ring().fraction_field()
    for multiplier in family.index_set():
        distances = family.value(multiplier)
        boundary = field(int(multiplier) ** 2) * field(square_bound)
        vectors = tuple(point for point in distances.index_set() if not exact_distance or distances.value(point) == boundary)
        if vectors:
            return Sets().product((NN, subsets))((multiplier, subsets(vectors)))
    return None
