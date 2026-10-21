"""First-success multiplier searches over exact definite close-vector shells."""

from __future__ import annotations

from dzack_research.preamble.all import NN, Lattices, ModulesOverIntegralDomains
from dzack_research.preamble.categories.sets.set_categories import (
    CartesianProductsOfSets,
    Sets,
)
from sage.structure.element import RingElement


def first_scaled_shell(
    lattice: Lattices.ParentMethods, target: Lattices.ElementMethods, square_bound: RingElement, max_multiplier: int, *, exact_distance: bool = False
) -> CartesianProductsOfSets.ElementMethods | None:
    r"""Return the least multiplier ``1 <= m <= max_multiplier`` with a nonempty shell or sphere.

    # source: Lorentzian.g:359-393 LORENTZ_FindPositiveVectors (OnlyShortest)

    For ``eVal = 1, 2, ...`` the source scales the affine target by ``eVal``
    and its square distance by ``eVal^2``, enumerates the fixed-target close
    vectors, and returns at the first ``eVal`` with a solution. The scaling,
    the ordered multiplier loop and the first-success stop are this port's;
    the fixed-target enumeration ``close_vectors(target, bound)`` is the
    definite lattice's. ``exact_distance`` keeps only the boundary
    ``q(x - m*target) = m^2*square_bound`` (the source's ``"isotrop"`` filter
    ``fNorm = 0``).
    """
    ambient = lattice.vector_space()
    if target.parent() is lattice:
        target = ModulesOverIntegralDomains(lattice.base_ring())(lattice).generic_fibre_map()(target).underlying_element()
    if target.parent() is not ambient:
        raise ValueError("the affine target must belong to the lattice or its rational span")
    field = ambient.base_ring()
    bound = field(square_bound)
    subsets = lattice.finite_subsets()
    for multiplier in range(1, int(max_multiplier) + 1):
        boundary = field(multiplier**2) * bound
        distances = lattice.close_vectors(ambient.scalar_multiple(field(multiplier), target), boundary)
        vectors = tuple(point for point in distances.index_set() if not exact_distance or distances.value(point) == boundary)
        if vectors:
            return Sets().product((NN, subsets))((NN(multiplier), subsets(vectors)))
    return None
