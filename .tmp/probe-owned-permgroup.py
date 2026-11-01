"""Can an owned permutation group be built from explicit permutations, and act on facets?

Lookup, not assertion: enumerate what the provider offers, then try the construction the
cell needs (a group on the ray positions, acting on *sets* of rays).
"""
from sage.all import Permutation
from sage.libs.gap.libgap import libgap

from sage_extensions.all import OwnedGroups
from sage_extensions.categories.group.g_sets import FiniteGSets

symmetric = OwnedGroups().S(4)
print("S(4):", symmetric)
print("natural_points:", symmetric.natural_points())
methods = [m for m in dir(symmetric) if not m.startswith("_")]
print("public methods:", methods[:14])

# element from a cycle tuple, as the provider's own tests do: group((1, 2))
try:
    g = symmetric((1, 2))
    print("symmetric((1,2)) ->", g, "| sign:", g.sign())
except Exception as error:
    print("symmetric((1,2)) raised:", type(error).__name__, str(error)[:120])

# build a subgroup from explicit permutations given as one-line images
try:
    p = Permutation([2, 1, 3, 4])
    elt = symmetric(*p.cycle_tuples()) if p.cycle_tuples() else symmetric(())
    print("element from cycle_tuples:", elt)
    subgroup = symmetric.subgroup((elt,))
    print("subgroup order:", subgroup.cardinality())
except Exception as error:
    print("subgroup construction raised:", type(error).__name__, str(error)[:160])

# acting on SETS of points: action_on only accepts the natural domain
facets = ((1, 2), (3, 4))
try:
    g_set = FiniteGSets(symmetric)(facets, lambda element, facet: tuple(sorted(element(point) for point in facet)))
    print("explicit-action GSet built:", g_set)
    print("orbits:", g_set.orbits().cardinality())
except Exception as error:
    print("FiniteGSets(...) raised:", type(error).__name__, str(error)[:160])
