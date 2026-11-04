"""The fifth site: transporting a *tuple of subsets*, as libgap.RepresentativeAction does."""
from itertools import combinations

from sage_extensions.all import OwnedGroups
from sage_extensions.categories.group.g_sets import FiniteGSets

symmetric = OwnedGroups().S(4)
pairs = [frozenset(p) for p in combinations([1, 2, 3, 4], 2)]
# points are tuples of two subsets, as `chain` and `image` are tuples of facets
tuples = tuple((a, b) for a in pairs for b in pairs if a != b)


def action(element, point):
    return tuple(frozenset(element(x) for x in sorted(subset)) for subset in point)


g_set = FiniteGSets(symmetric)(tuples, action)
print("tuple G-set size:", len(tuples))

a = ((frozenset({1, 2}), frozenset({3, 4})))
b = ((frozenset({1, 3}), frozenset({2, 4})))   # S_4 moves a to b
try:
    witness = g_set.transporter_witness(a, b)
    print("transporter_witness(a, b) ->", witness)
    print("   carries a to b:", action(witness, a) == b)
except Exception as error:
    print("transporter_witness raised:", type(error).__name__, str(error)[:130])

# a relatedness decision for a pair that cannot be related (sizes differ)
c = (frozenset({1, 2}), frozenset({1, 3}))
try:
    print("transporter_witness(a, c) ->", g_set.transporter_witness(a, c))
except Exception as error:
    print("transporter_witness(a, c) raised:", type(error).__name__, str(error)[:110])
print("orbit relation a ~ b:", g_set.orbits().orbit_of(a) == g_set.orbits().orbit_of(b))
