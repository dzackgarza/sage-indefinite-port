"""Resolve the two action-API calls on a tiny closed finite G-set.

Point set: all six 2-element subsets of {1,2,3,4}; S_4 acts transitively, so there is one
orbit and the stabilizer of a 2-subset is S_2 x S_2 of order 4.
"""
from itertools import combinations

from sage_extensions.all import OwnedGroups
from sage_extensions.categories.group.g_sets import FiniteGSets

symmetric = OwnedGroups().S(4)
points = tuple(frozenset(pair) for pair in combinations([1, 2, 3, 4], 2))
print("point set size:", len(points))


def action(element, point):
    return frozenset(element(x) for x in sorted(point))


g_set = FiniteGSets(symmetric)(points, action)
print("g_set:", g_set)

# (1) the G-set's point stabilizer -- called by groups.py and by orbit.stabilizer
point = points[0]
try:
    stab = g_set.stabilizer(point)
    print("g_set.stabilizer(point):", stab, "| type:", type(stab).__name__)
    for attr in ("cardinality", "group_generators", "natural_points", "action_on"):
        print(f"   has {attr}:", hasattr(stab, attr))
except Exception as error:
    print("g_set.stabilizer(point) raised:", type(error).__name__, str(error)[:140])

# (2) one representative per orbit
orbits = g_set.orbits()
print("orbit set:", orbits, "| iterable:", hasattr(orbits, "__iter__"))
listed = list(orbits)
print("orbit count:", len(listed))
for orbit in listed:
    print("   orbit representative:", orbit.representative(),
          "| |orbit|:", orbit.points().cardinality())
    try:
        print("   orbit.stabilizer() ->", orbit.stabilizer(), "| order:", orbit.stabilizer().cardinality())
    except Exception as error:
        print("   orbit.stabilizer() raised:", type(error).__name__, str(error)[:120])

print("\n=== (3) re-forming the stabilizer inside `symmetric` from its own generators ===")
stabilizer = g_set.stabilizer(point)
generators = tuple(stabilizer.group_generators())
print("stabilizer generators:", generators)
regrouped = symmetric.subgroup(generators)
print("regrouped:", regrouped)
print("   cardinality:", regrouped.cardinality())
print("   has natural_points:", hasattr(regrouped, "natural_points"))
print("   natural_points:", regrouped.natural_points() if hasattr(regrouped, "natural_points") else "n/a")
try:
    act = regrouped.action_on(regrouped.natural_points())
    print("   action_on natural points works:", act)
except Exception as error:
    print("   action_on raised:", type(error).__name__, str(error)[:110])
# and the facet action on the same closed set, restricted to the stabilizer
sub_g_set = FiniteGSets(regrouped)(points, action)
print("   facets under the regrouped stabilizer:", sub_g_set)
print("   orbits:", len(list(sub_g_set.orbits())))
