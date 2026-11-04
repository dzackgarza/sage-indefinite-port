# Draft: replace `_CellFaces`' direct GAP calls with the owned construction

**Status: READY TO APPLY at a clean source boundary. Verified by interface exercise; not applied,
because `/work` is bound live into a running acceptance sweep.**

Provider `c546077d`. Guard already in place: `tests/test_cell_faces_facet_orbits.py` (5 tests, 1.6 s).
Real gate: A_PERFECT's 92 cases.

## Verified semantics (exercised, not read-off-guesswork)

Test fixture: `S_4` acting on the six 2-element subsets of `{1,2,3,4}` — a closed finite set, one orbit,
stabilizer of a pair of order 4, and under that stabilizer the pairs split into 3 orbits.

| question | exact answer |
| --- | --- |
| build the ambient group | `symmetric = OwnedGroups().S(n)` |
| element from one-line images | `symmetric([tuple(cycle) for cycle in Permutation([i+1 for i in images]).cycle_tuples()])` — a **list** of cycles; `symmetric(c1, c2)` raises `TypeError` |
| subgroup | `symmetric.subgroup(generators)` |
| explicit action on an arbitrary finite set | `FiniteGSets(group)(points, action)` where `action(element, point)` returns the image |
| **one representative per orbit** | `g_set.orbits()` returns an iterable orbit set; iterating yields one orbit per class, each with `.representative()`, `.points()`, `.stabilizer()` |
| **the stabilizer subgroup retaining the action** | `g_set.stabilizer(point)` returns an owned subgroup with `group_generators()` and `cardinality()`, but **without** `natural_points()`/`action_on()`. Re-form it inside the ambient group: `symmetric.subgroup(tuple(stab.group_generators()))` — that object *does* have `natural_points()` and `action_on(...)`, and `FiniteGSets` can be built on it |

Measured on the fixture: orbit count 1, representative `frozenset({1, 2})`, `|orbit| = 6`, stabilizer order 4,
and 3 orbits of pairs under the regrouped stabilizer — the values the mathematics requires
(`S_2 x S_2` fixes `{1,2}` and `{3,4}` and moves the other four).

One detail to carry: `stabilizer.group_generators()` can include identity elements (the fixture returns
`1` several times). That is harmless for `subgroup`, but a caller must not assume they are non-trivial.

No new carrier is introduced: `symmetric.subgroup(...)` is the same constructor the cell group already uses.

## Sites to change

1. `_CellFaces.__init__`, `lorentzian_cells.py:3568` — `self.group = libgap.Group([...], libgap.eval("()"))`
2. `_CellFaces.flag_simplices`, `:3603` — `libgap.Orbits(stabilizer, facets, OnSets)` and
   `libgap.Stabilizer(stabilizer, representative, OnSets)`
3. `MarkedCellOrbitAlgorithm._timelike_global_orbits`, `:3458` — `libgap.Stabilizer(source.group, sorted(...), OnSets)`

## Shape of the replacement

```python
from sage.all import Permutation
from sage_extensions.all import OwnedGroups
from sage_extensions.categories.group.g_sets import FiniteGSets


def _owned_permutation(symmetric, images):
    cycles = Permutation([image + 1 for image in images]).cycle_tuples()
    return symmetric([tuple(cycle) for cycle in cycles]) if cycles else symmetric([])


def _subset_action(element, subset):
    """The action on sets of ray positions; `action_on` will not infer it."""
    return frozenset(element(point) for point in subset)


# __init__: the cell group
symmetric = OwnedGroups().S(len(self.rows))
self.group = symmetric.subgroup(
    tuple(_owned_permutation(symmetric, images)
          for images in backend.point_permutations_from_rows(cell, cell.rows))
)
self._facet_g_set = FiniteGSets(self.group)(
    tuple(sorted(facet) for facet in self.facet_incidences), _subset_action
)

# flag_simplices: descend one rank using representatives and a regrouped stabilizer
g_set = FiniteGSets(group)(tuple(sorted(facet + 1 for facet in f) for f in facets), _subset_action)
for orbit in g_set.orbits():
    representative = orbit.representative()
    subgroup = symmetric.subgroup(tuple(orbit.stabilizer().group_generators()))
    descend(..., subgroup)
```

`flag_simplices`' current signature passes a bare `group` (a GAP object) down the recursion; the
replacement must thread the *owned* subgroup plus the ambient `symmetric` so the regrouping is always
available. That is the one structural edit beyond the call substitution, and it is why the guard and the
A_PERFECT sweep are both required.

## Boundary order

1. Apply the substitution.
2. `tests/test_cell_faces_facet_orbits.py` must stay **5/5** unchanged.
3. A_PERFECT per-case sweep as the real gate; a regression means revert, not patch around.
4. Only then report the GAP calls as removed.
