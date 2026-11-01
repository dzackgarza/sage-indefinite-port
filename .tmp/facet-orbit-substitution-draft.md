# Unapplied draft: replace `_CellFaces`' direct GAP calls with the owned construction

**Status: DRAFT, NOT APPLIED, NOT RUN.** Prepared for the next clean source boundary. It must not be
applied while an acceptance sweep is running, because `/work` is bound live into it and an edit would
change the code under the running cases.

Verified prerequisites (from the lookup on `c546077d`, exercised):

- `OwnedGroups().S(n)` is the owned symmetric group on `n` points.
- An element comes from a **list** of cycles: `symmetric([(1, 3), (2, 4)])` → `(1 3)(2 4)`.
  `symmetric(c1, c2)` raises `TypeError`, so a list is mandatory.
- `symmetric.subgroup((element,))` builds the subgroup.
- `FiniteGSets(group)(points, action)` takes an **explicit action** on an arbitrary finite point set,
  which is precisely what `PermutationGroups.ParentMethods.action_on` refuses to infer ("specify a
  separate group action for another point set"). A `ValueError` results if the point set is not closed
  under the action.
- `PermutationGroups.ParentMethods.orbit(point)` and `stabilizer(point)` are defined through
  `natural_g_set()`, i.e. they act on natural points, not on sets of points.

Guard already in place: `tests/test_cell_faces_facet_orbits.py` pins the result — 4 facets, 4 distinct
incidence sets, group order 8, one facet orbit of size 4, facet stabilizer order 2, one flag — on
`Result_Enumeration_0`. The real gate is A_PERFECT's 92 cases.

## Sites to change

1. `_CellFaces.__init__`, `lorentzian_cells.py:3568` — `self.group = libgap.Group([PermList...], ())`
2. `_CellFaces.flag_simplices`, `:3603` — `libgap.Orbits(stabilizer, facets, OnSets)` and
   `libgap.Stabilizer(stabilizer, representative, OnSets)`
3. `MarkedCellOrbitAlgorithm._timelike_global_orbits`, `:3458` —
   `libgap.Stabilizer(source.group, sorted(...), OnSets)`

## Draft

```python
from sage.all import Permutation
from sage_extensions.all import OwnedGroups
from sage_extensions.categories.group.g_sets import FiniteGSets


def _owned_permutation(symmetric, images):
    """The owned element of ``symmetric`` with the given one-line images (0-based)."""
    cycles = Permutation([image + 1 for image in images]).cycle_tuples()
    return symmetric([tuple(cycle) for cycle in cycles]) if cycles else symmetric([])


def _facet_action(element, facet):
    """The action on *sets* of ray positions, which the provider will not infer."""
    return frozenset(element(point) for point in facet)


# in _CellFaces.__init__, replacing libgap.Group(...)
symmetric = OwnedGroups().S(len(self.rows))
generators = tuple(
    _owned_permutation(symmetric, images)
    for images in backend.point_permutations_from_rows(cell, cell.rows)
)
self.group = symmetric.subgroup(generators)
self._facet_g_set = FiniteGSets(self.group)(
    finite_ordered_set(tuple(sorted(facet) for facet in self.facet_incidences)),
    _facet_action,
)
```

### Two calls I have not resolved, and will not guess

- **Enumerating orbit representatives** from `self._facet_g_set.orbits()`. `groups.py` shows
  `action.orbits().orbit_of(point).points()` for a single orbit, so `orbit_of` exists; the accessor that
  yields one representative per orbit is not established from source.
- **The stabilizer of a facet** as a *subgroup of `self.group`*: `g_set.stabilizer(point)` returns a
  stabilizer of the G-set, but `flag_simplices` needs a group it can act on the next rank of facets with,
  i.e. a subgroup carrying the same permutation action on points.

Both are lookups of the same kind that resolved the earlier questions, and both must be read before the
patch is applied. Guessing them here would be exactly the unverified substitution I have refused twice.

## Order of work at the boundary

1. Resolve the two calls above against `g_sets.py`.
2. Apply the patch; run `tests/test_cell_faces_facet_orbits.py` (1.6 s) — it must stay 5/5 unchanged.
3. Run the A_PERFECT per-case sweep as the real gate; a regression there means revert, not patch around.
4. Only then report the GAP calls as removed from `lorentzian_cells.py`.
