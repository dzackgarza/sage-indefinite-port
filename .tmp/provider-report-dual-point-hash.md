# sage-extensions: elements of a dual module share one hash, so `finite_ordered_set` is O(n²) at ~17 ms per comparison

Reported against the revision the port pins, `sage-extensions` **`7b222e86214524a53802350f3cfda9c513707600`**
(`test: hold the isometry hash to equality and the subgroup to its generators`), installed at
`/repos/sage-extensions` in `ghcr.io/dzackgarza/sage:develop` @ revision `b0cd7d86`,
Sage 11.0.beta5 / CPython 3.14.10. Line numbers below are that revision's.

No port code and no polyhedral engine is needed to reproduce; the standard lattice `I6` suffices.

## Symptom

`finite_ordered_set` over finitely many **distinct** points of a `dual_module()` does not return in
practical time at the sizes polyhedral work produces. 570 distinct covectors did not finish in 560 s.
`RationalPolyhedralCone.__init__`, which calls it twice, inherited the cost: building one owned cone
from a 570-facet H-presentation did not finish in 779 s.

This is not a polyhedral cost. On the same cone and the same 202 rays, both exact engines answer in
milliseconds (measured separately, engine-only):

| engine | call | time | result |
| --- | --- | --- | --- |
| Normaliz | `NmzResult(cone, "SupportHyperplanes")` | 0.113 s | 570 facets, 0 equations |
| PPL | `Polyhedron(vertices=[0], rays=<202 rows>).inequalities()` | 0.052 s | 570 facets, 0 equations |
| PPL (reverse) | `Polyhedron(ieqs=<570 Normaliz rows>)` | 0.09 s | 112 extreme rays |

The two 570-facet families are the *same* rows: identical as raw signed primitive integer rows and
identical as sign-canonical sets, 570 distinct on each side.

## Bounded stack

`faulthandler.dump_traceback_later(90, repeat=True)` inside `_rational_polyhedral_cone`:

```
_rational_polyhedral_cone                 polyhedral_cones.py:865
  arrow_categories.object                 arrow_categories.py:1782, :968, :977
    objects._with_structure               objects.py:335
      owned_category._object_of           owned_category.py:1660
        RationalPolyhedralCone.__init__   polyhedral_cones.py:290   <- finite_ordered_set(...)
          finite_ordered_set              finite_ordered_sets.py:808
            Category.__call__             sage/categories/category.py:592
              FiniteOrderedSets._call_    finite_ordered_sets.py:292
                _on_points                finite_ordered_sets.py:333  <- dict.fromkeys(supplied)
                  SetMorphism.__eq__       set_categories.py:697
                    _compare_in_common_mor set_categories.py:707
                      ModuleMorphism._richcmp_  module_morphisms.py:1000
```

Line 333 is `distinct = tuple(dict.fromkeys(supplied))`. The hash table is not slow; it is
uninformative.

## Root cause

`hash()` of a point of `lattice.dual_module()` does not depend on the point. For 64 distinct dual
points of `I6` there is **1 distinct hash**; for 216 distinct dual points of a rank-6 Lorentzian
lattice, likewise **1**. Every insertion therefore collides with every extant key, `dict.fromkeys`
falls back to `==` once per pair, and each `==` is a `_richcmp_` of module morphisms
(`module_morphisms.py:1000`) costing ~17 ms.

Primal points are unaffected, which isolates the defect to the dual side:

| family (rank-6 Lorentzian lattice, same coordinates) | points | distinct hashes |
| --- | --- | --- |
| `lattice_` elements | 125 | 80 |
| `lattice_.dual_module()` elements | 125 | **1** |

The same collapse holds on the standard lattice: 64 distinct dual points of `I6` give 1 hash, and 216
distinct dual points of a rank-6 Lorentzian lattice also give 1 hash.

`dict.fromkeys` on the primal family is below timing resolution (0.001 s at n=64); on the dual family
it is quadratic:

| n | `dict.fromkeys`, dual | `finite_ordered_set`, dual |
| --- | --- | --- |
| 8 | 0.585 s | 0.647 s |
| 16 | 1.961 s | 2.446 s |
| 32 | 6.525 s | 9.527 s |
| 64 | 22.926 s | 49.805 s (larger family) |
| 96 | — | did not finish in ~200 s |
| 570 | ≈1 800 s extrapolated | did not finish in 560 s (observed) |

## What the contract should be

1. Distinct points must hash distinctly: `__hash__` is a function of the point, not of the parent,
   the domain, a shared defining morphism, or an identity that all points share.
2. Consequently `finite_ordered_set` over *n* distinct points is roughly *O(n)*, not *O(n²)* with a
   17 ms comparison. At n=570, seconds.
3. `RationalPolyhedralCone.__init__` given an H-presentation with hundreds of facets returns in
   seconds. It already does for cells with a small facet count — the port's own large-cell test
   builds an owned cone from Normaliz H-data for `Result_Enumeration_1` in 12 s — so the regression
   is the facet count crossing the hash-collision threshold, not the construction itself.
4. Hashing must stay consistent with the equality already in place: equal points equal hashes,
   unequal points need not differ, but at these sizes they must in practice.

## Reproducer (provider-only, no port, no engine)

```python
import time
from sage_extensions.all import Lattices, ZZ
from sage_extensions.categories.polyhedral_cones import finite_ordered_set

I6 = [[1 if i == j else 0 for j in range(6)] for i in range(6)]

lattice_ = Lattices(ZZ)(I6)
dual = lattice_.dual_module()
labels = tuple(dual.module_generating_set())
coordinates = [(a, b, c, 0, 0, 0) for a in range(1, 5) for b in range(0, 4) for c in range(0, 4)]
points = tuple(
    dual.linear_combination({l: ZZ(v) for l, v in zip(labels, row, strict=True) if v})
    for row in coordinates
)

print(len(dict.fromkeys(points)), "distinct points of", len(points))   # 64 of 64
print(len({hash(p) for p in points}), "distinct hashes")               # 1

for n in (8, 16, 32, 64):
    started = time.monotonic()
    dict.fromkeys(points[:n])
    print(n, f"{time.monotonic() - started:.3f}s")

for n in (8, 16, 32):
    started = time.monotonic()
    finite_ordered_set(points[:n])
    print(n, f"{time.monotonic() - started:.3f}s")
```

Expected output: 64 distinct points, 1 distinct hash, and the quadratic ladder above. To see the
primal contrast, build the same coordinate family in `lattice_` instead of `dual` and note 80
distinct hashes and sub-millisecond `dict.fromkeys`.

## Blast radius observed from the port

`Result_Enumeration_37` (total mode) is the specimen: a rank-6 perfect cell with 202 rays whose cone
has 570 facets. The port reaches the defect through the provider's public construction, from both
routing branches:

```python
# port: cell_cone(cell) -> RationalPolyhedralCones(lattice).from_rays(rays)
#   -> _cone_from_engine_polyhedron -> _rational_polyhedral_cone -> RationalPolyhedralCone.__init__

# port: cell_cone(cell) -> _rational_polyhedral_cone(lattice, normaliz_covectors, ...)
#   -> RationalPolyhedralCone.__init__
```

Both branches end in the same `finite_ordered_set` over ~570 dual covectors, so no routing choice
inside the port avoids it, and no change of routing is justified by this evidence.
