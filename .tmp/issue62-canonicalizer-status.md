# Issue #62: canonical representatives and duplicate-free orbit enumeration — what is implemented

Status: **implemented on the supported domain**; #62 itself stays open (it asks for the general
effectively-enumerable integral `Γ`-set, which is broader than what `pointed_isometry` decides).

## What was already there, and what was missing

`recursive.py` already had `marked_orbit_action` (the retained G-set), `marked_orbit` (one orbit as a lazy
predicate set), `marked_square_slice` (`L[n]` as a predicate set), `marked_orbit_quotient` (`L[n]/O(L)` as the
provider's `_OrbitSet`) and `marked_orbit_representatives` (a height-ordered walk that dedups against emitted
representatives), plus `_coordinate_box` — the well-order with finite initial segments that #62 part A asks for.
There was **no canonicalizer**: nothing turned a single vector into the least member of its orbit.

## What was added

* `marked_orbit_canonical(vector) -> MarkedOrbitCanonical(representative, transporter, height)` — #62 part C. It
  scans only the finite initial segment up to the vector's own height and returns the first equivalent element,
  which is the global minimum because any smaller equivalent would be inside the segment. Termination is
  guaranteed by the vector itself lying in its own box. The transporter satisfies
  `transporter(representative) == vector`, and `transporter` is a real isometry, not a coordinate relabelling.
* `marked_orbit_is_canonical(vector)` — the derived predicate, exact.
* `marked_orbit_enumeration(lattice, square, limit=None, height_bound=None)` — #62 parts B and D. It yields one
  `MarkedOrbitDecision(candidate, emitted, representative, transporter)` per square-`square` candidate: `emitted`
  exactly for the first member of each orbit in the declared order, with the isometry carrying the retained
  representative onto the candidate in every case. Because an emitted vector is its orbit's minimum, two
  emissions cannot be equivalent — duplicate-freedom is a consequence of minimality, not of remembering a list.
  `height_bound` bounds the **search** (bounded-height classification); the equivalence test stays global, so a
  transporter may leave the box and classes are never split by it.
* `marked_orbit_representatives` now delegates to that enumeration, so the two surfaces cannot drift. Its
  semantics are unchanged: same order, same emissions, `limit` still a cap on yields and not a termination claim.

## Contract kept

* **Declared-order relative.** The representative is canonical relative to the declared basis and coordinate
  order. `test_the_quotient_is_unchanged_under_a_permuted_presentation` checks exactly that split: coordinates of
  the representative depend on the presentation, while the orbit relation — hence the quotient — does not.
* **No completion signal.** A finite quotient stops producing representatives while the search continues; absence
  of further output is not a certificate. Tested as prefix stability: asking for fewer representatives gives a
  prefix of the longer request, and nothing in the protocol reports exhaustion.
* **Unknown is never inequivalence.** `UnfinishedCapability` from the backend propagates out of both surfaces.
  Reading it as "distinct" would emit one orbit twice, which #62 forbids.
* **No parallel framework.** Everything is on `EichlerRecursion` with `pointed_isometry` as the equivalence and
  the lattice's own `linear_combination` / `zero` / `q` / `content` for elements. Nothing generic about groups or
  sets was invented.

## Two findings from the work, both kept as tests

1. **The undecided subdomain is non-primitive radical vectors, and only those.** Measured on `L = <2> + 0`:
   `(0, 1)` versus `(0, -1)` gets a witness, `(0, 1)` versus `(0, 2)` and `(0, 3)` versus `(0, 1)` correctly
   return `None`, but `(0, 2)` versus `(0, -2)` reaches `primitive_isotropic_subobject`, which asserts that the
   span of `(2 e_1)` has torsion and must be saturated. `pointed_isometry` routes every square-zero vector to
   the isotropic-**line** machinery, which needs the line primitive and non-degenerate, and unlike
   `pointed_automorphism_generators` it does not first detour through `vector_lies_in_radical`.
   `_require_decidable_pointed_domain` therefore refuses exactly "radical with content > 1" and names it. This
   is also why #62 test 4's content family (`c e_1`) is not reachable through this route yet: those are radical
   vectors once `q` vanishes on them.
2. **`q` of the hyperbolic plane is `2xy`, not `xy`.** `b(v, v)` for `[[0,1],[1,0]]` is `2xy`, so square 2 is
   `{(±1, ±1)}` and square 4 is the `(±1, ±2)` / `(±2, ±1)` orbit. My first expectations were written against
   `xy` and were wrong; the tests now derive their expected minima from the orbit description and check
   minimality independently, rather than guessing coordinates.

## Cost

`tests/test_marked_orbit_canonical.py`: 9 tests, 75 s, peak tree RSS 3.75 GiB (Julia baseline 3.15 GiB, python
590 MiB). The cost is the algorithm itself — every candidate is compared to each emitted representative through
`pointed_isometry` — not test overhead. `tests/test_degenerate_pointed_markings.py` still passes 39/39, unchanged.

## Not claimed

The general #62 iterator for an arbitrary effectively-enumerable integral `Γ`-set: `pointed_isometry` decides
`O(L)`-transport of vectors, and that is the only equivalence wired in. An infinite quotient is not demonstrated
in-domain here, because the families that give one on these lattices (unbounded content) are exactly the
non-primitive radical vectors the route refuses.
