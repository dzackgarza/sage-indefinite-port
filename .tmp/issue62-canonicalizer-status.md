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

## The non-primitive radical gap is closed

The first version refused every non-primitive radical vector, because `pointed_isometry` routes square-zero
vectors to the isotropic-**line** machinery, which asserts that the span of `(2 e_1)` has torsion. That left the
unbounded-content family unreachable. It is now decided by `_radical_pointed_isometry`, at the existing owning
algorithm and through the retained radical splitting:

* `L` degenerate reduces (`indefinite_reduced_presentation`) to `C + R`, where every isometry is
  `[[A, 0], [B, D]]` with `A` in `O(C)`, `D` in `GL(R)` and `B` in `Hom(C, R)` arbitrary -- row convention, so
  `(0, r) M = (r B, r D)`. A radical marking is `(0, r)`, and `B = 0`, `A = 1` leaves exactly `D r = s`;
* `GL_k(Z)` acts on a column by its **content**, and is transitive on primitive vectors, so `r ~ s` iff
  `gcd(r) = gcd(s)`. The zero vector is the content-0 case and is fixed by everything, so it is never conflated
  with a nonzero class;
* the witness is constructed, not asserted: `V_r r_0 = e_1` and `V_s s_0 = e_1` from Smith normal form of the
  primitive rows, `D = V_r V_s^{-1}` unimodular, and `D r = s` is checked as it is built. The block is conjugated
  back to `L` with the retained `_conjugated`, and the result is verified to carry the marking onto its image.
* everything happens in the **reduced** presentation, so a permuted or mixed basis needs no special case: the
  radical is never assumed to occupy the last coordinates of the *given* basis. Tested on `<2> + 0`, on the
  radical-first presentation `[[0,0],[0,2]]`, on the mixed basis `[[2,2],[2,2]]` whose radical is spanned by
  `(1, -1)`, and at **radical rank 2** on `<2> + 0 + 0` (radical `span(e_1, e_2)`) and on `[[0,0,0],[0,0,0],[0,0,2]]`
  (radical first).
* integral invertibility is asserted where it is used: `|det| = 1` for the radical block and for the returned
  isometry, in the port's column convention, and the integral inverse is checked to have unit determinant;
* at radical rank 2 the transport `(0,1,0) -> (0,0,1)` **requires** an off-diagonal unimodular entry — no diagonal
  block can do it — so that specimen is the one that proves a genuine `GL_2(Z)` block is being constructed rather
  than a scaling. Pairing preservation is checked against a spread of arbitrary lattice vectors, not only against
  the marked one, which is what makes the result an isometry rather than a vector-matching map.

### Convention correction

The first version of this note had the block layout **transposed**: it wrote `(0, r) M = (r B, r D)` for
`M = [[A, 0],[B, D]]` while computing with row multiplication, and concluded that an arbitrary `B` could push a
radical vector out of the radical. Both are wrong. The retained form (`_degenerate_group_generators`, and the shear
`block[row, nonsingular_rank] = 1` in the code) is

    M = [[A, B],
         [0, D]]        (row blocks; the zero is bottom-left, B upper-right)

and in the row convention `(x, y) M = (x A, x B + y D)`. So a radical vector `(0, r)` maps to `(0, r D)`: it
**stays in the radical**, as it must, since the radical is characteristic — an isometry preserves the pairing and
hence `{v : b(v, L) = 0}`. `B` moves the nondegenerate part into the radical direction (the shear) and is
irrelevant to radical transport, which is decided by `D` alone. The implemented witness takes `A = 1`, `B = 0`:
that is a *choice* for simplicity, not a requirement, and it is now judged as such rather than being justified by a
false general statement.

`UnfinishedCapability` survives where support is genuinely not proved: cross-lattice radical transport between two
distinct loci would need an isometry between their reduced presentations carrying one radical coordinate onto the
other, and that is refused rather than guessed.

## Two findings from the work, kept

1. **The hyperbolic plane costs about 55 s per non-integral pointed transport, and that regression is kept.** `q`
   of `U = [[0,1],[1,0]]` is `2xy` (not `xy`), and transporting `(2, 1)` onto `(-2, -1)` does not take the
   already-integral branch: it reaches `_lattice_action(target).transporter(image)`, the correction locus recorded
   as correct-but-costly. An "already canonical, identity path" test does **not** cover this: the identity path
   returns before any complement, extension or correction machinery runs. So the transport regression lives in
   `tests/test_canonical_hyperbolic_transport.py`, marked `slow` (registered in `pyproject.toml`, with no `addopts`,
   so nothing is hidden by default and a rapid iteration deselects it with `-m "not slow"`). Measured: 55.2 s for
   the transport test, 0.01 s for the derived-minima test, peak tree RSS 3.47 GiB.
2. **The declared order is what makes the minimum negative.** `-v` is always in the orbit of `v`, so the least
   member is frequently a negative tuple; my first expectations were written against the positive representative
   and were wrong. The tests now derive the expected minima from the orbit description (`x = ±5` forced, `y mod 5`
   the invariant on the square-50 slice of `<2> + 0`) instead of guessing coordinates.

## Cost

`tests/test_marked_orbit_canonical.py`: 9 tests, 75 s, peak tree RSS 3.75 GiB (Julia baseline 3.15 GiB, python
590 MiB). The cost is the algorithm itself — every candidate is compared to each emitted representative through
`pointed_isometry` — not test overhead. `tests/test_degenerate_pointed_markings.py` still passes 39/39, unchanged.

## Not claimed

The general #62 iterator for an arbitrary effectively-enumerable integral `Γ`-set: `pointed_isometry` decides
`O(L)`-transport of vectors, and that is the only equivalence wired in. The infinite quotient demonstrated here is
square 0 of `<2> + 0`, one class per `|c|`, which the enumeration produces as a prefix-stable stream with no
completion signal; whether some other family of the port's lattices also gives an infinite quotient is not
investigated.
