# Complexity theory and test budgets

Each file in this directory states the cost of one ported algorithm: the pinned
source routine and its translation, the branches and the mathematical conditions
that select them, the work of each branch with the theorem or bound behind each
count, worst cases, and the unit costs of the backend primitives.

`ci/test_budget_models/` applies this theory to every collected test node, and
`ci/build_test_budgets.py` writes `ci/test_budgets.json`. The pytest plugin
`ci/test_budgets.py` sets each node's pytest-timeout from that table, refuses a
collected node without an entry, and logs one `BUDGET` line per test.
`ci/plan_acceptance.py` packs Acceptance jobs by the same budgets.

## Rules for a budget

- A budget states how long a correct implementation **should** take. It is
  derived by inspecting the mathematics and the code path. It is never read off
  a measured runtime: a measured runtime absorbs the slowness a budget exists to
  expose.
- Every parametrized case has its own budget. Two lattices of one family can
  differ in cost by orders of magnitude, so no per-function maximum and no
  sibling's value.
- Cost is not a function of input size. Classify the case into the branch the
  code takes, then price that branch from the invariants that drive it:
  reflective or not, isotropic or anisotropic and the Witt index, the
  discriminant group A_L and |O(q_L)| (a small-rank lattice with a large
  discriminant group can be the most expensive case), signature conditions,
  cusp, orbit, coset and cell counts, class numbers.
- Compute an invariant of the input exactly when it can be computed without
  running the computation under test; otherwise cite the literature, or state
  the assumption by name.
- No multiplicative safety factor. It scales the slack with the estimate and
  hides regressions in long tests. The budget is the derived maximum itself.
- No minimum and no default.
- Unit costs are what the primitive costs at the size used, in one shared Sage
  session. Anchor them to the reference implementation's own CI bounds where
  those exist: polyhedral_common's `ci_20_reflective.yml` runs all 8821
  ListReflect lattices within 120 minutes. A suite total in the thousands of
  hours means a model is wrong, not that the suite is large.
- Budgets catch order-of-magnitude regressions. Differences of seconds are noise.
- A budget above 20 minutes for a test or witness this repository writes needs a
  smaller input reaching the same branch, or an optimization issue under #64.
- Work shared by several cases (the same lattice's classes, generators or
  permutation model) is charged once per session, not to every case.
- No port factor and no "interpreted Python is slow" pricing. Sage's arithmetic
  is Cython over FLINT, PARI and GAP, and a correct implementation runs its hot
  loops compiled or vectorized. Price every step at the compiled cost of the
  same work the reference implementation does; a budget above the reference
  implementation's bound means the model is wrong. Per-element interpreted work
  in a hot loop is a defect the budget exposes, not a cost it includes.
- Every unit cost is a measured primitive from `ci/test_budget_models/unit_costs.json`
  (written by `unit_costs.py`, which times primitives and never a test), or an
  operation count from the theory times measured primitives. A unit cost typed in
  by hand is a guess and is not allowed. A primitive the models need and the
  benchmark lacks is added to the benchmark first.

## Classifying a step before pricing it

Every step is priced as a function of invariants of its input. A time stated
without that function is not a budget.

| Step | Class | Price |
| --- | --- | --- |
| Constructing any object (lattice, sublattice, module, morphism, homset, finite set, group, category placement) | negligible | Construction is lazy: it stores the defining datum. Derived data is computed on first request and cached. |
| Linear algebra over a field: rank, kernel, image, solving, products, determinants | negligible | Polynomial in the rank at compiled cost |
| Checking that a given morphism is an isometry, an automorphism or an equivariant map | negligible | Linear algebra: compare b(f x, f y) with b(x, y) on generators |
| Rational eigenspaces, invariant subspaces, fixed subspaces, orthogonal complements over a field | negligible | Linear algebra over the field |
| The same constructions over ZZ (saturation, integral kernels, Smith and Hermite forms, sublattice indices) | depends on the algorithm | Module theory over ZZ: price the named algorithm, polynomial in rank and log of the entries |
| Finding an isometry or an automorphism group; any search or enumeration (short vectors, orbits, cosets, perfect forms, roots) | hard | Price from the algorithm and the counts that bound its search: for example the number of vectors of a fixed norm in a bounded region, the size of the automorphism group, the class number of the genus, the orbit or coset count, the perfect-form count |
| A procedure not known to terminate on its inputs (Vinberg on a non-reflective lattice, an unbounded walk) | hard | Price only with a proved termination bound for the case; otherwise state that the case has none |
| A routine taken from the pinned upstream source | read it | Do the complexity analysis of that routine: its loops, their bounds, and the invariants that set the bounds |

The invariants that set prices include the rank, signature and definiteness, the determinant and its factorization, the discriminant group A_L and the size of O(q_L), the primes of bad reduction and ramification, the size of the genus and which theorems apply (strong approximation, Eichler's criterion, Kneser's results on spinor genera), the Witt index, lengths of words and flags, and counts of vectors of fixed norm in bounded regions. A profile dominated by a negligible step is a defect in the code.
