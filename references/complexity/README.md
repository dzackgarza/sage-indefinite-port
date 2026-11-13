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
