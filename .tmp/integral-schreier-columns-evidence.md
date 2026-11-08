# Full integral Schreier actions in a common frame

Port: `2deaefb8cb80e2cdaa85b91c84a70feeebcff923`.
Provider: `2dfd20b685cdb57361f280553345f207864bd75f`.

The unchanged original O_plus_N assertion at 2e36cff/7960e1ae is terminal under its original 2400-second bound while constructing a full Schreier composite as an owned formed-module isometry. The provider already supplies the exact orbit, generator permutations, letters and transversals.

The port now lowers each supplied factor through its actual module morphism and the specified integral-lattice inclusion. For a map f, endpoint comparisons give the common-frame column action C_target^(-1) M_f C_source. Every Schreier action is the exact rational product T_(s.p)^(-1) S T_p; admission over ZZ occurs after this complete product. The existing Nielsen reduction consumes every full integral action and raises its surviving actions in the lattice isometry homset. The one-point orbit retains its direct-letter branch.

Pinned GroupAction.g:869–886 constructs the integral stabilizer before simplifying its generators. MatrixGroupSimplification.h:650–751,1386–1487 supplies the product/inverse Nielsen route. The finite image is used only for the Schreier indices; reduction retains the full integral maps, including the kernel of that image.

| Selected mathematical assertions | Result |
| --- | --- |
| Every forward/inverse Schreier action in 2U and the rational-factor hyperbolic example | pass |
| Noncommuting axis action on a nontrivial integral frame, including its central sign | pass |
| Integral generating families, primitive unipotent powers, symmetric action and large exact coefficients | pass |
| Public pointwise maximal stabilizer in 4U equals the full skew-parameter module | pass |
| Selected assertions | 11 pass |
| Selected session seconds | 40.80 |

The noncommuting example retains the six-point lattice orbit and compares every source family action in the frame with axis indices 1,2,3. Both directions are checked on the full integral basis. The generating-family assertions retain unipotent and sign actions invisible in finite quotients, including coefficients beyond machine integer width.

The original unchanged O_plus_N line/plane assertion is active at this published pair under its original 2400-second bound in `sage-indefinite-proof-2deaefb`. Preserve it until terminal. Optional log/JUnit: `.tmp/o-plus-N-2deaefb-2dfd20b6.{log,xml}`. The small assertions prove the stated factor-realization and generating-family operations; the original classification remains open until its own result is available.
