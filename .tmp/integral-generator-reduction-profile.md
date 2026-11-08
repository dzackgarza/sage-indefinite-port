# Full integral-generator reduction

The unchanged I_1_1 assertion passes at port `2cdf004966bb1d16c79f0a19ddc79a0c0f7e63ae` / provider `bc6bd57262e576194e3c8764c83f3bb1d081ba52`, within its original 2400-second limit. JUnit case time is 2272.321 seconds and the required discriminant stabilizer index is 1. The exact report is in `.tmp/acceptance-2cdf004-evidence.md`; the interval sample below does not independently establish that assertion.

A ten-second main-thread sample reaches `_integral_generator_set` while it consumes `_integral_stabilizer_family`. Inclusive times overlap and must not be added.

| Sampled operation | Seconds |
| --- | ---: |
| Next integral generator at recursive.py:158 | 3.60 |
| Full Schreier-factor composition at recursive.py:106 | 3.55 |
| First product/inverse candidate at recursive.py:202 | 1.60 |
| Third candidate at recursive.py:204 | 1.40 |
| Second candidate at recursive.py:203 | 1.30 |
| Fourth candidate at recursive.py:205 | 1.10 |
| Direct L1 norm in normalized at recursive.py:143 | 1.30 |
| Inverse L1 norm at recursive.py:144 | 0.70 |

A separate live stack in the first pair product has direct-action norms 71 and 14, with the last supplied action of norm 342. These sampled actions therefore fit the current certified int64 branch. This input does not support changing the large-integer arithmetic route. The completed original I_1_1 assertion passes; this sample supplies no arithmetic failure.

The pinned source was read at [MatrixGroupSimplification.h](https://github.com/MathieuDutSik/polyhedral_common/blob/1592b2463bc7a04d94e900bf8efbada607405d8d/src_group/MatrixGroupSimplification.h): `f_reduce` at lines 685–750 computes the same four product/inverse pairs; the exhaustive V2 scheduler at lines 798–1074 retains untreated pair intervals and stops only after no pair improves. Lines 1454–1476 select the entrywise L1 complexity, and lines 1493–1605 bound fixed-width products. The existing port worklist retains complete integral actions and exact product/inverse witnesses. No new generator reduction or arithmetic change follows from this interval sample.

Optional local samples are `.tmp/enriques-2cdf004-bc6bd57-reduction.json` and `.tmp/enriques-2cdf004-bc6bd57-reduction-locals.txt`. Evidence for the full assertion at the latest port/provider pair remains required.


The unchanged O_plus_N line-and-plane assertion at a170e2e/0580f98d is terminal at its original 2400-second bound. The line assertion completed. The final plane stack reaches `subspace_extensions -> shortened_extension -> integral_flag_generators -> _integral_generator_set`, in `normalized(a_inverse @ b, b_inverse @ a)`. That proof includes native definite groups and direct selected-generator images, but predates the coefficient-family fractional-span lift and retained coset transitions. The unchanged proof is active at 66576eb/16c798db; preserve it until terminal. If it still fails in pair reduction, obtain its actual full integral action family and compare against the pinned exhaustive scheduler before changing the reduction. Finite-image duplicates cannot justify discarding integral kernel generators.


A one-shot snapshot of the active I_1_1 proof at 583b36d/16c798db supplies full integral basis actions seen at that execution point. It is a prefix of the lazy Schreier family, with raw-key insertion order unavailable; future family entries are not captured. It cannot establish the whole original stabilizer or the original input order. The retained normalized actions retain their full inverses. Snapshot creation takes 0.0169 seconds through Python's documented remote debugger interface, included in the proof's original bound. The optional local input is `.tmp/integral-state-2942513.pickle`. Use the exact matrix data to compare the pinned reduction before making another reduction change; finite-image duplicates do not certify redundancy of these full integral maps.
