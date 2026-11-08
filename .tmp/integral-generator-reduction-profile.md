# Full integral-generator reduction

The unchanged I_1_1 assertion is active at port `2cdf004966bb1d16c79f0a19ddc79a0c0f7e63ae` / provider `bc6bd57262e576194e3c8764c83f3bb1d081ba52`, with its original 2400-second limit. The sample below does not establish its final index assertion.

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

A separate live stack in the first pair product has direct-action norms 71 and 14, with the last supplied action of norm 342. These sampled actions therefore fit the current certified int64 branch. This input does not support changing the large-integer arithmetic route. The original I_1_1 assertion remains live and must finish before assigning a new failure.

The pinned source was read at [MatrixGroupSimplification.h](https://github.com/MathieuDutSik/polyhedral_common/blob/1592b2463bc7a04d94e900bf8efbada607405d8d/src_group/MatrixGroupSimplification.h): `f_reduce` at lines 685–750 computes the same four product/inverse pairs; the exhaustive V2 scheduler at lines 798–1074 retains untreated pair intervals and stops only after no pair improves. Lines 1454–1476 select the entrywise L1 complexity, and lines 1493–1605 bound fixed-width products. The existing port worklist retains complete integral actions and exact product/inverse witnesses. No new generator reduction or arithmetic change follows from this interval sample.

Optional local samples are `.tmp/enriques-2cdf004-bc6bd57-reduction.json` and `.tmp/enriques-2cdf004-bc6bd57-reduction-locals.txt`. Whole-consumer timing and the final frozen assertion remain required.
