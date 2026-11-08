# Stage-D mechanism comparison artifact (U+U(2)): distinct rational generators, their actions, and their conversions

Provider: sage-extensions `c2c81a3` (+ record `c62c79c`): `Tensor._tensor_hash` no longer raises over
unhashable owned scalars. Changed-source provenance for the specimen record: provider
`src/sage_extensions/tensors/tensor.py` sha256 `8da116e10f939e86eb3f3ef849d8f541b3c1f763f254042efd703543eaa1d029`
at `c2c81a3`, imported at runtime as `/repos/sage-extensions/src/sage_extensions/tensors/tensor.py`
(proot binding of that checkout, PYTHONPATH `/repos/sage-extensions/src`). The attribution probe's
`lattice_morphisms.py` digest does not name this repair; the tensor.py digest does. Port: `d5507b07`
(this checkout; its stage-D witness-scope correction, det +1 / proper subgroup / index unproved, is a
different specimen and stays evidence). Guard: PEAKRSS_LIMIT_MIB=6144, single heavy slot, launchers
`.scratch/launch-sign-loss{,2,3}.sh`.

- `probe-generator-action-comparison.py` (+ its log): specimens U (control), U(2) (control), U+U(2).
  Prints at every layer: each of the 7 distinct supplied rational-group generators' ambient tensor and
  its image of the line vector e; each stabilizer generator's ambient tensor, its integral restriction
  rows and its image of e (the converted output); the SIGN-CHECK membership of the -e element with its
  exact integral restriction; INDUCTION's per-generator line actions; form preservation of every
  integral generator. After the repair: 28 stabilizer generators with distinct honest tensors, among
  them diag(-1,-1,1,1) ambient and integral, plus integral shears; INDUCTION line signs include -1 on
  U+U(2); controls U/U(2) unchanged (1 generator, sign -1, form preserved).
- `probe-letter-layers.log`: the letter layer after the repair - `_signed_generator_image` and the
  augmentation evaluate every one of the 7 free letters to its own generator (pre-repair: s1..s6 all
  read s0's value). Equal-map consistency of separately built equal maps, with equality plus shared
  hash, is pinned provider-side in `tests/lattices/test_isometry_hash_consistency.sage`
  (`test_rational_scalar_isometry_hash_is_total_and_faithful`,
  `test_subgroup_augmentation_letters_evaluate_to_their_own_images`).

The original N = U+U(2)+E8(-2) stage-D measurement with this repair runs guarded under job jmw0qw891;
it is not duplicated here. Acceptance (line/plane/flag 5/10/19, |image| = 92897280) remains pending on
that run and the original call's verdict.
