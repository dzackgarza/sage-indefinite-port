# Native full integral definite-isometry witnesses

Provider computation: `58709aa58f6075219a192cab5d74b188d0f567f3`.
Published provider with the odd-lattice and shared-genus assertions: `7960e1ae8a7cca25abf689f37ab4745376cfeecf`.

`LatticeIsometryMor._witness_matrix` now uses Sage `QuadraticForm.is_globally_equivalent_to(return_matrix=True)` for full integral definite endpoints. Sage calls PARI `qfisom` and returns M with self(M)=other. Making both forms positive with the same sign, using Hessian 2G and querying the target form against the source form supplies a column map from source to target. The provider verifies unimodularity and `M^T G_target M=G_source` before raising the owned isometry.

The owner assertion in `tests/lattices/test_definite_isometry_witnesses.sage` changes an integral basis by replacing its first vector with itself plus twice the second. This transition is unimodular. Each returned witness has the exact endpoints, preserves every basis pairing and satisfies both inverse laws. The negative binary pair has determinant 23 and the same genus, but its independently proved minima are different, so neither isometry set can contain a witness.

| Owner assertion | Mathematical result | Call seconds | Existing owner gates |
| --- | --- | ---: | --- |
| Positive A2 basis change | Exact integral witness | 8.24 | pass |
| Negative D4 basis change | Exact integral witness | 3.41 | pass |
| Negative E8 basis change | Exact integral witness | 13.87 | pass |
| Negative standard rank-eight lattice basis change | Exact integral witness | 8.59 | pass |
| Same-genus binary pair with minima two and four | Both isometry sets empty | 7.14 | pass |

The initial selected root-lattice/binary session exits 0. The additional odd-lattice session and strengthened shared-genus session each exit 0. The former A2 path reaches the unchanged 60-second hang bound in the OSCAR route; its mathematical assertions are unreached. Its remaining Julia shutdown process was isolated and terminated after the test result, and is gone.

These assertions establish the definite-witness operation on the stated independent inputs. The original I_2_10 plane classification remains required under its unchanged 2400-second bound, with exactly E8(-1) and I_{0,8} quotient types. The full frozen suite still requires one final published port/provider pair.
