# Original Allcock primitive-plane assertion

Port: `2e36cff9b11acdc50f7576d5846c45a9dc02c7d9`.
Provider: `7960e1ae8a7cca25abf689f37ab4745376cfeecf`.

The unchanged frozen node `tests/acceptance/test_isotropic_orbits.sage::test_i_2_10_isotropic_classes_are_one_odd_and_one_even[rank1-isotropic_plane_orbits]` passes under its original 2400-second bound in the shared Sage container. The generated parameter label names the second parameter instance: its mathematical isotropic rank is two. The retained JUnit is `.tmp/acceptance-2e36cff-7960e1ae-allcock-planes.xml`.

The assertion computes all returned primitive-plane classes, takes each isotropic reduction and compares the complete parity multiset with the frozen oracle. Thus it establishes exactly two returned classes, one even and one odd. It does not separately assert an isometry of each returned reduction to a named lattice. Pinned Allcock arXiv/9905166/main.tex:232–241, Corollary 4, identifies the two required quotient types as E8(-1) and I_{0,8}.

| Executed evidence | Value |
| --- | ---: |
| Frozen test cases | 1 |
| Failures/errors/skips | 0 |
| Mathematical call seconds | 2007.555 |
| JUnit session seconds | 2016.811 |
| Original bound seconds | 2400 |

The pytest session exits successfully and its process group is empty. Retire local repetition at this pair; final closure still requires this unchanged assertion in the full frozen suite at the final published port/provider pair. O_plus_N, A7, Dawes and the broader source/public/arithmetic scope remain open.
