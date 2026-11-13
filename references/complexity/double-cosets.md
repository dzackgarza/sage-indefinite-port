# Double cosets of finite permutation groups: complexity

## Algorithm and pinned source

The consumer is `G.double_cosets(H, H)` on owned permutation groups
(sage-extensions `categories/group/groups.py:6665-6690`). When `G/H` is a
finite `G`-set, `H \ G / H` is the set of `H`-orbits on `G/H`
(`categories/group/g_sets.py:1032`). The pinned decomposition test is
`TestDoubleCosetPermDecomposition`
(`references/vendor/polyhedral_common@1592b246/CI_tests/01_RatIntAutomorphy/ProcessExamples.g:179-253`).

## Branches

| Branch | Condition | Route |
| --- | --- | --- |
| Finite right cosets | `G/H` realized as a finite `G`-set | `H`-orbits on `G/H` |
| Finite left cosets | `G/Q` finite, `G/H` not | inverted `H`-orbits on `G/Q` |
| None | neither finite | assertion |

All acceptance cases are finite permutation groups of degree 55 to 96, so
they take the first branch.

## Cost

- Group orders: Schreier–Sims for a group of degree `d` with `g`
  generators costs `O(d^2 g + d log^3 |G|)` (Seress, *Permutation Group
  Algorithms*, Thm. 4.2.4). For `d <= 96` this is milliseconds.
- Coset space `G/H`: `[G:H]` points. The `H`-orbits cost
  `[G:H] * g_H` coset images; each image is one membership test in `H`
  through its stabilizer chain, `O(d log |H|)`.
- The owned layer wraps each point once: `[G:H]` preamble objects.

Budget: `0.5 + 0.1 (1 + 50 d^2 (g_G + g_H) U_PERM) + [G:H](g_G + g_H)(d U_PERM + 4 U_ELT) + 20 [G:H] U_ELT`.

## Known worst cases

A large index `[G:H]`: the coset space is enumerated explicitly.

## Unit costs of backend primitives

The unit costs are those of `orthogonal-group-generators.md`, "Unit costs of
backend primitives".
