# #63: public consumer routes against actually accepted cases

Acceptance of this unit is: "every public consumer route's claimed availability is checked against the
accepted cases at one recorded port/provider pair, each failing exposed operation is filed at its owner,
and no readiness claim rests on the absence of `@unfinished` markers or xfails."

This note is that check. It deliberately does **not** reason from marker absence.

## The marker state, recorded so nothing leans on it

```
@unfinished decorators in src/        : 0
xfail markers in tests/acceptance/    : 0
```

`#63` says explicitly that this state is not evidence. It is recorded only because the unit's premise
is that the absence exists and tempts a readiness claim; the table below is what decides availability.

## The five declared consumer entry points

`tests/acceptance/consumer.py` names exactly five `require(...)` targets, and all five are defined once
in the source:

| entry point | defined in |
| --- | --- |
| `indefinite.marked_lattices:isometry` | `marked_lattices.py` |
| `indefinite.marked_lattices:automorphism_group_generators` | `marked_lattices.py` |
| `indefinite.marked_lattices:primitive_marking_classes` | `marked_lattices.py` |
| `indefinite.edgewalk:edgewalk_fundamental_domain` | `edgewalk.py` |
| `indefinite.lorentzian_cells:perfect_domain_traversal` | `lorentzian_cells.py` |

## Route → accepting cases → verdict, at recorded pairs

Each row records the run that accepted it: the port commit and provider revision are the ones the gate
recipe wrote into `.tmp/ci/results/<family>-0/{commit,research_commit}` next to that run's `junit.xml`.

| entry point | acceptance module(s) gating on it | accepted run | verdict |
| --- | --- | --- | --- |
| `isometry` | `test_isometry` | port `823b1075` + `cdafca84`: **166 passed** (1279.75 s) | **accepted at a recorded pair** |
| `isometry`, `automorphism_group_generators`, `primitive_marking_classes` | `test_centralizers` | port `73cded56` + `cdafca84`: **5 passed** (78.71 s) | **accepted at a recorded pair** |
| `primitive_marking_classes` (port side) | `test_boundaries` (dawes and dtower cases gate on `MARKED_CLASSES`) | none at a recorded pair | **not accepted at a pair**; the gates exist, the runs at this pair do not |
| `perfect_domain_traversal` | `test_perfect_domains` | port `c746b913` + `cdafca84`: **92 passed** (357.79 s) | **accepted at a recorded pair** |
| `automorphism_group_generators` via the Enriques/boundaries path | `test_boundaries` | none at a recorded pair | **not accepted**: the Enriques family is the known capped case, and `#64` owns it |
| `automorphism_group_generators` via orthogonal groups | `test_orthogonal_group` | none at a recorded pair | **not accepted at a pair**; the A7 node passed once at port `0d1fcb73` (662 s, peak 1616 MiB), which is not this pair |
| all three marked entry points via stabilizers | `test_vector_stabilizers` | none at a recorded pair | **not accepted at a pair**; per-node evidence exists at older pins only |
| `isometry`, `primitive_marking_classes` via vector orbits | `test_vector_orbits` | none at a recorded pair | **not accepted at a pair** |
| `edgewalk_fundamental_domain` | `test_reflective` | none at a recorded pair | **not accepted at a pair**; the family is branch-sampled by `#64` |

Independent of the five: the double-coset decomposition is provider-owned and its 18 cited instances
passed (**19 cases**) at port `d21b69ae` + `cdafca84`; it is not a port route.

## What this establishes, and what it refuses to

Established: three of the five declared routes have an accepting run at a named pair with per-case
JUnit — `perfect_domain_traversal`, the marked isometry/automorphism/classe routes through
`test_isometry` and `test_centralizers`. The source defines all five, so "exposed" is true for all five
while "accepted" is true for three.

Refused: that the two remaining routes (`edgewalk_fundamental_domain`, and the marked routes as
exercised by the boundaries/stabilizer/orbit families) are ready. Their files exist, their cases
collect, and no marker suppresses them — and none of that is evidence. `#64` owns bringing those
families to one pair, and this note does not pre-empt it.

Also refused: any claim that the *general* capability is delivered. `AGENTS.md` and `README.md` already
state that arbitrary-rank isotropic k-planes and partial flags remain required and incomplete
(`#44`), which is a scope statement this reconciliation does not soften.

## Failing exposed operations, filed at their owner

One is now identified, and it is filed.

`test_enriques_modular_groups_match_dutour_sikiric_hulek[case_87]` — which gates on
`automorphism_group_generators` and `primitive_marking_classes` — fails at port `1bedc833` + provider
`cdafca84` in ~140 s against a 1860 s cap, not by capping. `group.group_generators()` for the
intersection subgroup `{g in O(N) : g lies in every selected subgroup}` traces through the provider's
`predicate_subgroups.py` preimage/coset construction and dies with
`OverflowError: value too large to convert to int` inside a `_TransversalCoset` repr, after a
superseded `cached_method` `KeyError`.

It is routed to `sage-extensions` as `.tmp/handoff-owner-enriques-coset-overflow.md`, with the raw
`--tb=long` log beside it, and with the contract the consumer needs stated: generators or a typed
refusal, never an untyped `OverflowError` (`#33`'s typed-error invariant).

The rest of the reconciliation found no exposed operation failing *as a port route* at a recorded
pair: the four families run at a pair all pass. The remaining items are coverage gaps, which belong to
`#64`, not defects — reporting an untested route as a broken one would misdescribe it.

Attribution note: the failure was proven independent of the port's `#67` caching change by running the
same case against a genuine pre-hoist copy — identical failure at 136.76 s against 141.62 s.
