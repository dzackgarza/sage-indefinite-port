# Mathematical specification, implementation, and acceptance

## Governing principle

This project ports pinned upstream mathematical algorithms to Sage. **The
authority to establish what is mathematically correct is separate from the
authority to implement it.** No contributor may establish the acceptance
standard for their own mathematical implementation. This is a division of
responsibility, not a prohibition list or a declaration that paperwork
certifies mathematical truth.

### Jurisdiction

| Authority | Responsibility | Source of authority |
| --- | --- | --- |
| **Mathematical specification** | Defines the mathematical operations, domains, outputs, completeness, and equivalence criteria to be preserved. | Pinned upstream code and tests, their mathematical sources, and separately established project requirements. |
| **Independent test author** | Writes **every test**, its fixtures, and its expected results, including faithful upstream translations and nonmathematical tests. | Authoritative upstream code/tests, external mathematics, or the pre-existing behavioral specification; never the implementation's observed result as authority. |
| **Implementation contributor** | Translates the specified algorithms and interfaces, debugs failures, and demonstrates behavior against separately authored tests. | The specification and independently accepted tests. May propose test questions and run temporary probes, but **cannot write or amend committed test code or expected outcomes for the work they implement**. |
| **Independent reviewer** | Challenges the correctness of **every test** before acceptance and periodically afterward, with particular scrutiny for mathematically false assertions. | Direct inspection of the cited source, independent calculation or proof, and the executable assertion and fixtures. Must be distinct from the implementer and test author. |
| **Project owner** | Assigns test authors and reviewers independently of the implementation contributor, resolves disputes, and decides acceptance. | The project's stated requirements; neither passing CI nor an agent's own signoff substitutes for the owner's decision. |

An implementer cannot create independence by launching another agent, asking
another model to endorse its result, rewriting a proof under another name, or
selecting its own reviewer. A reviewer must actually exercise independent
mathematical judgment. Distinct identities alone are not a review.

## Requirements applying to **every test**

The following are universal, not merely rules for new mathematical tests. They
cover every committed test, upstream translation, regression, parameterized
case, fixture assertion, and engineering test. **No pre-existing test is
grandfathered into mathematical authority.**

1. **Cite the governing source.** Each test must identify an inspectable,
   specific and versioned authoritative source for its *expected behavior*:
   the pinned upstream test or algorithm and its precise location, the
   independent mathematical proof and its exact claim, or (for nonmathematical
   checks) an independently established protocol/interface requirement.
   The citation must cover the test's actual inputs, assumptions, and asserted
   result. Merely naming a paper, a related theorem, the implementation under
   test, or the test author's own output does not satisfy this requirement.
   Where no source establishes the assertion, it is not an accepted test.
2. **Independently author it.** Every accepted test and its mathematical
   expected result must be written by a test author independent of the
   contributor implementing the behavior tested. This includes direct
   translations of upstream tests: source authority comes from upstream, but
   the implementer still must not be the sole author of the translation or
   fixtures against which they are judged. The implementer can suggest
   questions or run temporary probes, never promote their probes as tests.
3. **Independently review and periodically re-review it.** A reviewer distinct
   from both the test author and implementer must inspect every test before it
   has acceptance authority. **Every test is reviewed again at each release
   and at least once every 90 calendar days**, whichever comes first. A
   changed assertion or cited source immediately invalidates the previous
   approval. The project owner selects the reviewers, not the implementer.
   Periodic reviews should use a different independent reviewer from the
   previous review wherever feasible.

**The object of review is to detect false tests.** For mathematical tests the
reviewer must actively attempt to falsify the assertion: check that the
claimed upstream computation or theorem genuinely applies, reconstruct the
independent calculation or proof, and examine counterexamples, representation
changes, hidden fixture assumptions, signs, integrality, negative decisions
and claims of exhaustiveness or generation. An executable test can pass while
its expected mathematics is false. For nonmathematical tests, review the
assertion against the cited contract and implementation-independent fixtures.
Reviewing only metadata, approving a fellow agent's explanation, or observing
green CI is not review.

An accepted review must identify the exact test(s), source/revision, author,
reviewer, mathematical reasoning or independently reproduced result, review
date, and acceptance decision. A single genuinely general argument can cover
a parameterized family, but a sample is not review of the remaining tests.
Reviews can be recorded in the existing review/issue history; this rule does
not require a new tracking service or form. A test lacking current independent
review may be run as a diagnostic but **must not support a mathematical
correctness or completion claim**.

## What tests are authorized

### 1. Translated upstream tests

Upstream tests are the default mathematical regression corpus. Their inputs,
assertions, hypotheses, and intended mathematical meaning must be preserved.
A translation may change representation, indexing, or syntax only insofar as
the source and target mathematical objects are identified correctly. Its
origin must be traceable to the pinned upstream revision, file, and test or
algorithm. A change that adds an expected mathematical result, strengthens or
weakens an assertion, changes a distinguishing input, or discards a case is
**not** merely a translation: it invokes the new-test authority below.

An independent test author translates upstream tests; the implementer may run
them but does not author the translation. Independent source-fidelity review
is required before they count as evidence of upstream equivalence. Where
upstream execution is possible, actual upstream outputs are preferred over an
agent's restatement of what upstream supposedly does.

### 2. Additional mathematical tests

New mathematical assertions may enter the accepted test corpus **only after
independent mathematical specification and independent mathematical review**.
The test author's work must establish:

* **The question being tested:** the specified upstream or separately
  authorized mathematical contract, exact input and hypotheses, and which
  required behavior the test distinguishes.
* **The correct answer:** the exact pinned upstream computation (including
  input/output identification), or a complete and checkable independent
  mathematical argument. Inferences, signs, normalization, integrality,
  subgroup restrictions, and negative or exhaustive claims must be justified
  at their actual scope.
* **The acceptance assertion:** the precise proposition the executable test
  checks, with no unstated inference from a sample result to completeness of a
  general algorithm.

The implementation's output, an expected value invented by its author, a
plausible citation, an autogenerated proof, and agreement between independently
*run* copies of the port are not independent mathematical specifications.
If the pinned upstream cannot answer a proposed question and no independent
proof establishes it, the expectation is **undetermined**. It does not become
an acceptance test by running it repeatedly.

The implementing contributor may identify an omitted case or propose a test
*question*. Another, independently appointed mathematical test author must
establish its expected answer without treating the implementer's output or
explanation as an oracle. The implementing contributor does not author, alter,
or approve that mathematical assertion.

### 3. Implementation diagnostics

Implementers may write and run disposable experiments, compare
representations, inspect candidate answers, and investigate counterexamples.
Independently authored nonmathematical engineering tests may also check
transport, parsing, error handling, and similar behavior. **These do not acquire mathematical
acceptance status merely by being committed, passing, or acquiring metadata.**
Any diagnostic proposed for promotion into a committed test must be independently
written and reviewed under the universal test requirements above.
When classification is disputed, mathematical acceptance authority decides.

## Mandatory mathematical review

**Admission precedes acceptance and is never permanent.** An independent test
author submits the precise assertion with its independent source/result or
proof. A reviewer who neither implemented the feature nor authored that assertion
then independently checks the source, input correspondence, hypotheses,
calculation or proof, and the exact assertion the test will execute. A reviewer
must reject a source that is merely relevant in subject matter but does not
imply the claimed result. For claims of generation, exhaustion, termination,
uniqueness, integral equivalence or nonexistence, the review must establish
that claim rather than just verify some returned examples.

**Review is an attributable mathematical decision.** The review must be
inspectable and identify the independently specified test or coherent test
family, its immutable source/proof reference, the accepted mathematical
assertions and hypotheses, the test author, the *distinct* reviewer, and the
reviewer's actual mathematical reasoning or independently reproduced result.
The same obligations apply to the mandatory periodic review, including an
explicit search for mathematically false expectations, not just a declaration
that the earlier approval still exists.
One mathematical argument can authorize an entire parameterized family when
the argument truly covers the parameters; repeated forms and templates are
not required. Rejected or unresolved assertions remain outside accepted
mathematical evidence.

Changing an accepted test's mathematical meaning (including expected values,
inputs that distinguish cases, skips, tolerances, or completeness criteria)
requires a new independent mathematical decision. Mechanical refactoring
of tests belongs to the independent test author; preserving the assertion
requires an independent source-fidelity check, not a new proof. The
implementer may fix the production code against accepted assertions but
has no authority to amend or redefine the tests.

## Acceptance of the port

Executable tests are witnesses to behavior on their inputs, **not proofs of
the complete mathematical contract**. Acceptance of an algorithm requires
review of the upstream-supported domain and every relevant dispatch branch,
including its completeness, stopping conditions, group-generation claims,
orbits, stabilizers, transporters, and public-consumer semantics. Additional
capabilities beyond upstream must be explicitly recognized as *new
mathematics*, independently specified and justified; they cannot silently be
presented as part of a faithful translation.

Previously introduced implementer-authored tests remain useful
for diagnosis but receive **no retrospective independent authority** from
their existence, passing status, citations, or provenance fields. They count
toward acceptance only after their assertions and executable form have been
independently authored or independently re-established, independently reviewed,
and included in the periodic review cycle.

CI, source-provenance scripts, checklists, and `ORACLE-*` fields may
check syntax or record references. **They cannot appoint an independent
mathematical authority, decide a theorem, or approve a test.** This document
governs contributor conduct and mathematical review; it does not claim
automated enforcement of mathematical truth.
