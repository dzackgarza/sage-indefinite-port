# Contributing to the Sage indefinite lattice port

## Mathematical authority is not held by the implementer

This is a port of pinned upstream mathematical algorithms. The upstream
algorithms, their stated domains and guarantees, and their existing tests are
the starting specification. An implementation contributor's job is to
translate that specification faithfully, not to invent the mathematical
statements against which their own implementation will be judged.

**The contributor implementing a mathematical capability must not introduce,
change, approve, or certify the expected mathematical result of a new test
for that capability.** This applies equally to agents, people, indirect
helpers, generated fixtures, and tests described as “independent oracles.”
Writing a purported proof, finding a citation, or asking another tool to
repeat one's calculation does not make the result independent.

* **Upstream tests:** Translate existing upstream test inputs and assertions
  without altering their mathematical claims. Record the exact pinned upstream
  location and any representation conversion. Where practical, compare actual
  upstream results against the public Sage consumer. No new expected result is
  authorized by the fact that the port happens to return it.
* **Additional mathematical tests:** A separate mathematical test author,
  independent of the implementation contributor and their proposed output,
  must establish the expected result from the pinned upstream computation or
  an independently checked proof **before** the test becomes an acceptance
  assertion. This author must check the actual hypotheses, normalization,
  domain, signs and deduction, not simply certify the implementer's account.
  The source and the independent decision must be inspectable by reviewers.
* **Implementation-side experiments:** The implementer may explore candidate
  inputs, run comparisons, discover counterexamples and use disposable
  diagnostic probes. Those experiments may help debug code. They are **not**
  mathematical acceptance tests, irrespective of how many runs pass or how
  detailed the accompanying justification appears. An implementer may
  nominate a proposed new test for separate consideration, but must not
  install its mathematical expectations as accepted assertions.
* **Unresolved claims:** If no independent mathematical source or proof
  establishes an expected answer, do not invent or guess one. Leave that
  mathematical expectation unaccepted, and investigate the pinned upstream
  behavior. The implementation's output cannot resolve the question.

## What constitutes acceptance

A reviewer must distinguish three questions: whether the port executes, whether
it agrees with the pinned upstream algorithm on the tested inputs, and whether
its general mathematical claims follow from the upstream contract. A test pass
answers neither of the last two questions by itself. Group generation needs
completeness, orbit enumeration needs exhaustion and inequivalence, and
transporters need the correct integral and subgroup constraints. An independent
review must address the relevant assertion, rather than merely endorse a
citation, provenance field, or executable test.

Tests introduced by an implementation contributor before this separation was
established are **not** retrospectively accepted as independent mathematical
evidence. Their results may be informative diagnostics, but cannot justify a
correctness or completion claim without independent mathematical scrutiny.

This is a **contributor responsibility and review rule**, not a claim that a
mechanical checker can decide mathematical truth. No oracle metadata format,
self-signed declaration, CI badge, or automated agent review substitutes for
independent mathematical judgment. Do not build a new verification harness to
replace that judgment.
