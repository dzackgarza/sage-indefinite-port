# Contribution boundaries for the Sage indefinite port

This repository **does not own a separate mathematical category framework**.
Its public API uses the formed-lattice objects, functors, morphisms and
tensor operations provided by [research](https://github.com/dzackgarza/research).
A reference to an upstream project's contribution policies is not an
instruction to copy its construction programme into this port.

The binding contract is [AGENTS.md](../AGENTS.md),
[port-design.md](port-design.md#implementation-assignments), the
[pinned source map](source-to-sage-translation-map.md), and the
[issue tree rooted at #4](https://github.com/dzackgarza/sage-indefinite-port/issues/4).
When a *research-owned* implementation genuinely needs repair, its governing
policies are [research/CONTRIBUTING.md](https://github.com/dzackgarza/research/blob/main/CONTRIBUTING.md).
The generic framework's separate architectural constraints live in
[sage-categories/CONTRIBUTING.md](https://github.com/dzackgarza/sage-categories/blob/main/CONTRIBUTING.md).
Those documents govern their respective owners, not a new port-local
reimplementation.

The port-specific consequences are:

- **Mathematical ownership.** Represent lattice isometries, marked objects,
  fixed pointwise markings, sublattices, flags, forms, discriminant modules
  and group actions using the existing preamble category/morphism objects.
  A convenience matrix, coordinate tuple or selected backend representation
  is never substituted for the public mathematical value.
- **Exact backend delegation.** Use existing exact FLINT, Sage, GAP,
  Normaliz/PPL or other appropriate arithmetic routines for assigned
  generic callees; implement the *specialized indefinite caller* according
  to the pinned upstream algorithm. Never replace a missing general
  computation with an unproved source-specific special case, duplicate
  subgroup search, arbitrary search budget or hand-written alternative.
- **Source and error provenance.** Specify the upstream routine, the exact
  failing public consumer and any necessary preservation theorem before
  changing its code. Preserve witnesses, hypotheses, boundary cases,
  normalizations, integral structures and termination/certification conditions.
- **Defect ownership.** Send genuinely reusable arithmetic/framework defects
  to their owning `research` interface with failing inputs and an independently
  meaningful contract. Do not edit the research checkout from a port worktree
  or impose an unrelated research redesign on the port.
- **Frozen acceptance.** Changes to oracle datasets or suite assertions
  require issue #33's recorded source-backed justification. Follow its
  accepted branch-representative sampling amendment; never drop cases to
  make a run green.
- **Publication and review.** Local commit/push checks cover syntax only.
  Detailed CI quality/type diagnostics are nonblocking during development;
  mathematical acceptance is the committed-port CI evidence specified in
  TODO.md. Tooling, formatting or type cleanup does not interrupt an
  unfinished port algorithm.

These are the relevant constraints for a remote checkout. Host-specific
absolute paths, local memory vaults, machine-level Docker installations and
project-specific `just` commands in *other repositories* are not required
port preparation.
