# Issue #53 pointwise stabilizer: what the port supplies, and what it refuses to claim

Ownership, as the issues state it: **research** owns the generic quotient/canonicalizer framework;
**this port** supplies the actual indefinite `Γ`-actions, exact equivalence/stabilizers where they exist, and the
concrete enumeration inputs. This note records exactly that boundary for the `L = U + U` specimen, so nothing here
reads as a claim on research's side.

## The specimen

`L = U + U` with basis `(e1, f1, e2, f2)`, `b(e_i, f_i) = 1`, orthogonal across the two copies.
`Γ = Stab_{O(L)}(e1)` **pointwise** -- not `O(L)`, not a finite-index subgroup.
`v_m = e2 + m·f1`, primitive isotropic, with `b(e1, v_m) = m`.

Measured on the specimen, not assumed: `q(v_m) = 0`, `content(v_m) = 1` and `b(e1, v_m) = m` for `m = 0..4`.

## What is supplied

1. **The action input, verified.** `EichlerRecursion.pointwise_stabilizer_generators(e1)` returns the retained
   `Aut((L, e1))` generators -- 6 of them, computed in 2.5 s -- and re-applies each one to `e1`, raising if any
   fails to fix it pointwise. That check is the difference between "a subgroup fixing the vector" and "some
   isometries": a quotient driven by these is a *pointwise* quotient only if every generator fixes the marking.

2. **A proved separating invariant.** `pointwise_stabilizer_invariant(marked, v)` is
   `(q(v), content(v), b(u_1, v), …, b(u_k, v))`: `q` and content are isometry-invariant, and `b(u_i, ·)` is
   preserved because `g` fixes `u_i` and preserves the pairing, so `b(u_i, gv) = b(g u_i, g v) = b(u_i, v)`.
   `pointwise_stabilizer_separates` therefore returns `True` only when it has **proved** two vectors are in
   different orbits. On the specimen it separates the whole family: `b(e1, v_m) = m`, so the `v_m` lie in
   infinitely many distinct `Γ`-orbits -- the #53 statement, now executable rather than quoted.

3. **Witnessed positive transport, by certificate.** `verify_pointwise_stabilizer_transport(marked, v, w,
   candidate)` checks that a supplied element fixes every marked vector and carries `v` onto `w`. A real instance
   is used in the tests: the isometry swapping `e2` and `f2` (fixing `e1`, `f1`) lies in `Γ` and carries
   `v_m ↦ f2 + m·f1`, a **different** vector with the same invariant, so the positive direction is witnessed by a
   non-identity element rather than by the identity.

4. **Based vector versus unoriented line, measured.** `line_stabilizer_generators(e1)` returns 7 generators and
   **does** contain an element acting as `-1` on `⟨e1⟩`; `pointwise_stabilizer_generators(e1)` contains none. The
   two quotients `L[0]/Stab(e1)` and `L[0]/Stab(Z e1)` are different problems, and a test asserts the
   distinguishing element exists on one side and not the other.

## What is refused, deliberately

* **No `Γ`-equivalence decision, and no `Γ`-orbit G-set.** A G-set needs an exact orbit relation; the only exact
  relation available here is the invariant, which is not complete. Supplying the invariant as an orbit relation
  would merge genuinely distinct orbits, and supplying "not separated" as equivalence would report unproved ones.
  The API therefore exposes no `pointwise_stabilizer_equivalent`, and a test asserts that absence so a later
  addition cannot quietly reinterpret silence as a decision.
* **No completeness claim for the generators** beyond the source's own statement, and in particular **no inference
  from decidability to finite generation** -- the two are independent, and #62 says so explicitly.
* **No parametrisation of the shell.** The `v_m` family is a test instrument for the invariant and the witness
  checker, not a substitute for `L[0]/Γ`. Fixed `m`, the orbit is not decided here beyond the witnesses actually
  produced.
* **No subgroup membership.** Given generators of `Γ`, "is this element in `Γ`" is membership in a subgroup of
  infinite index; the certificate checker takes a candidate element and verifies it, which is a different (and
  decidable) question.

## Numbers

`tests/test_pointwise_stabilizer_shell.py`: 7 tests, 7.2 s, peak tree RSS 433 MiB.
With the facet and canonical files: 36 passed in 23.4 s, 646 MiB.

## Shared-gap routing

Nothing here needs a new shared API. The one thing the port cannot express with maintained operations is
`Γ`-orbit *equivalence* for an infinite-index `Γ` given by generators -- deciding `∃γ ∈ Γ : γv = w`. That is exactly
research's generic quotient/canonicalizer territory under #62, and the port's inputs to it are now concrete: a
verified generator set, a proved separating invariant, and a certificate checker for witnesses.
