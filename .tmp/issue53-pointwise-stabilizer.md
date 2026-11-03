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

## The decision itself, now implemented

`indefinite/pointwise_shell.py` decides `Gamma`-equivalence **exactly** on the primitive isotropic shell and
returns an explicit element of `Gamma`. The theorem it rests on:

> `Gamma = N ⋊ Theta`, with `N = {u(rho, gam) : rho, gam in Z}` the kernel of `Gamma -> O(Q)`, `Q = e1-perp/<e1>`
> isometric to `U`, and `Theta` mapping isomorphically onto `O(Q) = {I, -I, swap, -swap}`.

The kernel computation is the substance: a `g` fixing `e1` and acting trivially on `Q` has `g(e1) = e1`,
`g(e2) = e2 + rho e1`, `g(f2) = f2 + gam e1`, and writing `g(f1) = lam f1 + z`, the four isometry conditions give
`lam = 1` and `z = -gam rho e1 - gam e2 - rho f2`, i.e. exactly `u(rho, gam)`. `O(U) = {I, -I, swap, -swap}` is a
two-line count from `M [[0,1],[1,0]] M^T = [[0,1],[1,0]]` over `Z`, and `sigma = diag(1,1,-1,-1)`, `tau = (e2 f2)`
realize the two generators inside `Gamma`.

On `v = a e1 + m f1 + c e2 + d f2` the action is `(c, d) -> (c - m gam, d - m rho)` and `a -> a + c rho + d gam -
m gam rho` under `u`, `(c,d) -> (-c,-d)` and `(c,d) -> (d,c)` under `sigma`, `tau`, with `m` fixed. **This was
checked numerically before being used**: `u(1,3)` sends `(-6,1,2,3)` to `(2,1,-1,2)`, exactly the predicted values.

So, on the primitive isotropic shell (`a m + c d = 0`, content one):

* `m != 0`: `a = -cd/m`, and the orbit is the class of `(c,d)` in `(Z/m)^2` modulo `m Z^2` and the four
  symmetries;
* `m = 0`: `c d = 0`, `a` moves by the ideal `gcd(c,d) Z`, and the orbit is the class of `(c,d)` plus
  `a mod gcd(c,d)` -- with `gcd(0,0) = 0`, so `e1` and `-e1` are different orbits, which is the based-vector
  statement made concrete.

`gamma_class` computes the key, `gamma_witness` returns the constructed element (`sigma`/`tau` composed with
`u(rho, gam)`, with `(rho, gam)` solved through Smith normal form), `None` means **proved** inequivalence.

Independence of the evidence, since a classification checked only against itself proves little: every retained
generator of `Aut((L, e1))` is applied to eleven shell vectors and must preserve the class key. That catches a key
that is too *fine* -- the plausible error -- and it passed on all six generators.

## What is refused, deliberately

* **No `Γ`-orbit G-set over the *whole* shell.** `pointwise_shell.gamma_equivalent` is an exact decision on the
  primitive isotropic shell of this specimen; what is still missing for a G-set is the same decision for the
  imprimitive and non-isotropic parts, which the classification does not cover and which `shell_coordinates`
  therefore refuses. The generic enumerator built on a decidable relation is research's framework.
* **No `pointwise_stabilizer_equivalent` on `EichlerRecursion`.** The invariant-based helpers there stay
  non-deciding by design, because a general lattice has no analogue of this specimen's normal form; the exact
  decision lives in `pointwise_shell` where the theorem holds, and a test asserts the absence of a general
  equivalence method so the invariant cannot be misread as one.
* **No completeness claim for the generators** beyond the source's own statement, and in particular **no inference
  from decidability to finite generation** -- the two are independent, and #62 says so explicitly.
* **No parametrisation of the shell.** The `v_m` family is a test instrument for the invariant and the witness
  checker, not a substitute for `L[0]/Γ`. Fixed `m`, the orbit is not decided here beyond the witnesses actually
  produced.
* **No subgroup membership.** Given generators of `Γ`, "is this element in `Γ`" is membership in a subgroup of
  infinite index; the certificate checker takes a candidate element and verifies it, which is a different (and
  decidable) question.

## Numbers

`tests/test_pointwise_stabilizer_shell.py`: 7 tests, 7.2 s, peak 433 MiB.
`tests/test_pointwise_shell_decision.py`: 8 tests, 6.0 s, peak 426 MiB -- the exact decision, its witnesses, the
generator cross-check, and the refusals outside the shell.
With the facet and canonical files: 36 passed in 23.4 s, 646 MiB.

## Shared-gap routing

The exact decision is now supplied for the specimen's primitive isotropic shell, by normal-form classification
rather than by a generic subgroup-membership algorithm: `(Z/m)^2` up to four symmetries for `m != 0`, and
`gcd(c,d)`-ideals for `m = 0`. What remains genuinely open, and what research's framework needs from #54's
arithmetic decision source in general, is the same question for an arbitrary infinite-index `Γ` where no such
normal form is available. The port's inputs there are the verified generator set, the proved separating invariant
and the certificate checker.
