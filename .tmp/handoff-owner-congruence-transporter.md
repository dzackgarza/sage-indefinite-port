# Resolved record: the congruence transporter is GAP's orbit / action-homomorphism / word chain

**Status: not an owner gap.** An earlier version of this file claimed the transporter in the image of `GL_k(Z)`
"is not supplied by any maintained contract". That was wrong, and the wrong conclusion came from three calls made
with ill-typed arguments, not from missing capability. Recorded here so the mistake is not repeated.

## What the failed calls actually proved

1. `H.orbit(...)` on a `FinitelyGeneratedMatrixGroup_gap_with_category` — the attribute does not exist. That says
   nothing about the group's action: GAP's `Orbit` is reached through `H.gap()`.
2. `libgap.Orbit(G, [1,0], libgap.OnRight)` — `GAPError: Action not well-defined`, and
   `libgap.RepresentativeAction(...)` returned `fail`. Cause: the group's matrices are over `GF(5)` while the point
   was a list of **integers**. GAP's `Z(5)` is a primitive root (`Z(5)^2` is the field element 4), so even a
   hand-built point has to be a field element, not a small integer.
3. `PreImagesRepresentative(EpimorphismFromFreeGroup(H), rep)` — `no method found ... The 2nd argument is 'fail'`.
   The word map was fine; its second argument was the `fail` produced by (2).

## The contract that works, with typed arguments

    ring  = libgap.ZmodnZ(modulus)          # GF(p) when the modulus is prime
    one   = libgap.One(ring)
    point = libgap([libgap(entry % modulus) * one for entry in row])     # a row of the ring's own module
    orbit = libgap.Orbit(group.gap(), point, libgap.OnRight)             # rows act on the right here
    action = libgap.ActionHomomorphism(group.gap(), orbit, libgap.OnRight)
    permutation = libgap.RepresentativeAction(libgap.Image(action), source_position, target_position,
                                              libgap.OnPoints)
    element = libgap.PreImagesRepresentative(action, permutation)         # back in the matrix group
    word    = libgap.PreImagesRepresentative(libgap.EpimorphismFromFreeGroup(group), element)
    syllables = [int(entry) for entry in libgap.ExtRepOfObj(word)]        # [generator, exponent, ...]

Minimal trace, `k = 2`, modulus 5, `(1,0) -> (2,0)`: `orbit` has 24 elements, `Image(action)` has order 240,
`permutation = (1,20,3,17)(2,8,7,10)...`, `element = [[Z(5), 0], [Z(5)^2, Z(5)^3]]`,
`word = x1^-1*x2^-1*x1^-2`, syllables `[1,-1,2,-1,1,-2]`.

**The lift-word obligation** stays in the port: the word is composed from the *integral* generator lifts in
textual order, and the result is verified to be unimodular and to carry the source residue onto the target before it
is returned. For that trace the product is `E12^-1 * E21^-1 * E12^-2`, of determinant 1, sending `(1,0)` to
`(2,0)`. The reversed order sends it to `(3,0)`, so the order is not a convention to be guessed at: the
verification in `RadicalCongruenceImage.transporter` raises rather than returning a wrong carrier.

`Stabilizer(group.gap(), point, OnRight)` and `Size` of it work under the same typing, which is what
`RadicalCongruenceImage.stabilizer_order` uses (1, 2 and 10 for the three #41 specimens).

`Orbit`, `ActionHomomorphism`, `Image`, `RepresentativeAction`, `Stabilizer` and the word map were verified for
moduli `2, 4, 5, 6, 8` — prime and composite alike, `ZmodnZ(m)` being `GF(p)` exactly when `m` is prime.

## What is still genuinely port-side

Nothing in the action layer. The port supplies its own integral lifts and composes the word; that composition and
its verification are the port's obligation and remain in `RadicalCongruenceImage.transporter`.
