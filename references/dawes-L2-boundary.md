# The boundary for L₂ and Γ₂

## Oracle conflict

For L = 2U ⊕ ⟨−2⟩ ⊕ ⟨−6⟩, the correct numbers of rational
isotropic line, plane and flag orbits under Γ₂ are (2, 2, 3).
The frozen `dawes_L2_gamma2` fixture requires (3, 2, 4).
Its cited theorem conflicts with the preceding lemmas in the same source.

The pinned source is `vendor/arxiv/2108.06236/main.tex`:

| Source | Claim |
| --- | --- |
| Lines 204–212 | Definition of L₂ and Γ₂ |
| Lemma `gamma2gplem`, lines 984–998 | Γ₂ = O⁺(L) |
| Lemma `L2planeslem`, lines 999–1014 | Two primitive null-plane orbits |
| Lemma `L2lineslem`, lines 1015–1026 | Two primitive null-vector orbits |
| Theorem `L2boundarythm`, lines 1044–1047 | Three boundary points and four incidences |

The [arXiv v1 source](https://arxiv.org/html/2108.06236v1#S3.SS5)
has the same conflict. The argument below determines the orbit counts
independently of the port's orbit enumeration.

## Lines

Use the ordered basis e₁,f₁,e₂,f₂,v,w, with (eᵢ,fᵢ)=1,
v²=−2 and w²=−6. Then

    D(L) = Z/2 ⊕ Z/6,
    q(a,b) = −a²/2 − b²/6 modulo 2Z.

Enumeration of these residue classes gives exactly the isotropic elements
(0,0) and (1,3). The lattice contains 2U. The Eichler criterion therefore
gives at most one stable SO⁺ orbit for each of these classes.
Both occur: e₁ has divisor 1, and

    l = 2e₂ + 2f₂ + v + w

is primitive, null and has divisor 2. Their divisors distinguish their
O(L) orbits. The stable SO⁺ group is contained in Γ₂. Thus Γ₂ has
exactly two line orbits. Passing from primitive vectors to unoriented
lines preserves this count.

## Planes and flags

The source normal-form lemma and `L2planeslem` give the two plane
representatives E₁=⟨e₁,e₂⟩ and E₂=⟨e₁,l⟩. Their sublattices are
primitive: the displayed basis coordinate matrices have a unit minor.

The SL₂(Z) action on E₁ extends by the contragredient action on
⟨f₁,f₂⟩ and the identity on ⟨v,w⟩. These extensions lie in stable
SO⁺(L). They act transitively on primitive lines in E₁. Hence E₁
contributes one flag orbit, incident with the divisor-one line orbit.

For E₂, take the integral Eichler transvections

    u = t(e₁,−e₂),  z = t(l,−f₁).

Their actions on the ordered basis (e₁,l) are

    u = [[1,2],[0,1]],  z = [[1,0],[1,1]].

Both preserve L, act trivially on D(L), and lie in the identity
component over R. Thus both lie in Γ₂. Their inverses replace the
second transvection parameter by its negative.

A primitive vector a e₁+b l has gcd(a,b)=1 and divisor gcd(a,2b).
There are exactly two possible divisors, distinguished by the parity
of a. The operations induced by u and z are

    (a,b) ↦ (a+2kb,b),  (a,b) ↦ (a,b+ka),  k∈Z.

They reduce any primitive pair to one of (±1,0) or (0,±1): first
reduce b modulo a to |b|≤|a|/2, then reduce a modulo 2b to
|a|≤|b|. Repeat while both coordinates are nonzero. The largest
absolute coordinate strictly decreases after each complete cycle.
The parity of a fixes which terminal pair occurs. The word
(u z⁻¹)² acts as −I and identifies the two signs. Consequently
the stabilizer of E₂ has exactly two line orbits, with divisors
1 and 2. E₂ contributes two flag orbits.

The quotient incidence graph therefore has a divisor-one point joined
to both planes and a divisor-two point joined only to E₂. Its
numbers of points, curves and edges are (2,2,3).

## Executable witnesses

`tests/test_dawes_L2_boundary_witnesses.py` checks the complete finite
discriminant calculation, the primitive null vector and its divisor,
the two transvection actions, preservation of all basis pairings,
and inverse identities. It uses the owned lattice and isometry
operations. The plane classification uses the cited source normal
form; the finite witnesses alone do not replace that classification.

The frozen fixture remains subject to the correction procedure in #33.
