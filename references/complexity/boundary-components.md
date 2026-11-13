# Baily–Borel boundary components and modular groups: complexity

## Algorithm and pinned source

For a lattice `L` of signature `(2, n - 2)` and a finite-index subgroup
`Γ <= O(L)`, the Type III points, Type II curves and their incidences are the
`Γ`-orbits of primitive isotropic lines, planes and line-in-plane flags.

| Step | Port or provider | Pinned source |
| --- | --- | --- |
| `O(L)`-classes of isotropic rank-`k` sublattices and flags | `INDUCTION.isotropic_flag_bases` (`isotropic_induction.py:471`), `primitive_marking_classes` (`marked_lattices.py:617-660`) | `INDEF_FORM_GetOrbit_IsotropicKplane` (`polyhedral_common@1592b246/src_indefinite/INDEF_FORM_GetOrbit_IsotropicKplane.cpp`) |
| Stabilizers of a line and a plane | `line_automorphism_generators`, `plane_automorphism_generators` (`isotropic_induction.py:744-750`) | `INDEF_FORM_StabilizerIsotropicPlane.cpp` |
| Split along `Γ` | `split_orbit` (sage-extensions `categories/lattice_morphisms.py:4638`) | double cosets of finite quotients |
| Enriques groups `Γ_h` | `coinvariant_extension_subgroup` (sage-extensions `categories/lattice_centralizers.py:1110`) | Dutour Sikirić and Hulek, arXiv:2302.01679 |

## Isotropic classes

- Lines. The stabilizer of a primitive isotropic line `I = Ze` is computed
  through `e^⊥/e`, of rank `n - 2` and signature `(1, n - 3)`. That lattice
  is Lorentzian, so its cost is one perfect-domain traversal
  (`orthogonal-group-generators.md`, "Lorentzian route").
- Planes. For a primitive isotropic plane `E`, `E^⊥/E` is negative definite
  of rank `n - 4`. The `O(L)`-classes of planes map onto classes in the genus
  of `E^⊥/E` (Scattone, Mem. AMS 374, §5.4; for squarefree `k` the map
  `E -> E^⊥/E` is a bijection onto the genus of `<-2k> ⊕ 2E8(-1)`,
  arXiv:2502.04301, lines 273–282). The cost is the genus enumeration:
  Kneser `2`-neighbours until the mass formula is met, `h` classes,
  `h^2` definite isometry tests of rank `n - 4`.
- Flags. A flag `I ⊂ E` is a line class inside each plane class; the cost is
  the line computation inside each `E^⊥`.

## Splitting

`split_orbit` realizes the `Γ`-orbits in the `O(L)`-orbit of a marking `m`
as the double cosets `ρ(Γ) \ ρ(O(L)) / ρ(O(L)_m)`, where `ρ` is the
character map to `O(q_L)` (`lattice_morphisms.py:4638-4700`). It computes the
`ρ(Γ)`-orbits on the finite set `ρ(O(L))/ρ(O(L)_m)`. The size of that set is
at most the number of isotropic configurations of the same type in `A_L`:
`<= |A|` for lines and `<= |A|^2` for planes and flags.

`splitting_isometries_from_stabilizer_image` (sage-extensions
`categories/orthogonal_quotients.py:1600-1630`) runs GAP `DoubleCosets` in the
permutation model of the image. It then lifts exactly one element for each
double coset, that is, for each output representative. It does not lift a
transversal of `[O(q) : ρ(Γ)]`. One lift has these parts:

- `PreImagesRepresentative` of the free-group epimorphism (`_GapWordSection`,
  `orthogonal_quotients.py:60-191`) on the existing stabilizer chain. This
  costs `O(w + deg)` permutation operations.
- `_evaluate_letters`: `w` products of `n x n` integer isometries. Each
  product is one FLINT product, `U_PRODUCT = 1e-5` s at `n = 12`.
- `_checked_lift`: one character image, `ell * n` element operations.
- One apex action on the marking.

The word length follows from the stabilizer chain of the image `Q` in degree
`deg`. The base length is `b = ceil(log2|Q| / log2 deg)`. Each level adds at
most one Schreier-tree path of depth `<= 2 ceil(log2 deg)`, and a strong
generator is on average a word of length `<= 2.5` in the supplied family. So
`w <= 5 b ceil(log2 deg)`.

`T_lift = (w + deg) U_PERM + w U_PRODUCT + (ell + 2) n U_ELT`

`T_split = |G/G_m| * g * U_PERM + (number of Γ-orbits) * T_lift`.

For the Enriques lattice `N` the image is of `O^+(10, 2)` type on the
`deg = 527` singular points, with `log2|Q| = 45.4`. This gives `b = 6`,
`w = 300` and `T_lift = 0.0032` s.

The theory of the character quotient is in
`discriminant-groups-and-congruence-orbits.md` (owner F), "Orbit splitting
across markings".

## Dawes buildings and the D-tower

The Dawes lattices (`2U ⊕ A2`, `2U(2) ⊕ A2`, `2U ⊕ <-2> ⊕ <-6>`, `L_2`) have
`|A| <= 48`, so the split sets are small and the cost is the isotropic
classes. The D-tower lattices `Λ_N = U^2 ⊕ D_{N-2}(-1)` have `|A| = 4`
(`|A| = 4` for `D_m`, `m >= 2`; `<-4>` for `N = 3`). Their cost is dominated
by the Lorentzian traversal of `e^⊥/e = U ⊕ D_{N-2}(-1)` and by the genus of
`D_{N-2}`, whose class number is 1 except for the extra classes at
`N = 14` and the even unimodular classes at `N ≡ 2 mod 8` (Laza and O'Grady,
arXiv:1801.04845, lines 4460 and 4604–4616).

## K3 lattices of degree 2 and 4

`Λ_{2k} = 2U ⊕ 2E8(-1) ⊕ <-2k>` has rank 21 and `|A| = 2k`. The plane
classes are the genus of `<-2k> ⊕ 2E8(-1)` of rank 17: `h = 4` for `k = 1`
and `h = 9` for `k = 2` (Scattone, §§5.4, 6). Each test costs about `h^2`
rank-17 definite isometry tests.

## Enriques modular groups

`N = U ⊕ U(2) ⊕ E8(-2)` has signature `(2, 10)` and `A_N = (Z/2)^10` with the
form of `O^+(10, 2)`. The setup (primitive extension of the Enriques
involution and an isometry `E10(2) -> M`) is cached once per process. Per
case:

- `Γ_h`: the stabilizer of `h` in `O(M)` (a rank-9 definite complement), its
  image in `O(q_M)`, transport by the glue, and the preimage in `O(N)`. The
  stable subgroup is generated by Eichler transvections, so no enumeration of
  `O(q)` (order `2^21 3^5 5^2 7 17 31 ≈ 4.6 * 10^13`) is needed.
- Orbit counts. `O(N)` acts on the singular points of `A_N`
  (`(2^5 - 1)(2^4 + 1) = 527`), on its totally singular lines
  (`527 * 135 / 3 = 23715`) and on its point-line flags (`527 * 135 = 71145`).
  The `ρ(Γ_h)`-orbits on these sets give the recorded line, plane and flag
  counts. The number of output representatives is the recorded sum, up to
  `96969` for the case of smallest `Γ_h`. The section lifts dominate there.

## Shared work across the Enriques cases

All 87 cases use the same lattice `N`. Its `O(N)` generators, its `O(N)`-classes
of lines, planes and flags, the permutation model of `O(q_N)` and the transport
setup are computed once per session and charged to `case_1`, the first case in
collection order. Each case pays only its own work:

- `Stab(h)` in `O(M)` through the definite rank-9 complement;
- its image in `O(q_M)`;
- at most 20 section lifts for the preimage generators;
- the image of `Γ_h` from at most 40 generators, each fixed by the images of
  the `ell = 10` generators of `A_N`;
- the GAP double cosets and one `0.0032` s lift per output representative.

The cost per case is therefore linear in the recorded number of orbits. It is
largest for the case of smallest `Γ_h`, which has 96969 representatives.

## Known worst cases

- A small `Γ`: the number of output representatives grows as
  `[O(L) : Γ]` up to the size of the configuration sets in `A_L`.
- Planes with a definite complement genus of large class number.

## Unit costs of backend primitives

The unit costs are those of `orthogonal-group-generators.md`, "Unit costs of
backend primitives". Definite complements over `ZZ` use PARI `qfauto` and
`qfisom`, priced by `T_definite` of "Definite kernel" there.
