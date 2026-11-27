# Lorentzian perfect domains: mathematical and implementation certificates

This note records the proof obligations for `indefinite/lorentzian_cells.py`.
The algorithm authority is `polyhedral_common@1592b246`, especially
`src_lorentzian/lorentzian_perfect.h`, together with `Lorentzian.g` and the
translation map, section 12.  Mertens, *Automorphism Groups of Hyperbolic
Lattices*, arXiv:1303.3478v3 (`references/vendor/arxiv/1303.3478/hyperbolic.tex`)
supplies the dual-cone theorems.  Source code and independently recorded
fixtures, not agreement with this note's notation alone, decide correctness.

## 1. The two discrete sets and finiteness of the quotient

Let $b$ be integral, nondegenerate and of signature $(1,n)$, let $C$ be the
positive sheet, and let $G$ be the group of integral isometries preserving
it.  Identify the dual cone with $C$ using $b$.  In frame coordinates the
covectors of $L$ form a full lattice $M$ in $\mathbf Z^r$, invariant under
the adjoint action of $G$.  The total mode uses primitive nonzero points of
$M$ in the closed cone; isotropic mode uses only the primitive null points.
Primitivity changes no minimum and no cell.

### 1.1 Total mode

Mertens's set in (4), all nonzero integral points of the closed dual cone, is
admissible by his admissibility proposition (`propDadm`).  The same holds after intersecting with $M$:
if $m\mathbf Z^r\subseteq M$, first choose a point with pairing below
$\varepsilon/m$ and then multiply it by $m$.  Discreteness and invariance are
immediate.  Mertens's Theorem 3.15 and the forward implication of his
Theorem 2.3 give condition (2) for the full coordinate lattice and the same
cone action.  This condition concerns a single reduced fundamental domain,
not the chosen discrete subset.  It remains valid for $M\subseteq\mathbf Z^r$.
The reverse implication of Theorem 2.3 therefore gives finiteness for the
port's total set, and indeed for any discrete admissible invariant subset of
$M$ in the closed cone.  Primitive reduction needs a separate admissibility
argument, not merely inclusion.  For $x$ in the open cone and $d=kd_0$ with
$d_0$ primitive, strict duality gives $\sigma(x,d_0)=\sigma(x,d)/k\le
\sigma(x,d)$, because the fixed content $k$ is at least one.  Once
$\sigma(x,d)<\varepsilon$, the primitive representative works too.  No
arbitrary subset is claimed to inherit admissibility.

### 1.2 Null-only mode

The finiteness theorem for the full set does not directly include the
null-only set; admissibility must be proved.  Assume the lattice is
rationally isotropic.  Choose rational null $e,f$ in the same closed sheet
with $b(e,f)=1$: given a null $e$ and rational $g$ with $b(e,g)\ne0$, scale
$g$ to pairing one and subtract $b(g,g)e/2$.  Put $W=\{e,f\}^\perp$, negative
definite of dimension $d=r-2$, and use an integral basis of the saturated
rank-$d$ sublattice $W\cap L$.  All cube coordinates below are in that basis,
so every integral coefficient tuple $p$ is itself an element of $L$.

Every boundary ray other than $\mathbf R_{>0}e$ has a unique representative

$$v(t)=f+t-\frac{b(t,t)}2e,\qquad t\in W_{\mathbf R},$$

with $b(e,v(t))=1$, and

$$b(v(t),v(s))=-\frac{b(t-s,t-s)}2.$$

Given real $t$, a positive integer $Q$, and $d>0$, partition the unit cube
into $Q^d$ equal cubes.  Pigeonhole on $Q^d+1$ multiples $jt$ gives integers
$q,p$ with $1\le q\le Q^d$ and $\|qt-p\|_\infty\le1/Q$.  Choose one positive
integer $c$ clearing all coefficients of the rational quadratic
parametrization in the original lattice basis.  Then $cq^2v(p/q)$ is an
integral null vector of $L$ in the selected sheet; its covector consequently
lies in $M$.  It is nonzero because every $v(s)$ has coefficient $1$ on $f$ in the direct
decomposition $\langle f\rangle\oplus(\langle e\rangle\oplus W)$, while $s$ and
the subtracted term lie in spaces independent of $f$.  For $a>0$ the exact
pairing is

$$b(av(t),cq^2v(p/q))=-\frac{ac}{2}b(p-qt,p-qt),$$

which is nonnegative because $W$ is negative definite.  Since
$\|p-qt\|_\infty\le1/Q$, its value is at most $acK/Q^2$, where $K$ is the fixed
operator norm of $b$ from $\ell_\infty$ to $\mathbf R$ on $W$: the factor $q^2$ in the
scaled vector cancels exactly the factor $1/q^2$ in the squared approximation
error.  Dividing by its positive content
in $L$ preserves membership in $L$, hence membership of its covector in $M$,
and only decreases this nonnegative pairing.  At the ray of $e$ its own
primitive integral representative gives pairing zero; Mertens's admissibility
inequality has no strict-positivity requirement.  If $d=0$, both boundary rays
are rational and the same argument suffices.  Null-only points consequently
have arbitrarily small pairings with every boundary vector: this is exactly
admissibility.  Theorem 2.3, with condition (2)
already obtained above, now gives a finite quotient in isotropic mode too.

An anisotropic lattice has no such null set.  Its isotropic mode must refuse
input, as `_search_initial_vectors` does before any shell iteration; searching
forever is not an acceptable alternative.

## 2. Exact positive-vector shells

Let $t$ be an integral timelike direction, $d=\gcd_i b(t,e_i)$, and let $h$ be
a Bezout partner with $b(t,h)=d$.  The multiplier $m$ enumerates the fibre
$b(v,t)=md$.  Write $K=t^\perp$ and

$$w=(d/b(t,t))t-h,\qquad v=mh+x,\quad x\in K.$$

Then

$$b(v,v)=\frac{m^2d^2}{b(t,t)}+b(x-mw,x-mw).$$

The second term is nonpositive because $K$ is negative definite.  For the
positive-definite form $-b|_K$, a candidate is precisely a point of the affine
CVP ball about $mw$ of squared radius $m^2d^2/b(t,t)$.  Its norm is zero exactly
on the sphere and nonnegative exactly in the ball.  This is the formula in
`LORENTZ_FindPositiveVectorsKernel`.  `_PositiveVectorEnumerator` retains the
source's integral kernel basis, dual-LLL frame and candidate order, and checks
the reconstructed norm and pairing for every candidate.  Admissible input
ensures $d\ne0$; an arbitrary smaller pairing bound may legitimately return
the empty set.

Both modes use exact rational recursions, not a floating-point CVP decision.
The dual transforms are integral and inverse-transposes.  An exact LDL
factorization is certified before enumeration.  On a definite form, every
partial square contribution is nonnegative, so discarding a coordinate only
when its current contribution exceeds the remaining radius loses no solution;
the last coordinate is solved exactly (sphere) or over the exact integer
interval (ball).  `_fixed_distance_dual_lll_solutions` thereby enumerates the
whole sphere and `_bounded_distance_dual_lll_solutions` the whole ball, with
no PARI `qfcvp` precision or input-size hypothesis.  The admitted square locus
is then applied: it is always exactly `ZERO_SQUARES` or
`NONNEGATIVE_SQUARES`, as checked by the constructors.  Other loci are not
silently admitted.  All returned candidates are independently checked for
norm, pairing, and square locus.

With `OnlyShortest=false`, `vectors` concatenates every shell through the
pairing bound.  With `OnlyShortest=true`, it returns immediately after the
first nonempty shell.  All preceding shells are empty, and every vector in a
shell has the same pairing.  Therefore the legacy GAP wrapper, which asks for
the full ball and then selects the minimum positive pairing positions, returns
the same ordered shell as the C++ wrapper, which asks only for the shortest.
The initial search resumes with `next_multiplier` and the integral successor
of the previous pairing bound rather than rescanning an exhausted prefix.
It sets `max_scal=q(t)`.  Since $d\mid q(t)$, the vector $t$ itself lies in the
shell at the explicit multiplier $q(t)/d$, so total mode reaches a nonempty
shell within that many multipliers.  In isotropic mode, one of $v,-v$ for any
primitive null $v$ pairs positively with $t$ because $t^\perp$ is definite,
giving termination once the bound reaches that pairing.  This is a
semidecision with no a-priori size bound on the pairing; rational anisotropy is
therefore decided by the provider's exact Hasse--Minkowski/Simon procedure,
not inferred from an unsuccessful shell prefix.

## 3. The initial cell, flip bounds and flips

`initial_cell` follows `LORENTZ_GetOnePerfect`: find a timelike direction and
its first shell; while its rank is below the lattice rank, construct the affine
nullspace, select a negative outside direction, and cross the corresponding
wall.  Each returned candidate is accepted only after the augmented rank and
square-locus checks in `LorentzianPerfectCell`.

The source's sign hypotheses are checked, not assumed.  An initial outside
ray is negative by explicit certification.  For an interior facet, its span
contains a timelike point: if it contains a positive-square ray this is
immediate; if all its rays are null but two pair positively, their sum is
timelike; otherwise all rays lie in one null line, which is exactly the
exterior boundary case in rank two.  The Gram-dual direction normal is thus
negative definite.  For an algorithm-produced full-dimensional cell, the
supporting functional is its interior $D$-perfect point, so the base normal is
positive definite.  These are the preconditions $eNSPbas^2>0$ and
$eNSPdir^2<0$, checked in `_upper_bound` before either flipping scheme.  A
hand-built `LorentzianPerfectCell` is not covered by that theorem and will
instead fail loudly if these sign conditions are absent.

For independent affine normals $B,D$ write $N(s)=B+sD$ and let $u,v$ be their
Gram-dual directions.  `GetUpperBound` uses three possible exits:

* $N(s)_0=0$ when this gives a positive bound;
* the first doubled integer $s$ for which $N(s)$ is strictly negative;
* the least positive rational null parameter in their plane.

If both the constant and null bounds coincide, the ray is an inadmissible
isotropic endpoint and the direction is rejected.  With the displayed sign
hypotheses, the quadratic coefficient in $s$ is strictly negative, so the
doubling step terminates.  The selected upper value is checked to lie outside
the open positive cone.

`_kernel_flipping` preserves both source schemes.  In the first, bisection uses
continued-fraction middle values.  Equal shells move the lower bound, a shell
containing the critical set strictly is the successful flip, and a shell not
containing it enters the second scheme.  That scheme takes `ListTotal[0]` in
source CVP order, computes the exact wall parameter, and repeats until the new
shell contains the critical set and the selected vector (the stronger legacy
condition).  No noncritical-substitute or norm-sorted first vector is allowed.

Termination uses the actual bracket, not an appeal to the whole closed
interval.  Write $f(s)=A+2Bs+Cs^2$, with $A>0$ and $C<0$, and let $s_+$ be
its larger root.  When the trial square or pairing is nonpositive, the trial
parameter lies beyond $s_+$ and `upper` alone moves down.  Otherwise the
lower endpoint stays at or below $s_+$.  Every continued-fraction midpoint
lies in the middle third of its bracket, so the width contracts by at least
one third on every step.  Before any successful wall is found, the trial must
therefore become a point below $s_+$; that root and all comparisons are
computed exactly.

A non-blind direction has a neighbouring perfect point by Opgenorth's
Proposition 1.8.  It is enough to inspect the compact subsegment between this
first successful trial and the first wall.  That subsegment lies in the open
cone, where the face-to-face tessellation is locally finite, so the minimum
shell can change at only finitely many cells.  No boundary accumulation is
used.  At each point the candidate list is the minimum of finitely many
affine pairings.  The midpoint is the first continued-fraction convergent in
the middle third of the bracket midpoint.  The Stern--Brocot best-
approximation property makes it the least-denominator rational there, so any
fixed rational wall parameter left in the bracket must eventually be returned
exactly.

Each iteration therefore reaches, after finitely many steps, either a strict
shell extension (the successful flip) or the first wall where the old
critical set is no longer retained.  In the latter case scheme two starts at
that wall.  Entry into scheme two requires a noncritical first candidate,
equivalently a nonzero denominator; the source's enumerator order supplies
this precondition and the code checks it rather than skipping the candidate.
A new wall is then computed exactly from an integral expanded point.  Only
finitely many walls can occur in the compact subsegment, and each is
constructed rather than approached to within a tolerance, so that scheme
terminates too.  Global quotient finiteness is not used as an excuse for any
of these local loops.
`flip_across` computes the unique facet normal and orients it with an outside
ray.  The only totally isotropic exterior facet occurs in a Lorentzian plane;
there it is a cone-boundary face and has no neighbour, so it returns `None`.

## 4. Cell stabilizers and transporters

A cell's vectors span the ambient rational space.  Choose an ordered basis from
them.  Every integral cell transporter maps that basis to an $r$-tuple of target
vectors.  There are finitely many such tuples.  `cell_transporter`'s invariant
checks and colour refinements only reject candidates or select the order of a
complete backtracking enumeration; they never supply a negative answer by
themselves.  A leaf is accepted only after exact Gram congruence and a complete
configuration-membership check.  Consequently every transporter is found, and
`None` follows only after every possible basis image has been excluded.

Bucket keys impose only necessary conditions and never decide nonisometry by
hash inequality.  `_perfect_form_hash_key` contains multisets of diagonal and
off-diagonal Gram pairings, coordinate and Gram-covector contents, Smith
factors, and in the large-cell route row-sum and second-moment covector data.
An isometry maps each family to the corresponding one.  In particular, the
second moment changes by $M\mapsto A^{\mathrm t}MA$, Smith factors are
unchanged by a unimodular right action, and covector contents and norms are
preserved.  `_class_sum_pairing_signatures` sorts the existing invariant
colours before pairing each row with sums of colour classes, so its further
refinements are also canonical; with more than sixteen classes it simply
declines to refine.  Consequently every cell equivalent to the target's
initial cell has the target bucket key.  `isometry_to` tests that key for an
actual cell transporter at every bucket match and returns `None` only after
its complete source quotient has been exhausted without a witness.  The key
cannot discard a genuine match, and it cannot merge a false one without the
subsequent transporter certificate.

For a stabilizer, every image of every retained generating family is similarly
realized.  The pointwise stabilizer of an already assigned basis is computed by
the maintained finite-action operation, and the resulting Schreier generators
are checked as integral isometries and as configuration permutations.  This is
the complete setwise integral stabilizer, not merely a finite quotient or a
subgroup supplied by a graph automorphism.  The same transporter exhaustion
is used on every attempted traversal merge, so an equivalent neighbour cannot
be inserted as a duplicate orbit representative; that is the second input,
after flip effectivity, to traversal termination.

The diagonal-moment shortcut has its own exhaustive certificate.  If
$M=X^{\mathrm t}X=\operatorname{diag}(a,c,\ldots,c)$ with $a>c>0$, its minimum
vectors are the $r-1$ independent signed coordinate pairs.  Any automorphism
preserves them and their orthogonal complement, so the full definite group is
$C_2\times(C_2\wr S_{r-1})$; conversely every such signed permutation preserves
$M$.  Only generators additionally verified on the original indefinite Gram and
on every cell row are retained.  Failure of any certificate falls back to the
complete basis-image search, never to a subgroup assumption.  Coordinate-swap
reduction discards only an edge of an already connected transposition graph,
which cannot change the generated group.

## 5. Facets, their orbits and side pairings

The local H-description computes every supporting covector.  Its zero
incidence has rank $r-1$ exactly when it is a facet; the traversal checks that
condition.  The cell action on extreme rays is a finite permutation action.
Orbit construction by generator closure is therefore exhaustive, with exact
Schreier words.  Special crosspolytope and symmetric-chamber routes are accepted
only after their structure and generator certificates; otherwise the general
exact dual/facet walk runs.

For each facet-orbit representative, the flip constructs the neighbouring
cell.  An exact cell transporter identifies it with an existing representative
or establishes a new one.  The quotient stores that transporter.  Schreier
words recover one side pairing for every member of every retained facet orbit,
and their induced matrices are checked for integrality, unimodularity and Gram
preservation.  `perfect_domain_traversal` independently pairs the minimized
H-covectors with their incidences.  It rejects an unpaired facet unless the
span of its incident rays is totally isotropic, the exterior case above.  No
codimension-two crossing is needed for the generation theorem, so avoiding
ridges does not discard a facet required by the algorithm.

Opgenorth's tessellation theorem covers the entire cone, with a connected,
locally finite Voronoi graph.  Every cell class is therefore reachable from
the initial cell: no disconnected quotient component can be omitted.  Local
flip effectivity and complete transporter search make each expansion step
effective, while the invariant bucket key cannot hide an existing class.
Theorem 2.3 gives only finitely many orbits, hence the loop terminates.
Starting from one cell and adding cells only across quotient adjacencies makes
the retained transversal connected, as required by Mertens's Theorem 2.1.

## 6. Complete integral group generation and actual dispatch

Mertens's Theorem 2.1 (Opgenorth's graph-of-groups theorem) and the connected,
locally finite Voronoi graph show that the cell stabilizers and connecting
isometries generate the full sheet stabilizer $G=O^\Omega(L)$.  Equivalently,
a generic path in $C$ can be perturbed to cross only interiors and facets,
reducing an arbitrary displacement along the recorded quotient adjacencies.
The central $-I$ exchanges the two sheets, so adjoining it gives all of
$O(L)$, not just a discriminant quotient.

The nonreflective public dispatch reaches this construction through
`automorphism_group_generators` and the recursive Lorentzian leaf.  Indefinite
pre-reduction is an integral isometry $R\to L$; conjugation by it carries a
complete generating set of $O(R)$ to one of $O(L)$.  The signed-view functor
is also an integral isometry, so untwisting changes no generation assertion.
The rank-one and split-binary public branches are separate exact special
cases, and a definite rank-one input is not a Lorentzian cell-search input.

The reflective public dispatch instead uses the ordinary Vinberg normalizer
argument: $W(L)$ is normal in $O(L)$; a chamber can be returned by a wall word,
after which the residual is a chamber symmetry.  That gives a presentation of
the reduced sheet stabilizer and **still requires $-I$ for public
$O(L)$**.  See section 10 for the measured provider omission, which must not be
hidden by the cell-only proof above.

## 7. Primitive null marks and their global orbits

A future null ray of a polyhedral cone in $C$ is extreme: writing it as a sum
of two future causal points forces both points onto that same ray by
Lorentzian strict convexity.  Consequently the null extreme rays occur among
the certified configuration rows.  `cell_primitive_rays` checks the same
future component, filters exact norm zero and takes primitive representatives.
`local_stabilizer_orbits` then computes all orbits under the complete cell
stabilizer action.

On an interior shared facet, the side pairing bijectively identifies the null
marks incident to that facet with those of its image.  `transport_marks`
checks that this identification exists uniquely.  Gluing all local orbits by
all quotient adjacencies and taking connected components therefore identifies
exactly primitive null vectors in the same $O(L)$-orbit.  No other
identification is possible: every path from one local mark to another is a
word in cell stabilizers and side pairings, and conversely every isometry
carries cells and shared facets to cells and shared facets.

For the rank-13 and rank-17 specimens, the algorithm may stop earlier only
under additional, independently verified hypotheses.  The five rank-15
unimodular definite classes and their distinct norm-one counts exhaust the
possible quotient invariants in rank 17.  In rank 13, the unique admissible
definite genus and the Smith--Minkowski--Siegel mass certify exhaustion: the
sum of reciprocal full definite automorphism orders is $1891/367873228800$.
Distinct theta counts separate the displayed classes.  If any hypothesis or
mass identity fails, the catalogue path declines and the general complete
cell-orbit algorithm continues.

## 8. Timelike and other public marking classes

The source itself stops for nonzero norm in
`LORENTZ_GetOrbitRepresentative_Kernel`; the public timelike algorithm is a
new proof obligation, not an inference from the source's exception.  It uses
the barycentric subdivision of each cell.  That subdivision is invariant under
the finite cell stabilizer, and one maximal face flag per stabilizer orbit
covers the cell modulo that group.  A point's supporting-face chain and its
normalized coefficients in the associated simplex are intrinsic, so two
points are equivalent exactly when an element of the cell stabilizer carries
one chain and its coefficients to the other.

In signature $(1,n)$ the barycentre Gram entries of a simplicial piece have
one sign.  A positive fixed square bounds every coefficient on a ray of
positive square.  Descend to the support of the positive coefficients.  A
distinct causal ray pairs strictly positively with a null ray, so if the
null coefficient is nonzero then another supported coefficient is at least
one; their nonnegative cross term bounds it.  Any other positive coefficient
is likewise bounded by the fixed square.  `_simplicial_vectors_of_square`
enumerates precisely the resulting finite integral coefficient set,
including the saturation index of a nonsaturated simplicial sublattice, and
certifies every reconstructed norm.
The same construction on each adjacency facet, modulo the facet stabilizer,
connects all local mark orbits.  Connected components are again exactly the
global orbits.  Only primitive marks are retained.  A negative-square mark in
the normalized positive cone correctly raises `NoLocalMarkTheoremError`
rather than receiving an invalid local enumeration.

## 9. Independent executable evidence

At provider `e2803c36` and the port revision recorded in `results/commit`, the
actual frozen pytest nodes were run under `qc_sage_session` with JUnit output:

* `test_perfect_domains.sage`: 80 passed, 0 failed/skipped (53.31 s).
* `test_isometry.sage`: 155 passed, 0 failed/skipped (45.92 s), including all
  145 Lorentzian pairs and their public transporter checks.
* `test_vector_orbits.sage`: 63 passed, 0 failed/skipped (152.33 s).
* `test_vector_stabilizers.sage`: 7 passed, 0 failed/skipped (15.03 s).
* `test_orthogonal_group.sage`: 25 passed and 7 failed (16.32 s); every
  failure is the reflective public defect isolated in section 10.
* `IsotropicCases_8812` and `IsotropicCases_8819`: both frozen public nodes
  passed together in 19.34 s.  The rank-13/17 genus computations required
  Sage's supported external Maxima interface because the image's ECL FAS
  cache is broken.
* The exact historical PR #90 `test_cusp_catalogue_certificates.py` was read
  from commit `264e2fe5` and run without amendment at the current head: both
  nodes passed in 2.44 s.  These independently check the five rank-17 theta
  invariants and the rank-13 genera, four reciprocal automorphism orders, and
  the exact mass identity, rather than merely the one-bit existence nodes.

Local run files are under `.tmp/ci/results/issue15-e2803c36/` with port and
provider revision records.  Direct checks also established the three
cell-traversal mod-3 images by generator-by-generator mutual membership, one
primitive $U+E_8$ null orbit, one primitive class for `IsotropicCases_55`,
and the published setwise stabilizer image indices.  On Mertens's
nonreflective det-$155$ lattice, public `L.O()` returned 24 generators with
mutual generator membership against Mertens's complete $\Omega\cup\{-I\}$
family and equal mod-3 image order 48.  These are inputs to the proof, not
replacements for its completeness arguments.

## 10. Measured public reflective-full-group defect

The public boundary described in section 6 is not yet correct at provider
`e2803c36`.  `_edgewalk_orthogonal_group_generators` returns only simple-root
reflections and the reduced report's chamber symmetries, omitting $-I$.
There is no legitimate substitution of this reduced family for public $O(L)$.

For `LORENTZ_OrthogonalGroup0`, public `L.O().group_generators()` returns ten
root reflections; its mod-3 image has order $2579025599882610278400$, exactly
half the independently recorded image.  Every public generator belongs to the
reference image, but not every reference generator belongs to the public one.
`LORENTZ_OrthogonalGroup2` has the same defect.  Case 1's root reflections
happen already to generate the full reference image, so it does not test the
sheet-exchange branch.  The direct perfect-domain families match the reference
in all three cases, proving that this is specifically the provider's
reflective-to-public dispatch, not the cell algorithm.

The same defect appears in seven actual frozen OSCAR discriminant-image
nodes, all independently confirmed to use the reflective shortcut.  Their
computed/expected cardinalities are $8/16$, $32/64$, $48/96$, $12/24$,
$8/16$, $144/288$, and $8/16$ (`from_sage_3`, `10`, `11`, `12`, `14`, `24`,
and `25`).  The complete frozen `orthogonal_group` report is therefore 25
passed and 7 failed, not an OSCAR-fixture disagreement.  By contrast, the
nonreflective Mertens det-$155$ public dispatch has equal generator-by-
generator membership in both directions against Mertens's complete
$\Omega\cup\{-I\}$ reference set.

The provider must adjoin the actual `-L.Aut().identity()` to its public
reflective generator family (duplicate generators are harmless).  It must not
insert $-I$ into the reduced Allcock report, whose chamber-symmetry field has
a different meaning.  For both defective stabilizer specimens, augmenting the
public family by this exact element restores both directions of the original
mod-3 membership test.  The port provides no falsified report and no
dependency edit in its private runtime.  This omission blocks a completion
claim for the public contract until the maintained provider boundary is
repaired.
