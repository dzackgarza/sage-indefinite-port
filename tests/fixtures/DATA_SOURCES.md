# Data Sources and Provenance

This document records the exact provenance, source files, authors, publications, and extraction methods for every dataset in [`tests/fixtures/`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures).

* * *

## Datasets

### 87 Numerical Enriques Polarizations

- **File**: `tests/fixtures/enriques_87_polarizations.json`

- **Primary Source**: Mathieu Dutour Sikirić and Klaus Hulek, *Moduli of polarised Enriques surfaces — computational aspects*, J. London Math. Soc. 2024, [arXiv:2302.01679](https://arxiv.org/abs/2302.01679), Tables `table_subgroups1` and `table_subgroups2`; vendored at `references/vendor/arxiv/2302.01679/Enriques_compu_rev.tex`.

- **Extraction**: `uv run references/extract/dh23_polarizations.py` parses both tables from the vendored TeX. Each record carries the TeX line it came from.

- **Content**: for each of the 87 conjugacy classes of groups $\Gamma_h$: the subset $S$ of the $E_{10}$ diagram, $\#S$, $|\bar\Gamma_h|$, the isotropic line, plane and flag counts $\#I_1, \#I_2, \#I_{12}$, the degree of the smallest realization $h_{\min}$, and $\phi(h_{\min})$.

- **Check**: the script asserts cases 1–87 in order and that case 87 has 528 lines and 24242 planes, which the paper also states independently in its text (lines 1168–1170).

### Orbits of Primitive Vectors of Norm 2..72 in U + E8(-1)

- **File**: `tests/fixtures/enriques_polarization_orbits.json`

- **Source**: Dutour Sikirić–Hulek, arXiv:2302.01679, Table `ListNumberPolarizationsModuli`, vendored TeX. For each $2d = 2, \dots, 72$ the table gives $\#h$, the number of orbits of primitive vectors $h$ with $h^2 = 2d$ in $\mathrm{Num}(S) \cong U \oplus E_8(-1)$. The paper notes these agree with the CDGK appendix. It also gives the number of conjugacy classes of the groups $\bar\Gamma_h$, and cumulative counts for both.

- **Extraction**: `uv run references/extract/dh23_norm_orbits.py`. It checks that $2d = 2g-2$, that the cumulative row is the running sum of $\#h$ (ending at 312), and that $\#\bar\Gamma_h \le \#h$.

### Unpolarized Enriques Boundary Strata and Stabilizers

- **File**: `tests/fixtures/unpolarized_enriques.json`

- **Primary Source**: Dutour Sikirić & Hulek, arXiv:2302.01679, subsection "The Tits building" (vendored TeX `references/vendor/arxiv/2302.01679/Enriques_compu_rev.tex`, lines 1159–1170), and Table `table_subgroups2` case 87 (line 814) for the flag count.

- **Extraction**: `uv run references/extract/dh23_unpolarized.py` builds $N = U \oplus U(2) \oplus E_8(-2)$, takes the representatives from the text, and asserts $\det N = 1024$ and that every line and plane representative is totally isotropic.

- **Content**:

  - Gram matrix of $N$, with $(e_1, e_2)$ and $(e_3, e_4)$ the standard bases of $U$ and $U(2)$.

  - Isotropic line representatives $\mathbb{Z}e_1$, $\mathbb{Z}e_3$, with discriminant-image stabilizer indices 1 and 527.

  - Isotropic plane representatives $\mathbb{Z}e_1 + \mathbb{Z}e_3$ and $\mathbb{Z}(2e_1 + 2e_2 + w) + \mathbb{Z}e_3$, with indices 527 and 23715. Here $w = \alpha_0 + \alpha_2$ has norm 4 in $E_8$, i.e. $w^2 = -8$ in $E_8(-2)$, as the paper requires.

  - The stable group $\widetilde O^+(N)$ has 528 lines and 24242 planes (text) and 72199 flags (table).

### 8,821 Reflective Lorentzian Forms and Their Isotropy

- **Files**: `tests/fixtures/reflective_forms_8821.json`, `tests/fixtures/isotropic_cases_8821.json`

- **Upstream**: `references/vendor/polyhedral_common@1592b246/CI_tests/20_Reflective/ListReflect` (Gram matrix and number of simple roots) and `references/vendor/polyhedral_common@1592b246/CI_tests/DATA/IsotropicCases` (whether an isotropic vector exists), the same 8,821 lattices in the same order; the script asserts the pairing row by row.

- **Note**: the rank-3 part is Allcock's classification of reflective Lorentzian lattices of rank 3 (vendored separately with its explicit simple roots under `references/vendor/GeometryDatabase_Rank3_Lorentzian_lattices@63a9067a/`).

- These are the reference implementation's own recorded outputs (polyhedral_common, Mathieu Dutour Sikirić); agreement with them is acceptance for this port. Extracted by `references/extract/polyhedral_common_ci.py` (run under Sage's Python; GAP literals are evaluated with libgap). Each record cites its vendored file and index.

### Allcock's Reflective Lorentzian Lattices of Rank 3

- **File**: `tests/fixtures/allcock_rank3_reflective.json`

- **Source**: Allcock's classification of the reflective Lorentzian lattices of rank 3 (arXiv:1010.0486, 1111.1264, both vendored), as published in `references/vendor/GeometryDatabase_Rank3_Lorentzian_lattices@63a9067a/RK3_all`. It has 8,595 lattices, each with its elementary divisors, explicit simple roots and Weyl-group id. The convention is signature (2,1) with roots of positive norm, the same as `ListReflect`.

- **Extraction**: `uv run references/extract/allcock_rank3.py`. It asserts that every simple root has positive norm and gives an integral reflection. It also checks that every Gram matrix appears in `reflective_forms_8821.json` with the same number of simple roots, which confirms that `ListReflect`'s rank-3 part is exactly this classification.

### 145 Lorentzian Equivalence Instances

- **File**: `tests/fixtures/lorentzian_equivalence_145.json`

- **Upstream**: `references/vendor/polyhedral_common@1592b246/CI_tests/28B_LorentzianPerfStabEqui/TestCasesEqui.tar.gz` (extracted copy under `x/TestCasesEqui/`, checked against the archive). 145 pairs (`mat1`, `mat2`) of rank 10, signature (1,9), each with a recorded transporter; the script asserts `witness * mat1 * witness^T == mat2` for all 145.

- These are the reference implementation's own recorded outputs (polyhedral_common, Mathieu Dutour Sikirić); agreement with them is acceptance for this port. Extracted by `references/extract/polyhedral_common_ci.py` (run under Sage's Python; GAP literals are evaluated with libgap). Each record cites its vendored file and index.

### Lorentzian Stabilizer Generators

- **File**: `tests/fixtures/lorentzian_stabilizers_cases.json`

- **Upstream**: `references/vendor/polyhedral_common@1592b246/CI_tests/28B_LorentzianPerfStabEqui/TestCasesStab.tar.gz` (extracted under `x/TestCasesStab/`): 3 rank-10 lattices with generator sets of 15, 30 and 14 matrices; the script asserts every generator preserves its form.

- These are the reference implementation's own recorded outputs (polyhedral_common, Mathieu Dutour Sikirić); agreement with them is acceptance for this port. Extracted by `references/extract/polyhedral_common_ci.py` (run under Sage's Python; GAP literals are evaluated with libgap). Each record cites its vendored file and index.

### 103 Root Systems of Reflective Forms

- **File**: `tests/fixtures/root_systems_103.json`

- **Upstream**: `references/vendor/polyhedral_common@1592b246/CI_tests/01_RatIntAutomorphy/ListSimpleRootSystem_4_56_X_5_47`, the simple roots of the rank-4 (56) and rank-5 (47) entries of `ListReflect`, in the same order. Each record now carries its Gram matrix from `ListReflect`; the script asserts the root count matches `n_simple` and that every root defines an integral reflection.

- These are the reference implementation's own recorded outputs (polyhedral_common, Mathieu Dutour Sikirić); agreement with them is acceptance for this port. Extracted by `references/extract/polyhedral_common_ci.py` (run under Sage's Python; GAP literals are evaluated with libgap). Each record cites its vendored file and index.

### Classified Simplices (Dimensions 5, 6, 7)

- **File**: `tests/fixtures/classification_simplices.json`

- **Upstream**: `references/vendor/polyhedral_common@1592b246/CI_tests/DATA/ClassificationSimplices{5,6,7}`: 16 simplices (2 in dimension 5, 3 in 6, 11 in 7), the inputs of `01_RatIntAutomorphy/ProcessExamples.g`. Upstream records no outputs for them; its check is exact and mathematical (the double cosets partition the rational group: disjoint, sizes summing to its order).

- Extracted by `references/extract/polyhedral_common_ci.py`; each record cites its vendored file.

### 18 Finite Double-Coset Instances

- **File**: `tests/fixtures/double_coset_cases.json`

- **Upstream**: `references/vendor/polyhedral_common@1592b246/CI_tests/DoubleCosets/DBL/`: 18 inputs (degrees 55–96). Each file holds two permutation groups (0-based points) and a set of 0/1 vectors; the two group orders and the vector count are encoded only in the filename. The script parses each file exactly and assigns `group_g` (the big group) and `group_h` (the subgroup) by computing both orders with libgap and matching the filename. The previous fixture had the two groups swapped in all 18 cases.

- **Recorded outputs**: none. Acceptance uses the mathematical identities of a double-coset decomposition (disjoint union, sizes summing to the group order).

- Extracted by `references/extract/polyhedral_common_ci.py`; each record cites its vendored file.

### 40 Lorentzian Perfect Domain Forms

- **File**: `tests/fixtures/lorentzian_perfect_domains.json`

- **Upstream**: `references/vendor/polyhedral_common@1592b246/CI_tests/28B_LorentzianPerfStabEqui/Result_Enumeration`: for 40 lattices of ranks 3–6, the exact number of perfect-domain orbits in the "isotropic" and "total" modes.

- These are the reference implementation's own recorded outputs (polyhedral_common, Mathieu Dutour Sikirić); agreement with them is acceptance for this port. Extracted by `references/extract/polyhedral_common_ci.py` (run under Sage's Python; GAP literals are evaluated with libgap). Each record cites its vendored file and index.

### K3 Lattice and Degree-2 / Degree-4 K3 Baily–Borel Boundaries

- **File**: `tests/fixtures/k3_modular_strata.json`

- **Sources** (vendored TeX; the extraction asserts that every transcribed label appears on the cited lines):

  - K3 lattice $3U \oplus 2E_8(-1)$: one orbit of primitive vectors per represented norm, by Eichler's criterion, Gritsenko–Hulek–Sankaran arXiv:1012.4155, Lemma `lem:eichler` (`references/vendor/arxiv/1012.4155/main.tex`).

  - Degree 2: four Type II curves meeting in one Type III point (Laza arXiv:1205.3144, theorem citing Scattone §6.2), with Type II root types $A_{17}$, $E_8^2+A_1$, $D_{16}+A_1$, $E_7+D_{10}$ (Shah's table as given there).

  - Degree 4: nine Type II curves with generalised types $A_{11}+E_6$, $A_1+A_1+A_{15}$, $D_8+D_8+\langle-4\rangle$, $D_{12}+D_5$, $D_{16}+\langle-4\rangle$, $E_7+E_7+A_3$, $E_8+E_8+\langle-4\rangle$, $E_8+D_9$, $D_{17}$ (Jones arXiv:2502.04301, theorem citing Scattone §6.3).

- **Extraction**: `uv run references/extract/k3_strata.py`.

- **Removed**: the "$2^{\omega(d)-1}$ zero-dimensional cusps" family attributed to Attwell-Duval. No source for it is vendored or known, and for squarefree $d$ the Eichler-criterion count is one cusp, which contradicts it. The previous degree-4 labels included $A_{17}$ and $D_{10}+E_7$ (degree-2 labels) and used $A_1$ where the classification has $\langle-4\rangle$.

### Baily–Borel Boundaries of the D-Tower (Laza–O'Grady)

- **File**: `tests/fixtures/dtower_boundaries.json`

- **Source**: Laza–O'Grady, arXiv:1801.04845 (vendored TeX). $\Lambda_N = U^2 \oplus D_{N-2}$ with $D$ negative definite (line 625); $\Gamma(N) = O^+(\Lambda_N)$ for the cases recorded (line 642). The diagram `bbpicture` (lines 4629–4631) gives the Baily–Borel boundaries of $\mathcal F(9)$, $\mathcal F(10)$ and $\mathcal F(11)$:
  - $N=9$: 1 Type III point and 1 Type II curve ($D_7$), 1 incidence;
  - $N=10$: 2 points and 2 curves ($E_8$, $D_8$), 3 incidences;
  - $N=11$: 1 point and 2 curves ($E_8\oplus D_1$, $D_9$), 2 incidences.

- **Not included**: Type II counts for other $N$. The paper's theorem for $3 \le N \le 20$ refers to a table under the label `tabletype2`, but that label belongs to a different table, the Type II components of $\widehat{\mathcal F}$. No per-$N$ counts are stated, so none are reconstructed.

- **Extraction**: `uv run references/extract/laza_ogrady_dtower.py`. It asserts every label on its cited line.

### Orbits in I_{2,10} (Allcock)

- **File**: `tests/fixtures/allcock_i_2_10_orbits.json`

- **Source**: Allcock, math/9905166 (vendored TeX). The Enriques period lattice $\hat K \cong I_{2,10}$ (lines 177–179) and $\Gamma = \operatorname{Aut}\hat K$. Corollary 3: one orbit of norm $-1$ vectors. Corollary 4: two orbits of primitive isotropic vectors ($v^\perp/v \cong I_{1,9}$ or $II_{1,9}$) and two orbits of isotropic planes ($V^\perp/V \cong E_8(-1)$ or $I_{0,8}$).

- **Extraction**: `uv run references/extract/allcock_enriques_period.py`. It asserts every statement on its cited lines.

### Dawes Tits Buildings and Index Chain

- **File**: `tests/fixtures/dawes_buildings.json`

- **Sources** (vendored TeX; the extraction asserts each fact on its cited lines). Dawes assumes root lattices negative definite (arXiv:2205.10601, line 158):

  - arXiv:2205.10601, section "Examples":
    - the stable orthogonal group of $2U \oplus A_2$ has a building with 1 point, 1 curve and 1 edge;
    - $O^+$ and the stable group of $2U \oplus \langle-2\rangle \oplus \langle-6\rangle$ share a building with 2 points, 2 curves and 3 edges;
    - the stable groups of $2U(2)\oplus A_2 \subset U\oplus U(2)\oplus A_2 \subset 2U\oplus A_2$ and $O^+(2U\oplus A_2)$ form a chain with indices 20, 27 and 2 (total 1080).

  - arXiv:2108.06236, Theorem `L2boundarythm`: for $L_2 = 2U\oplus\langle-2\rangle\oplus\langle-6\rangle$ and $\Gamma_2$, the boundary has 3 points, 2 curves and 4 incidences.

- **Not included**: the building of $2U(2)\oplus A_2$ (Figure `2u2a2building`). The text states no counts for it, and reading them off the TikZ drawing would mean interpreting the figure.

- **Extraction**: `uv run references/extract/dawes_buildings.py`.

- **Replaced**: the previous fixture recorded no outputs, and its third lattice $U\oplus U(2)\oplus A_2$ was not a published example.

### Vinberg Roots, Discriminant Images and Isometry Groups (OSCAR)

- **File**: `tests/fixtures/oscar_lattice_oracles.json`

- **Sources**: OSCAR's own tests (an independent implementation), vendored at `references/vendor/Oscar.jl@d135b70b/test/`:

  - `NumberTheory/vinberg.jl`: four reflective lattices with their numbers of simple roots, plus 2 cusps for one case, 6 Coxeter-diagram edges for another, and the explicit roots of the test's chamber for a third.

  - `Groups/spinor_norms.jl`: 28 rank-3 lattices from the `from_sage` table (10 run in OSCAR's CI, 18 commented out there but still recorded). Each gives $|O(q_L)|$ and the orders of the images in $O(q_L)$ of $O(L)$ and of its real-spinor-norm subgroup. Two further cases have a recorded bijective image map.

  - `Groups/isometry_group.jl`: $|O(U)| = 4$, and a vector stabilizer of order 2.

- **Extraction**: `references/extract/oscar_lattice_tests.py`, run under Sage's Python. It asserts each value on its source line, that every lattice is indefinite, and that explicit roots give integral reflections.

### Centralizers of Finite-Order Isometries (OSCAR)

- **File**: `tests/fixtures/isometry_centralizers.json`

- **Source**: OSCAR's own tests, an independent implementation, vendored at `references/vendor/Oscar.jl@d135b70b/test/NumberTheory/QuadFormAndIsom/` (`lattices_with_isometry.jl` and `enumeration.jl`, line ranges per record).

- **Content**:

  - Five indefinite lattices with an isometry of order 4, 5 or 6. For each, OSCAR records the order of the centralizer's image in $O(q_L)$ (72, 2, 96, 24192), or that the image is all of $O(q_L)$.

  - For a signature-(1,9) genus, OSCAR records 11 classes of pairs $(L, f)$ with characteristic polynomial $(x-1)^4(x+1)^6$ (9 locally).

- **Extraction**: `references/extract/oscar_centralizers.py`, run under Sage's Python. It parses the Julia matrix literals mechanically and reads each recorded value from its `@test` line. It computes $B G B^T$ and the isometry in the lattice basis, and asserts integrality, form preservation and indefiniteness. The definite order-600 case is excluded by that check.

- **Replaced**: the previous `centralizer_involutions.json` held derived identities, not data.

### Indefinite Isometry Pairs (Indefinite.jl and Hecke)

- **File**: `tests/fixtures/indefinite_isometry_pairs.json`

- **Sources**: `references/vendor/Indefinite.jl@374a5ebb/TestLor/` (three pairs that `test_gap.jl` tests for equivalence: $U\oplus3\langle-1\rangle$, $U\oplus E_8(-1)$, $U\oplus U(2)\oplus3\langle-2\rangle$), and `references/vendor/Hecke.jl@e2ab5716/test/QuadForm/Quad/ZLatticeAutIso.jl`, testset "isometry testing":
  - $U(2)\oplus A_2$ with the recorded witness $u$;
  - $U$ against $U(2)$ (non-isometric);
  - a recorded isometric $3\times3$ pair.

- **Certificates** (checked by the extraction, independently of the port):
  - a verified witness $W\,G_1\,W^T = G_2$;
  - or equal genera with no spinor generators, which for indefinite rank $\ge 3$ means one isometry class (Eichler), computed with Sage's genus code;
  - or, for the non-isometric pair, different determinants.

- **Extraction**: `references/extract/isometry_pairs.py`, run under Sage's Python.

### 6 Metamorphic Indefinite Forms

- **File**: `tests/fixtures/ci_indefinite_comp.json`

- **Upstream**: `references/vendor/polyhedral_common@1592b246/CI_tests/19_IndefiniteComp/AllTests.g` (`FullTest`), whose six families are built by upstream's `GetGramMatrixFromList` in `references/vendor/polyhedral_common@1592b246/CI_tests/common.g`. The script loads `common.g` into libgap and calls that function, so each Gram matrix is upstream's own (root lattices positive definite); signatures are computed from it. The previous fixture recorded invented signatures such as (2,4) for U+2U+A2, whose upstream Gram has signature (4,2).

- **Recorded outputs**: none; the upstream suite is metamorphic (stabilizer generators preserve the form, a random conjugate is found isometric, and vector, plane and flag orbit counts agree between a form and its conjugate).

- Extracted by `references/extract/polyhedral_common_ci.py`; each record cites its vendored file.
