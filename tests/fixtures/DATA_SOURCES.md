# Data Sources and Provenance

This document records the exact provenance, source files, authors, publications, and extraction methods for every dataset in [`tests/fixtures/`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures).

* * *

## 1. Published Mathematical Tables and Theorem Oracles

### A. 87 Numerical Enriques Polarizations

- **File**: [`tests/fixtures/enriques_87_polarizations.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/enriques_87_polarizations.json)

- **Primary Source**: Mathieu Dutour Sikirić and Klaus Hulek, *On the classification of numerical Enriques surfaces and their moduli spaces*, [arXiv:2302.01679](https://arxiv.org/abs/2302.01679) [math.AG], 2023, Tables 1 and 2.

- **Content**: 87 rows containing case numbers, minimal degrees $d_{\min}$, polarization orbit lengths, zero-dimensional cusp orbit counts $\#\mathcal{I}_1$, one-dimensional boundary component counts $\#\mathcal{I}_2$, and flag orbit counts $\#\mathcal{I}_{12}$.

- **Key Sentinel Cases**:

  - Case 1: $(\#\mathcal{I}_1, \#\mathcal{I}_2, \#\mathcal{I}_{12}) = (5, 9, 18)$

  - Case 87: $(\#\mathcal{I}_1, \#\mathcal{I}_2, \#\mathcal{I}_{12}) = (528, 24242, 72199)$ (stable component-preserving endpoint)

### B. Unpolarized Enriques Boundary Strata and Stabilizers

- **File**: [`tests/fixtures/unpolarized_enriques.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/unpolarized_enriques.json)

- **Primary Sources**:

  - Dutour Sikirić & Hulek (arXiv:2302.01679), Sections 3.2 and 3.3.

  - V. V. Nikulin, *Surfaces of type K3 with finite automorphism group and Enriques surfaces*, J. Soviet Math.
    22 (1983).

  - A. I. Sterk, *Finiteness results for automorphy groups of 2-reflective lattices and Enriques surfaces*, Math.
    Ann.
    272 (1985), 237–264.

- **Content**:

  - Gram matrix of $N = U \oplus U(2) \oplus E_8(-2)$ in standard basis.

  - Primitive isotropic line representatives $I_{1,1} = \mathbb{Z}e_1, I_{1,2} = \mathbb{Z}e_3$.

  - Primitive isotropic plane representatives $I_{2,1} = \mathbb{Z}e_1 + \mathbb{Z}e_3, I_{2,2} = \mathbb{Z}(2e_1 + 2e_2 + w) + \mathbb{Z}e_3$.

  - Exact discriminant orthogonal group stabilizer image indices: $1, 527$ (lines), $527, 23715$ (planes).

  - Stable group counts: $528 = 1 + 527$, $24242 = 527 + 23715$, and flag count $72199$.

### C. K3 Polarized Modular Boundary Strata

- **File**: [`tests/fixtures/k3_modular_strata.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/k3_modular_strata.json)

- **Primary Sources**:

  - *Degree 2 K3 Boundary*: Gianni Scattone, *On the compactification of moduli spaces for algebraic K3 surfaces*, Memoirs of the AMS, Vol. 70, No. 372 (1987).

  - *Degree 4 K3 Boundary*: Kathleen Jones, *The boundary of the moduli space of degree 4 K3 surfaces*, Ph.D. thesis, University of Bath (2000).

  - *Degree $2d$ Cusp Formula*: Simon Attwell-Duval, *The boundary of moduli spaces of polarized K3 surfaces*, Ph.D. thesis, University of Bath (2021).

  - *Unimodular Transitivity*: Martin Eichler, *Quadratische Formen und orthogonale Gruppen*, Springer-Verlag (1952).

- **Content**:

  - Degree 2: Baily–Borel counts $(1, 4, 4)$ and four Type II root types $E_8 \oplus E_8 \oplus A_1, E_7 \oplus D_{10}, D_{16} \oplus A_1, A_{17}$.

  - Degree 4: 9 Type II components with generalized root types $A_{11} \oplus E_6, A_{15} \oplus 2A_1, A_{17}, 2D_8 \oplus A_1, D_{10} \oplus E_7, D_{12} \oplus D_5, D_{16} \oplus A_1, D_{17}, 2E_8 \oplus A_1$.

  - Attwell-Duval formula test cases for squarefree degrees ($d = 2, 3, 5, 6, 10, 14, 15, 30, 42, 70, 105, 210$).

### D. Independent Tits Buildings (Dawes)

- **File**: [`tests/fixtures/dawes_buildings.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/dawes_buildings.json)

- **Primary Source**: Matthew Dawes, *The geometry of the boundary of orthogonal modular varieties*, Ph.D. thesis, University of Bath (2020), and implementation `buildings.sage`.

- **Content**: Independent line orbit, plane orbit, and incidence multigraph computations for $2U \oplus A_2$, $2U \oplus \langle-6\rangle \oplus \langle-2\rangle$, and $U \oplus U(2) \oplus A_2$.

### E. Conway–Sloane: Spinor-Genus Pair and Root-Lattice Automorphism Orders

- **File**: [`tests/fixtures/conway_sloane_cases.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/conway_sloane_cases.json)

- **Source**: J. H. Conway and N. J. A. Sloane, *Sphere Packings, Lattices and Groups*, 3rd ed., Springer, 1999 [CS99].

- **Content**:

  - `spinor_genus_pair_determinant_minus_128`: the ternary forms (51a) and (51b) of Chapter 15 §11 ("Computational complexity"). The text states that both forms lie in the genus I_{2,1}(2 × 64), that this genus contains two spinor genera and hence two classes, and that (51b) represents the second class; Corollary 22 (n = 3, d_0 = 128) is the cited reason.
    Every field of the record transcribes those sentences.

  - `root_lattice_automorphism_orders`: the automorphism-group orders of A2, A3, D4, E6, E7, E8 from Chapter 4 (§4 for the convention g = g_0 g_1; §6 for A_n; §7 for D_n; §8 for E6, E7, E8), each with the quoted sentence.
    The Gram matrices are the Cartan matrices of the fundamental roots (Chapter 4, Table 4.1). The record `indefinite_jl_E8` pairs the Gram matrix shipped in Indefinite.jl's `TestCases/LATT_AUTOMORPHISM_case1_ListMat_E8` with the same order; the identification with E8 rests on Chapter 2, Table 2.2 (one even unimodular lattice in dimension 8).

- **Extraction**: transcribed from the Zotero extraction of the book (item T2WVLTDB). The `quote` fields are verbatim up to ASCII rendering of the mathematics.

* * *

## 2. Upstream Algorithmic Regression Corpora

### A. 8,821 Isotropic Decision Cases

- **File**: [`tests/fixtures/isotropic_cases_8821.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/isotropic_cases_8821.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/DATA/IsotropicCases`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: 8,821 exact integer Gram matrices with ground-truth boolean flags indicating the existence of a non-zero integral isotropic vector.

### B. 8,821 Reflective Lorentzian Forms

- **File**: [`tests/fixtures/reflective_forms_8821.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/reflective_forms_8821.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/20_Reflective/ListReflect`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: 8,821 Lorentzian Gram matrices with exact number of simple roots ($n_{\mathrm{simple}}$).

### C. 145 Lorentzian Equivalence Instances

- **File**: [`tests/fixtures/lorentzian_equivalence_146.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/lorentzian_equivalence_146.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/28B_LorentzianPerfStabEqui/TestCasesEqui.tar.gz`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: 145 pairs of equivalent Lorentzian Gram matrices $(M_1, M_2)$ with exact unimodular transporter matrices $U \in \mathrm{GL}_n(\mathbb{Z})$ satisfying $U M_1 U^{\mathsf{T}} = M_2$.

### D. Lorentzian Stabilizer Generators

- **File**: [`tests/fixtures/lorentzian_stabilizers_cases.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/lorentzian_stabilizers_cases.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/28B_LorentzianPerfStabEqui/TestCasesStab.tar.gz`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: 3 full Lorentzian orthogonal group generator sets.

### E. 103 Root Systems of Reflective Forms

- **File**: [`tests/fixtures/root_systems_56.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/root_systems_56.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/01_RatIntAutomorphy/ListSimpleRootSystem_4_56_X_5_47`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: 103 root systems of reflective forms in dimensions 4 and 5.

### F. Classified Simplices (Dimensions 5, 6, 7)

- **File**: [`tests/fixtures/classification_simplices.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/classification_simplices.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/DATA/ClassificationSimplices5`, `ClassificationSimplices6`, `ClassificationSimplices7`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: 12 classified simplices across dimensions 5, 6, and 7.

### G. 18 Finite Double-Coset Instances

- **File**: [`tests/fixtures/double_coset_cases.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/double_coset_cases.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/DoubleCosets/DBL/`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: 18 frozen finite group double coset instances across permutation degrees 55–96 with orders up to 92,160.

### H. 40 Lorentzian Perfect Domain Forms

- **File**: [`tests/fixtures/lorentzian_perfect_domains.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/lorentzian_perfect_domains.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/28B_LorentzianPerfStabEqui/Result_Enumeration`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: 40 exact Gram matrices across ranks 3, 4, 5, 6 with exact isotropic and total perfect domain orbit counts.

### I. 6 Metamorphic Indefinite Forms

- **File**: [`tests/fixtures/ci_indefinite_comp.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/ci_indefinite_comp.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/19_IndefiniteComp/AllTests.g`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: Forms `["U", "2U"]`, `["U", "2U", "A2"]`, `["U", "2U", "A3"]`, `["U", "2U", "A2", "A2"]`, `["U", "U", "E7"]`, `["U", "2U", "2E8"]`.
