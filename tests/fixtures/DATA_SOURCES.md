# Data Sources and Provenance

This document records the exact provenance, source files, authors, publications, and extraction methods for every dataset in [`tests/fixtures/`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures).

* * *

## 1. Published Mathematical Tables and Theorem Oracles

### A. 87 Numerical Enriques Polarizations

- **File**: `tests/fixtures/enriques_87_polarizations.json`

- **Primary Source**: Mathieu Dutour Sikirić and Klaus Hulek, *Moduli of polarised Enriques surfaces — computational aspects*, J. London Math. Soc. 2024, [arXiv:2302.01679](https://arxiv.org/abs/2302.01679), Tables `table_subgroups1` and `table_subgroups2`; vendored at `references/vendor/arxiv/2302.01679/Enriques_compu_rev.tex`.

- **Extraction**: `uv run references/extract/dh23_polarizations.py` parses both tables from the vendored TeX. Each record carries the TeX line it came from.

- **Content**: for each of the 87 conjugacy classes of groups $\Gamma_h$: the subset $S$ of the $E_{10}$ diagram, $\#S$, $|\bar\Gamma_h|$, the isotropic line, plane and flag counts $\#I_1, \#I_2, \#I_{12}$, the degree of the smallest realization $h_{\min}$, and $\phi(h_{\min})$.

- **Check**: the script asserts cases 1–87 in order and that case 87 has 528 lines and 24242 planes, which the paper also states independently in its text (lines 1168–1170).

### B. Unpolarized Enriques Boundary Strata and Stabilizers

- **File**: `tests/fixtures/unpolarized_enriques.json`

- **Primary Source**: Dutour Sikirić & Hulek, arXiv:2302.01679, subsection "The Tits building" (vendored TeX `references/vendor/arxiv/2302.01679/Enriques_compu_rev.tex`, lines 1159–1170), and Table `table_subgroups2` case 87 (line 814) for the flag count.

- **Extraction**: `uv run references/extract/dh23_unpolarized.py` builds $N = U \oplus U(2) \oplus E_8(-2)$, takes the representatives from the text, and asserts $\det N = 1024$ and that every line and plane representative is totally isotropic.

- **Content**:

  - Gram matrix of $N$, with $(e_1, e_2)$ and $(e_3, e_4)$ the standard bases of $U$ and $U(2)$.

  - Isotropic line representatives $\mathbb{Z}e_1$, $\mathbb{Z}e_3$, with discriminant-image stabilizer indices 1 and 527.

  - Isotropic plane representatives $\mathbb{Z}e_1 + \mathbb{Z}e_3$ and $\mathbb{Z}(2e_1 + 2e_2 + w) + \mathbb{Z}e_3$, with indices 527 and 23715. Here $w = \alpha_0 + \alpha_2$ has norm 4 in $E_8$, i.e. $w^2 = -8$ in $E_8(-2)$, as the paper requires.

  - The stable group $\widetilde O^+(N)$ has 528 lines and 24242 planes (text) and 72199 flags (table).

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
