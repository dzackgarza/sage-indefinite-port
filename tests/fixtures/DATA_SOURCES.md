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

- **File**: [`tests/fixtures/classification_simplices.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/classification_simplices.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/DATA/ClassificationSimplices5`, `ClassificationSimplices6`, `ClassificationSimplices7`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: 12 classified simplices across dimensions 5, 6, and 7.

### 18 Finite Double-Coset Instances

- **File**: [`tests/fixtures/double_coset_cases.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/double_coset_cases.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/DoubleCosets/DBL/`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: 18 frozen finite group double coset instances across permutation degrees 55–96 with orders up to 92,160.

### 40 Lorentzian Perfect Domain Forms

- **File**: `tests/fixtures/lorentzian_perfect_domains.json`

- **Upstream**: `references/vendor/polyhedral_common@1592b246/CI_tests/28B_LorentzianPerfStabEqui/Result_Enumeration`: for 40 lattices of ranks 3–6, the exact number of perfect-domain orbits in the "isotropic" and "total" modes.

- These are the reference implementation's own recorded outputs (polyhedral_common, Mathieu Dutour Sikirić); agreement with them is acceptance for this port. Extracted by `references/extract/polyhedral_common_ci.py` (run under Sage's Python; GAP literals are evaluated with libgap). Each record cites its vendored file and index.

### 6 Metamorphic Indefinite Forms

- **File**: [`tests/fixtures/ci_indefinite_comp.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/ci_indefinite_comp.json)

- **Upstream Location**: `references/polyhedral_common/CI_tests/19_IndefiniteComp/AllTests.g`

- **Upstream Author**: Mathieu Dutour Sikirić

- **Repository**: [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

- **Content**: Forms `["U", "2U"]`, `["U", "2U", "A2"]`, `["U", "2U", "A3"]`, `["U", "2U", "A2", "A2"]`, `["U", "U", "E7"]`, `["U", "2U", "2E8"]`.
