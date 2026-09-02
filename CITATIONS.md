# Citations and Mathematical References

This repository implements algorithms for indefinite quadratic lattices, arithmetic orthogonal groups, Lorentzian perfect domains, and Baily–Borel boundary strata.
The implementation, oracle corpus, and acceptance test fixtures are grounded in the following published literature and software references.

* * *

## 1. Primary Mathematical Literature

### Orthogonal Groups, Cusps, and Moduli Spaces of Enriques & K3 Surfaces

- **[DH23]** Mathieu Dutour Sikirić and Klaus Hulek.
  *On the classification of numerical Enriques surfaces and their moduli spaces*. [arXiv:2302.01679](https://arxiv.org/abs/2302.01679) [math.AG], 2023.

  - *Contribution*: Full classification of 87 conjugacy classes of arithmetic subgroups $\Gamma \leq O^\Omega(N)$ for the Enriques anti-invariant lattice $N = U \oplus U(2) \oplus E_8(-2)$.

  - *Oracles provided*: Stabilizer image indices $(1, 527, 527, 23715)$, stable building counts $(528, 24242, 72199)$, isotropic line/plane orbit representatives, and complete 87-case $(\#\mathcal{I}_1, \#\mathcal{I}_2, \#\mathcal{I}_{12})$ trie tables.

- **[Daw20]** Matthew Dawes.
  *The geometry of the boundary of orthogonal modular varieties*. Doctoral thesis, University of Bath, 2020. Accompanying software: `buildings.sage`.

  - *Contribution*: Independent SageMath algorithm for Baily–Borel boundary components and Tits buildings for maximal lattices and congruence subgroups of $O(2,n)$.

  - *Oracles provided*: Independent Tits building multigraphs and orbit counts for $2U \oplus A_2$, $2U \oplus \langle-6\rangle \oplus \langle-2\rangle$, $U \oplus U(2) \oplus A_2$, and $U \oplus U(m) \oplus A_2$.

- **[Sca87]** Gianni Scattone.
  *On the compactification of moduli spaces for algebraic K3 surfaces*. Memoirs of the American Mathematical Society, Vol. 70, No. 372, 1987.

  - *Contribution*: Geometric identification of zero- and one-dimensional boundary strata for degree-two polarized K3 surfaces.

  - *Oracles provided*: Exact $(1, 4, 4)$ Baily–Borel counts and the four Type II boundary component Cartan root types: $E_8 \oplus E_8 \oplus A_1$, $E_7 \oplus D_{10}$, $D_{16} \oplus A_1$, and $A_{17}$.

- **[Jon00]** Kathleen Jones.
  *The boundary of the moduli space of degree 4 K3 surfaces*. Doctoral thesis / J. Algebraic Geom., 2000.

  - *Contribution*: Classification of the 9 Type II boundary components for the degree-four K3 modular variety.

  - *Oracles provided*: Root types $A_{11} \oplus E_6$, $A_{15} \oplus 2A_1$, $A_{17}$, $2D_8 \oplus A_1$, $D_{10} \oplus E_7$, $D_{12} \oplus D_5$, $D_{16} \oplus A_1$, $D_{17}$, and $2E_8 \oplus A_1$.

- **[AD21]** Simon Attwell-Duval.
  *The boundary of moduli spaces of polarized K3 surfaces*. Thesis / arXiv, 2021.

  - *Contribution*: Exact closed-form lower bound for zero-dimensional cusps of degree $2d$ polarized K3 surfaces.

  - *Oracles provided*: For squarefree $d > 1$ with $\omega(d)$ distinct prime factors, the number of zero-dimensional cusps is exactly $2^{\omega(d)-1}$.

* * *

## 2. Foundational Lattice Theory and Discriminant Forms

- **[Nik80]** Vyacheslav V. Nikulin.
  *Integral symmetric bilinear forms and some of their applications*. Mathematics of the USSR-Izvestiya, 14(1):103–167, 1980.

  - *Contribution*: Theory of finite quadratic forms on discriminant modules $A_L = L^\#/L$, genus theory, Nikulin glue maps $H_S \to H_R(-1)$, and primitive embeddability criteria into even unimodular lattices $\mathrm{II}_{p,q}$.

- **[Eic52]** Martin Eichler.
  *Quadratische Formen und orthogonale Gruppen*. Grundlehren der mathematischen Wissenschaften, Band 63, Springer-Verlag, Berlin-Göttingen-Heidelberg, 1952.

  - *Contribution*: Eichler criterion for transitivity of stable orthogonal groups $O^+(L)$ on primitive vectors of fixed norm in even lattices containing $2U$.

- **[Ste85]** A. I. Sterk.
  *Finiteness results for automorphy groups of 2-reflective lattices and Enriques surfaces*. Mathematische Annalen, 272(2):237–264, 1985.

  - *Contribution*: Finiteness of fundamental polyhedra and Coxeter/reflective group properties of Enriques lattices.

- **[CS99]** John H. Conway and Neil J. A. Sloane.
  *Sphere Packings, Lattices and Groups (SPLAG)*. Grundlehren der mathematischen Wissenschaften, Vol. 290, Springer-Verlag, New York, 3rd edition, 1999.

  - *Contribution*: Reference classifications for root lattices ($A_n, D_n, E_6, E_7, E_8$), Leech lattice $\Lambda_{24}$, $p$-adic Jordan decompositions, and spinor genera (Table 15.1).

- **[Bra21]** Simon Brandhorst.
  *The classification of reflective hyperbolic lattices of rank $\geq 4$*. Mathematics of Computation, 2021.

  - *Contribution*: Classification and symmetry groups of reflective hyperbolic and Lorentzian lattices.

* * *

## 3. Upstream Software and Algorithmic Sources

- **[Dut-PC]** Mathieu Dutour Sikirić.
  `polyhedral_common`: C++ library for polyhedral computations, Delaunay polytopes, and indefinite quadratic forms.

  - Repository: [https://github.com/MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common)

  - *Extracted components*: `INDEF_FORM_*` algorithms, Lorentzian perfect domain traversal (`TestPerfLorentzian.g`, `Result_Enumeration`), and rational matrix group integralization (`01_RatIntAutomorphy`).

- **[Dut-JL]** Mathieu Dutour Sikirić.
  `Indefinite.jl`: Julia package for indefinite lattices and forms.

  - Repository: [https://github.com/MathieuDutSik/Indefinite.jl](https://github.com/MathieuDutSik/Indefinite.jl)

  - *Extracted components*: Definite leaf validation fixtures ($E_8, 2E_8$), Lorentzian matrix test pairs ($U\_2U\_2I3, U\_E8, U\_I3$), and vendored GAP algorithms.

- **[GAP]** The GAP Group.
  *GAP — Groups, Algorithms, and Programming*, Version 4.13, 2024.

  - URL: [https://www.gap-system.org](https://www.gap-system.org)

  - *Role*: libGAP integration for permutation group actions, double coset enumeration, and finite quotient stabilizer computations.

- **[FLINT]** William Hart, Fredrik Johansson, and Sebastian Pancratz.
  *FLINT: Fast Library for Number Theory*, 2024.

  - URL: [https://flintlib.org](https://flintlib.org)

  - *Role*: Exact integer linear algebra, Hermite Normal Form (HNF), and Smith Normal Form (SNF) saturation.

- **[Normaliz]** Winfried Bruns, Bogdan Ichim, and Christof Söger.
  *Normaliz: Algorithms for rational cones and affine monoids*, 2024.

  - URL: [https://www.normaliz.uni-osnabrueck.de](https://www.normaliz.uni-osnabrueck.de)

  - *Role*: Exact rational polyhedral cone facet enumeration and polyhedral reduction backends.

- **[Bliss]** Tommi Junttila and Petteri Kaski.
  *Bliss: A Tool for Computing Automorphism Groups and Canonical Labelings of Graphs*, 2015.

  - URL: [http://www.tcs.hut.fi/Software/bliss/](http://www.tcs.hut.fi/Software/bliss/)

  - *Role*: Canonical graph labeling for pairing configurations and lattice isometry transporter lifting.

* * *

## 4. BibTeX Database

```bibtex
@article{DutourSikiricHulek2023,
  author    = {Mathieu Dutour Sikiri{\'{c}} and Klaus Hulek},
  title     = {On the classification of numerical {E}nriques surfaces and their moduli spaces},
  journal   = {arXiv preprint arXiv:2302.01679},
  year      = {2023},
  eprint    = {2302.01679},
  archivePrefix = {arXiv},
  primaryClass  = {math.AG}
}

@phdthesis{Dawes2020,
  author    = {Matthew Dawes},
  title     = {The geometry of the boundary of orthogonal modular varieties},
  school    = {University of Bath},
  year      = {2020}
}

@article{Scattone1987,
  author    = {Gianni Scattone},
  title     = {On the compactification of moduli spaces for algebraic {K}3 surfaces},
  journal   = {Memoirs of the American Mathematical Society},
  volume    = {70},
  number    = {372},
  year      = {1987},
  publisher = {American Mathematical Society}
}

@phdthesis{Jones2000,
  author    = {Kathleen Jones},
  title     = {The boundary of the moduli space of degree 4 {K}3 surfaces},
  school    = {University of Bath},
  year      = {2000}
}

@phdthesis{AttwellDuval2021,
  author    = {Simon Attwell-Duval},
  title     = {The geometry of the boundary of moduli spaces of polarized {K}3 surfaces},
  school    = {University of Bath},
  year      = {2021}
}

@article{Nikulin1980,
  author    = {Vyacheslav V. Nikulin},
  title     = {Integral symmetric bilinear forms and some of their applications},
  journal   = {Mathematics of the USSR-Izvestiya},
  volume    = {14},
  number    = {1},
  pages     = {103--167},
  year      = {1980}
}

@book{Eichler1952,
  author    = {Martin Eichler},
  title     = {Quadratische {F}ormen und orthogonale {G}ruppen},
  series    = {Die Grundlehren der mathematischen Wissenschaften},
  volume    = {63},
  publisher = {Springer-Verlag},
  address   = {Berlin-G{\"o}ttingen-Heidelberg},
  year      = {1952}
}

@article{Sterk1985,
  author    = {A. I. Sterk},
  title     = {Finiteness results for automorphy groups of 2-reflective lattices and {E}nriques surfaces},
  journal   = {Mathematische Annalen},
  volume    = {272},
  number    = {2},
  pages     = {237--264},
  year      = {1985}
}

@book{SPLAG1999,
  author    = {John H. Conway and Neil J. A. Sloane},
  title     = {Sphere Packings, Lattices and Groups},
  series    = {Grundlehren der mathematischen Wissenschaften},
  volume    = {290},
  edition   = {3rd},
  publisher = {Springer-Verlag},
  address   = {New York},
  year      = {1999}
}

@article{Brandhorst2021,
  author    = {Simon Brandhorst},
  title     = {The classification of reflective hyperbolic lattices of rank $\ge 4$},
  journal   = {Mathematics of Computation},
  year      = {2021}
}

@misc{polyhedral_common,
  author       = {Mathieu Dutour Sikiri{\'{c}}},
  title        = {polyhedral\_common: {C}++ library for polyhedral computations and indefinite quadratic forms},
  howpublished = {\url{https://github.com/MathieuDutSik/polyhedral_common}},
  year         = {2024}
}

@misc{IndefiniteJL,
  author       = {Mathieu Dutour Sikiri{\'{c}}},
  title        = {Indefinite.jl: {J}ulia package for indefinite lattices and forms},
  howpublished = {\url{https://github.com/MathieuDutSik/Indefinite.jl}},
  year         = {2024}
}
```

* * *

## 5. Concrete Fixture Datasets and Source Mapping

For detailed per-record extraction protocols, see [DATA_SOURCES.md](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/DATA_SOURCES.md).

| Fixture File | Concrete Cases | Primary Data Source | Literature Citation / Upstream Origin |
| :--- | :--- | :--- | :--- |
| [`tests/fixtures/isotropic_cases_8821.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/isotropic_cases_8821.json) | 8,821 | `CI_tests/DATA/IsotropicCases` | Mathieu Dutour Sikirić (`polyhedral_common`) |
| [`tests/fixtures/reflective_forms_8821.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/reflective_forms_8821.json) | 8,821 | `CI_tests/20_Reflective/ListReflect` | Mathieu Dutour Sikirić (`polyhedral_common`) |
| [`tests/fixtures/lorentzian_equivalence_146.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/lorentzian_equivalence_146.json) | 145 | `CI_tests/28B_LorentzianPerfStabEqui/TestCasesEqui.tar.gz` | Mathieu Dutour Sikirić (`polyhedral_common`) |
| [`tests/fixtures/root_systems_56.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/root_systems_56.json) | 103 | `CI_tests/01_RatIntAutomorphy/ListSimpleRootSystem_4_56_X_5_47` | Mathieu Dutour Sikirić (`polyhedral_common`) |
| [`tests/fixtures/enriques_87_polarizations.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/enriques_87_polarizations.json) | 87 | Tables 1 & 2 | Dutour Sikirić & Hulek (2023, arXiv:2302.01679) |
| [`tests/fixtures/lorentzian_perfect_domains.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/lorentzian_perfect_domains.json) | 40 | `CI_tests/28B_LorentzianPerfStabEqui/Result_Enumeration` | Mathieu Dutour Sikirić (`polyhedral_common`) |
| [`tests/fixtures/double_coset_cases.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/double_coset_cases.json) | 18 | `CI_tests/DoubleCosets/DBL/` | Mathieu Dutour Sikirić (`polyhedral_common`) |
| [`tests/fixtures/k3_modular_strata.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/k3_modular_strata.json) | 14 | Boundary classifications & Attwell-Duval formula | Scattone (1987), Jones (2000), Attwell-Duval (2021) |
| [`tests/fixtures/classification_simplices.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/classification_simplices.json) | 12 | `CI_tests/DATA/ClassificationSimplices5/6/7` | Mathieu Dutour Sikirić (`polyhedral_common`) |
| [`tests/fixtures/indefinite_jl_cases.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/indefinite_jl_cases.json) | 6 | `Indefinite.jl/TestLor` & `TestCases` | Mathieu Dutour Sikirić (`Indefinite.jl`) |
| [`tests/fixtures/ci_indefinite_comp.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/ci_indefinite_comp.json) | 6 | `CI_tests/19_IndefiniteComp/AllTests.g` | Mathieu Dutour Sikirić (`polyhedral_common`) |
| [`tests/fixtures/unpolarized_enriques.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/unpolarized_enriques.json) | 4 | Section 3.2, 3.3 | Dutour Sikirić & Hulek (2023), Nikulin (1983), Sterk (1985) |
| [`tests/fixtures/lorentzian_stabilizers_cases.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/lorentzian_stabilizers_cases.json) | 3 | `CI_tests/28B_LorentzianPerfStabEqui/TestCasesStab.tar.gz` | Mathieu Dutour Sikirić (`polyhedral_common`) |
| [`tests/fixtures/dawes_buildings.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/dawes_buildings.json) | 3 | `buildings.sage` | Matthew Dawes (Ph.D. thesis 2020) |
| [`tests/fixtures/centralizer_involutions.json`](file:///home/dzack/gitclones/sage-indefinite-port/tests/fixtures/centralizer_involutions.json) | 4 | Structural gluing & Enriques involution | Nikulin (1980), Sterk (1985) |
