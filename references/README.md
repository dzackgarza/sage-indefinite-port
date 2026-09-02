# Reference Subtrees and Software Sources

This directory contains upstream reference repositories for algorithm extraction and differential verification:

- `polyhedral_common`: Clone of [MathieuDutSik/polyhedral_common](https://github.com/MathieuDutSik/polyhedral_common). C++ implementation of polyhedral computations, Delaunay polytopes, and indefinite quadratic form algorithms.

- `Indefinite.jl`: Clone of [MathieuDutSik/Indefinite.jl](https://github.com/MathieuDutSik/Indefinite.jl).
  Julia package for indefinite lattices and forms with vendored GAP algorithms.

## Pinned revisions

The mirrors are ignored by git; each checkout is pinned to the revision the translation plan was written against.

| Mirror | Revision | Notes |
| --- | --- | --- |
| `polyhedral_common` | `a55fcb7b71af48c88d7abbeab7889e9347916e43` | corpus layout: `CI_tests/01_RatIntAutomorphy`, `CI_tests/28B_LorentzianPerfStabEqui`, `CI_tests/DoubleCosets/DBL`, `CI_tests/19_IndefiniteComp` |
| `Indefinite.jl` | `374a5ebb5ae33e5688fd95052a4c68cbc91ae6aa` | GAP sources in `indef/lib/` |

## Mathematical Literature

For complete academic citations, bibliography, and BibTeX entries, see [CITATIONS.md](file:///home/dzack/gitclones/sage-indefinite-port/CITATIONS.md).

Key papers:

- **Dutour Sikirić & Hulek (2023)**: *On the classification of numerical Enriques surfaces and their moduli spaces*, [arXiv:2302.01679](https://arxiv.org/abs/2302.01679).

- **Dawes (2020)**: *The geometry of the boundary of orthogonal modular varieties*, Ph.D. thesis, University of Bath (`buildings.sage`).

- **Scattone (1987)**: *On the compactification of moduli spaces for algebraic K3 surfaces*, Memoirs of the AMS, Vol. 70, No. 372.

- **Jones (2000)**: *The boundary of the moduli space of degree 4 K3 surfaces*, Ph.D. thesis, University of Bath.

- **Attwell-Duval (2021)**: *The boundary of moduli spaces of polarized K3 surfaces*, Ph.D. thesis, University of Bath.

- **Nikulin (1980)**: *Integral symmetric bilinear forms and some of their applications*, Math.
  USSR Izv.
  14, 103–167.

- **Eichler (1952)**: *Quadratische Formen und orthogonale Gruppen*, Springer-Verlag.

- **Sterk (1985)**: *Finiteness results for automorphy groups of 2-reflective lattices and Enriques surfaces*, Math.
  Ann.
  272, 237–264.

- **Conway & Sloane (1999)**: *Sphere Packings, Lattices and Groups (SPLAG)*, 3rd ed., Springer.

- **Brandhorst (2021)**: *The classification of reflective hyperbolic lattices of rank $\ge 4$*, Math.
  Comp.
