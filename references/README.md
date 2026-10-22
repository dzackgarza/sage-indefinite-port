# Pinned algorithm and literature sources

The authoritative reference implementations for this port are **committed**
under `references/vendor/`, not local machine-only reference checkouts:

- [`polyhedral_common@1592b246`](vendor/polyhedral_common@1592b246/) — pinned
  C++ Lorentzian and indefinite algorithms (read for translation; do not compile).
- [`Indefinite.jl@374a5ebb`](vendor/Indefinite.jl@374a5ebb/) —
  pinned Julia/GAP routines (including `Lorentzian.g`).
- [Allcock reference and rank-three cases](vendor/arxiv/2110.03784/) and
  [rank-three lattice table](vendor/GeometryDatabase_Rank3_Lorentzian_lattices@63a9067a/).
- [Complete provenance and revision pins](vendor/README.md).

Use the [source-to-Sage translation map](source-to-sage-translation-map.md)
and [port-design assignments](port-design.md) to locate the required
routine and its exact generic callee. Some historically used directories
`references/polyhedral_common/` and `references/Indefinite.jl/` were
ignored local clones; they are **not** prerequisites for a cloud checkout.

## Mathematical Literature

For complete academic citations, bibliography, and BibTeX entries, see [CITATIONS.md](../CITATIONS.md).

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

- **Brandhorst (2021)**: *The classification of reflective hyperbolic lattices of rank $\ge 4$*, Math.
  Comp.
