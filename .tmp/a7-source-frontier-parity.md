# Source adjacency decomposition: pending facets and termination

The original A7 assertion at bb69c92/aa75796e is terminal under its original bound. The isolated walk on the actual 49,685-row coupled cell is also terminal under its diagnostic bound. The retained slow native comparison is positive: its witness maps one supplied facet incidence to the other. Point-orbit occupancy agrees, so that invariant cannot reject this comparison.

## Concrete source differences

The port's recursive walk selects `pending.pop()` and continues until that worklist is empty. The pinned source selects an orbit of minimum facet incidence and uses an exact Balinski completion certificate:

- `Indefinite.jl@374a5ebb/indef/lib/SetFunctionality.g:304-323` selects minimum incidence. Its `OrbitIntersection` and `ComputeIntersectionUndone` give the intersection across all translates of pending facets.
- `DualDescriptionAndADM.g:133-145` stops after a completed orbit if the remaining known facet count is at most rank(EXT)-2, or if all pending facets across their full group orbits share a ray.
- Pinned [POLY_RecursiveDualDesc.h](https://github.com/MathieuDutSik/polyhedral_common/blob/1592b2463bc7a04d94e900bf8efbada607405d8d/src_dualdesc/POLY_RecursiveDualDesc.h):837-855,1392-1401 likewise selects minimum incidence and checks termination before expansion.
- Pinned [Balinski_basic.h](https://github.com/MathieuDutSik/polyhedral_common/blob/1592b2463bc7a04d94e900bf8efbada607405d8d/src_dualdesc/Balinski_basic.h):55-60 requires a processed orbit and the small-frontier or common-ray certificate.

These are source stopping and scheduling rules, not a bounded search substitution. Every already discovered record is returned, including those whose further expansion is certified unnecessary.

## Certificate

Intersecting all translates of a set under a finite group leaves exactly the full point orbits contained in that set. Thus the source's common-ray test is equivalent to asking whether the intersection of the pending representative incidences contains a whole point orbit. The group owner supplies these point orbits.

A pointed full-dimensional rank-d cone has a (d-1)-dimensional polytope section. The graph of its dual is (d-1)-vertex-connected. A nonempty processed region cannot be separated from an undiscovered region by at most d-2 pending facets. If the full pending family shares a nonzero ray, those dual vertices lie in a proper face. The graph outside a proper face is connected: choose a generic increasing linear objective whose value on every vertex outside the face exceeds its value on that face, and follow increasing edges to its unique maximum. Every vertex on a proper face also has an edge leaving that face. Complete processed adjacency and these facts exclude undiscovered vertices in either region.

Orbit sizes count all facets in each disjoint orbit. The shared-ray test uses the entire orbit family, not only its representatives. The completed-orbit condition applies after all ridge orbits of the selected facet have been expanded. Existing exact rank, normal, incidence, integral action and transport-word checks stay in force.

## Mathematical checks

The coordinate facets of Delta(8,2) form complete families of incidence cardinalities seven and twenty-one. Full coordinate permutations give two orbits; fixing one coordinate separates its two facets from the moving families. The assertion checks every family, its exact orbit cardinality and every retained transport word in both word modes. All four cases pass in 7.53 seconds. The six signed-cube and cube-coordinate cases and the complete signed six-cube transport assertion also pass (nine selected cases with the original hypersimplex pair, 9.95 seconds).

The source `OrbitIntersection` is evaluated through the finite group owner's point orbits: the intersection of all translated pending sets is the union of the point orbits contained in their ordinary intersection. This is an exact callee substitution with the same full quantifier. The record and witness interfaces are unchanged. The large actual A7 facet computation and the original public consumer remain required.

## Actual large configuration

The complete facet operation on the retained large configuration finishes at published port 616a7fe/provider cfee1892. The input is the actual stable ray family and complete generator permutations from the earlier coupled cell. Each discovered facet retains its exact normal/incidence certificate and full orbit cardinality.

| Complete source operation | Value |
| --- | ---: |
| Input rays | 49685 |
| Facet orbits | 19 |
| Total facets from the disjoint orbit sizes | 17611202 |
| Wall seconds | 284.201973 |
| CPU seconds | 287.707900 |

This measures the finite facet core supplied with its actual stabilizer action. The earlier bounded isolated whole-cell walk also included stabilizer construction and diagnostic capture, so its time is not a like-for-like benchmark ratio. The original full A7 assertion at 616a7fe/cfee1892 is terminal at its unchanged original 1800-second thread bound. Its terminal stack is in exact facet certification and integer row conversion. The finite core result does not establish the whole orthogonal group or its discriminant image.

Optional data: `.tmp/a7-bb69c92-facet-comparison/base-0.pickle`, `.tmp/a7-frontier-source-large-cell.json` and `.tmp/a7-frontier-source-large-cell-records.pickle`.

## Retained matrix certificate measurement

The complete retained facet representative family on the same ray configuration was checked with the existing Python row arithmetic and with one native integer ray matrix retained across all certificates. Every certificate is true in both executions. Matrix construction is charged separately.

| Operation | Seconds |
| --- | ---: |
| One native matrix construction | 0.299812 |
| Existing complete certificate family | 3.158803 |
| Retained native matrix certificate family | 4.523677 |

The retained matrix substitution is slower and is rejected. Production certificate arithmetic stays unchanged. Sage's matrix constructor also does not accept `copy=False`; the measurement passes the retained matrix directly. Optional reproduction: `.tmp/a7-retained-matrix-certificate-benchmark.py`, `.json` and `.log`.

## Exact certificate row retention

The existing certificate converted all cone entries to integers twice: once when rebuilding the full row family, then again in its exact pairing expression. Retaining the sequence without the first conversion preserves all pairings and zero sets; the incidence rank operation still converts its own selected rows. An input iterable is materialized once, so both scans use the same rays.

| Complete retained family | Seconds |
| --- | ---: |
| Original certificate | 5.236482 |
| Retained row sequence | 3.282984 |

Every certificate agrees on the complete actual facet representative family. The two executions use the same process and inputs; other independent local work was running, so this wall measurement is a scoped comparison. The analytic cube certificate, engine facet recovery and complete coordinate-hypersimplex orbit/word assertions pass together in 9.67 seconds. No assertion about whole A7 completion follows from this local improvement. Optional paired input script/results: `.tmp/a7-retained-sequence-certificate-benchmark.{py,json,log}`.
