# Exact facet representatives for orthogonal-group assembly

Port computation: `96c762f03ebc8cd4276a1741e7e85c774b43a63c`.
Provider: `01cfd2db26e5bbd4c794819781d80d3417575dfe`.

Pinned `Indefinite.jl@374a5ebb/indef/lib/SetFunctionality.g:174–201,234–260` compares facet incidences using `RepresentativeAction(G, F, F', OnSets)` and computes orbit sizes as `|G|/|Stab_G(F)|`. `Lorentzian.g:1286–1308` supplies that formalism to adjacency decomposition. These are exact group operations; an invariant only filters candidates before the full equivalence decision.

Orthogonal-group assembly requests one representative per facet orbit. That request now reaches every recursive facet computation. The exact GAP set transporter decides whether a candidate belongs to an existing orbit, and the set stabilizer gives its full size. The source walk still processes every ridge orbit and certifies every accepted facet. Callers requesting full facet transports retain every member and generator word through the existing enumeration.

The independent cube assertions use the standard coordinate/sign action. Signed coordinate permutations act transitively on the 2d facets. Coordinate permutations without sign changes have two orbits, each of size d. The assertions certify the returned normals and distinguish these two orbits for d=2,3,6, exercising direct and recursive enumeration. Existing A2/A3 complete facet partitions and generator-word transports, the six-cube full-word action and the null-quadrant swap also pass.

| Selected assertions | Result |
| --- | --- |
| Cube representative and full-word assertions | 10 pass |
| Selected session total seconds | 22.56 |

A read-only capture from the unchanged A7 consumer at bf43db1/01cfd2db supplies a real finite action for a cell with 49,861 rays. A one-shot Python monitoring callback copied only its incidence and generator permutations, then released its monitoring tool. A subsequent state inspection confirms that the tool is free. The callback uses Python's documented per-code PY_START event; the running mathematical operations are unchanged.

| Captured finite action | Value |
| --- | ---: |
| Incidence cardinality | 20 |
| Supplied generator permutations | 11 |
| Full permutation-group order | 185794560 |
| Exact set-stabilizer order | 5760 |
| Orbit size from orbit-stabilizer | 32256 |
| Independently enumerated orbit size | 32256 |
| Set-stabilizer computation seconds | 2.092 |
| Set transporter and exact image certificate seconds | 2.397 |
| Full orbit enumeration seconds | 4.805 |

The selected set transporter is independently checked by applying it to every incident position and comparing the resulting set. These measurements concern this captured orbit only. They do not establish the original A7 discriminant-image assertion. The immutable bf43db1/01cfd2db original proof remains active under its 1800-second bound. Preserve it until terminal, then exercise the representative-only source revision on the unchanged consumer.
