# Dawes boundary for 2U(2) + A2

The named A2 is negative definite, with Gram [[-2,1],[1,-2]]. Let e1,f1 and e2,f2 be the canonical bases of the two U(2) summands. Let Gamma be the discriminant kernel intersected with O_plus. Every element of Gamma fixes the discriminant module pointwise.

## Independent integral witnesses

The following primitive vectors are null and have divisor two:

- e1, f1, e2, f2;
- e1+e2, e1+f2, f1+e2, f1+f2;
- e1+f1+e2-f2.

Their divided discriminant classes are distinct and have order two. Negation fixes each such class. Thus their primitive unoriented lines belong to distinct Gamma orbits. This proves a lower bound for the line vertices without invoking the port's orbit classifier.

The following primitive null planes give every maximal isotropic plane of the hyperbolic 2-primary discriminant quadric:

| Plane basis | Nonzero divided classes |
| --- | --- |
| e1, e2 | [e1/2], [e2/2], [(e1+e2)/2] |
| e1, f2 | [e1/2], [f2/2], [(e1+f2)/2] |
| f1, e2 | [f1/2], [e2/2], [(f1+e2)/2] |
| f1, f2 | [f1/2], [f2/2], [(f1+f2)/2] |
| e1+e2, f1-f2 | [(e1+e2)/2], [(f1+f2)/2], [(e1+f1+e2+f2)/2] |
| e1+f2, f1-e2 | [(e1+f2)/2], [(f1+e2)/2], [(e1+f1+e2+f2)/2] |

Each plane is primitive: its displayed basis has an integral minor of determinant one. Distinct discriminant planes cannot belong to one Gamma orbit. Each nonzero point belongs to exactly two of these planes. Their point-in-plane pairs also give distinct Gamma flag classes.

| Independent discriminant lower bound | Value |
| --- | ---: |
| Line vertices | 9 |
| Plane vertices | 6 |
| Flag incidences | 18 |

These are lower bounds for the integral arithmetic quotient. The published diagram supplies the complete quotient and the same incidence graph.

## Pinned diagram and frozen oracle

`references/vendor/arxiv/2205.10601/orbits_in_lattices.tex:803-809` defines black vertices as rational isotropic lines and white vertices as planes. The diagram labelled `2u2a2building`, at lines 1450-1495, has the following vertex colours:

| Pinned diagram | Value |
| --- | ---: |
| Black line vertices | 9 |
| White plane vertices | 6 |
| Edges | 18 |

The frozen `dawes_2U2_A2_stable` fixture assigns the opposite vertex types. Its stored edge labels retain the diagram, with the two types exchanged. The independent assertion `tests/test_dawes_scaled_planes_discriminant_flags.py` verifies every null, primitive and divisor equation, all the distinct divided classes, the complete finite discriminant-plane list and the coloured graph isomorphism after transposing those types. It passes at provider aa75796e (13.79 seconds with the normal pytest gate).

The proposed correction is `references/dawes-2U2-A2-oracle-correction.patch`. It exchanges the fixture's point and curve vertex sets and transposes every incidence pair. It also corrects the building's point and curve cardinalities. The patch applies cleanly to the current frozen file. Applying it requires the decision on issue #33 under the unchanged-suite requirement.

## Required implementation result

The captured old-pair failure has line e2 and returned representatives e1, f2, e1+e2, f1+f2 and f1. All primitive-isotropic universe checks pass. Their divided discriminant classes omit that of e2. The classifier must return a complete family, with exact subgroup transporters and the full incidence graph. The independent lower bound and diagram correction do not establish that implementation result.
