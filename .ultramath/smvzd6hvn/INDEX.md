# UltraMath run `smvzd6hvn` — retained artifacts

Problem: an orbit-representative facet enumeration for the perfect-cell cones, replacing the complete
H-description that `cell_cone` asks an engine for. 8 members, 240 min, reasoning effort max. Two directions were
accepted; the rest ended without a submitted solution, but four left material worth keeping.

| member | direction | state | artifact worth reading |
| --- | --- | --- | --- |
| m6 | proof / contract layer | **accepted** | `m6/THEOREM_PACKAGE.md` — facet-normal canonicity, ridge normal cone is a 2-cone with the two adjacent facets as extreme rays, ridge-closure completeness, stabiliser = set stabiliser of the maximal incidence, transport predicate P1-P5 |
| m7 | exact certificate primitives | **accepted** | `m7/solution.md`, `m7/code/primitives.py`, `m7/measurements.md`, `m7/results/bench.json` — 499 facets over 5 cones against Sage's enumeration |
| m3 | the adjacency walk | unsubmitted | `m3/spec.md`, `m3/prototype/orbit_addm.py`, `m3/prototype/validate_specimens.py` — 12 specimens agreeing with an exact full dual description |
| m5 | independent ground truth | unsubmitted | `m5/orbit_facets.py`, `m5/run.log` — a `Fraction`-only double description with facet and orbit counts on 14 small cones |
| m2 | source-faithful `f_adj` reconstruction | unsubmitted | `m2/orbit_facets_proto.py`, `m2/report.md` |
| m4 | polar-LP separation with an exact completeness certificate | unsubmitted | `m4/proof.md`, `m4/proto/` |
| m8 | arithmetic candidates / streaming certification | unsubmitted | `m8/facet_orbit_spec.md`, `m8/check_lemmas.py` |

What the port uses, and where: the accepted contract layer and the certificate layer were ported to
`src/sage_indefinite_port/indefinite/facet_certificates.py`; the walk is
`src/sage_indefinite_port/indefinite/facet_orbits.py`, validated against `cell_cone` on `A2+<1>+<-1>` and
`A3+<1>+<-1>` and against the closed-form cube cones by `tests/test_facet_orbits.py`.

Two things this run established that must not be re-derived: the reduction cannot be done by feeding one ray per
orbit to an engine (`cone{(1,0),(2,0)}` under the swap), and the dual slice `{u : <u,r> >= 0, <h,u> = 1}` is
bounded only when `h` is an interior point of the **primal** cone -- a functional positive on the rays leaves it
unbounded, which is how the first seed attempt failed.
