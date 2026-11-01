"""Baseline of the libGAP facet-orbit computation, for a future substitution test.

Constructs `_CellFaces` on a small frozen fixture cell and records what the current
libGAP path produces: facet count, incidence classes, the cell's permutation-group order,
and the facet orbits under the full group. A characterisation test for replacing these
GAP calls with the provider's maintained operations must reproduce exactly these numbers.
"""
import json
import time

from sage_extensions.all import Lattices, ZZ

from sage_indefinite_port.indefinite.lorentzian_cells import (
    NONNEGATIVE_SQUARES,
    LorentzianPerfectLocalBackend,
    _CellFaces,
)
from sage.libs.gap.libgap import libgap

fixtures = json.load(open("/work/tests/fixtures/lorentzian_perfect_domains.json"))
case = next(c for c in fixtures if c["id"] == "Result_Enumeration_0")
backend = LorentzianPerfectLocalBackend()
lattice_ = Lattices(ZZ)(case["gram"])
cell = backend.initial_cell(lattice_, NONNEGATIVE_SQUARES)
print(f"cell rows: {len(cell.rows)}")

started = time.monotonic()
faces = _CellFaces(backend, cell)
print(f"_CellFaces constructed in {time.monotonic() - started:.2f}s")

facets = sorted(tuple(sorted(index + 1 for index in facet)) for facet in faces.incidences)
print("facet count:", len(facets))
print("distinct facet incidence sets:", len(set(faces.facet_incidences)))
print("libgap group order:", libgap.Size(faces.group))
orbits = libgap.Orbits(faces.group, facets, libgap.OnSets)
print("facet orbits under the full cell group:", len(orbits))
print("orbit sizes:", sorted(len(orbit) for orbit in orbits))
started = time.monotonic()
flags = faces.flag_simplices(faces.full_face, faces.group)
print(f"flag_simplices(full_face) -> {len(flags)} flags in {time.monotonic() - started:.2f}s")
