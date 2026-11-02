"""What exactly is a primitive ray's coordinate, and how does it convert exactly to QQ?"""
from sage.all import QQ, ZZ, matrix

from sage_extensions.categories.polyhedral_cones import _rational_polyhedral_cone
from sage_indefinite_port.indefinite.lorentzian_cells import (
    ZERO_SQUARES,
    LorentzianPerfectCell,
    LorentzianPerfectLocalBackend,
)
from tests.acceptance.consumer import lattice

gram = matrix(ZZ, [[2, -1, 0], [-1, 2, 0], [0, 0, -2]])
roots = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1))
cell = LorentzianPerfectCell(
    lattice([list(row) for row in gram.rows()]), tuple((x, y, 1) for x, y in roots), ZERO_SQUARES
)
backend = LorentzianPerfectLocalBackend()
hyperplanes, equations = backend._normaliz_cell_h_presentation(cell, include_equations=True)
cone = _rational_polyhedral_cone(
    cell.lattice,
    backend._iter_normaliz_owned_covectors(cell, hyperplanes),
    equation_covectors=backend._iter_normaliz_owned_covectors(cell, equations),
)
rays = tuple(cone.primitive_rays())
print("rays:", len(rays), "rows:", len(cell.rows), flush=True)
ray = rays[0]
coordinates = ray.to_vector()
label = next(iter(cell.lattice.module_generating_set()))
value = coordinates(label)
print("coordinate type:", type(value).__module__, type(value).__name__, "repr:", repr(value), flush=True)
for name in ("numerator", "denominator", "lift", "to_rational", "rational", "__int__"):
    print(f"  has {name}:", hasattr(value, name), flush=True)
print("int(value):", int(value), flush=True)
try:
    print("QQ(value):", QQ(value), flush=True)
except Exception as error:
    print("QQ(value) raised:", type(error).__name__, str(error)[:120], flush=True)
try:
    print("QQ(int(num)) / int(den):", QQ(int(value.numerator())) / QQ(int(value.denominator())), flush=True)
except Exception as error:
    print("numerator route raised:", type(error).__name__, str(error)[:120], flush=True)
print("all ray coordinates:", [[str(coordinates(l)) for l in cell.lattice.module_generating_set()] for coordinates in
                               (r.to_vector() for r in rays)], flush=True)
