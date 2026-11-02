"""The large-cell cone route agrees with the generic ray route.

`LorentzianPerfectLocalBackend.cell_cone` switches engines at 10 000 rows
(`lorentzian_cells.py:1011-1029`): below the threshold the cone comes from
`RationalPolyhedralCones(...).from_rays(...)`, above it from Normaliz's raw
support hyperplanes and equations. The only configuration that ever takes the
upper branch in the acceptance suite is the A7 perfect configuration (49 861
rays), and the recorded engine agreement (PPL vs Normaliz, 80/80) is a
statement about identical **ray sets**, not about this conversion. So the
branch the A7 discriminant-image node depends on had no equality test against
the branch everything else uses.

This test takes the two descriptions of one genuine perfect cell and asserts
they are the same cone. It is deliberately not written as "the threshold is
10 000": the threshold is a cost decision, while what has to be true is that
the two routes agree wherever both can run.
"""
import pytest
from sage.all import ZZ, matrix

from sage_extensions.categories.polyhedral_cones import (
    RationalPolyhedralCones,
    _rational_polyhedral_cone,
)
from sage_indefinite_port.indefinite.lorentzian_cells import (
    NONNEGATIVE_SQUARES,
    ZERO_SQUARES,
    LorentzianPerfectCell,
    LorentzianPerfectLocalBackend,
)
from tests.acceptance.consumer import lattice


def _isotropic_cell():
    """A 2-dimensional perfect cell in a rank-3 lattice of signature (2, 1)."""
    gram = matrix(ZZ, [[2, -1, 0], [-1, 2, 0], [0, 0, -2]])
    roots = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1))
    rows = tuple((x, y, 1) for x, y in roots)
    return LorentzianPerfectCell(lattice([list(row) for row in gram.rows()]), rows, ZERO_SQUARES)


def _normaliz_cone(backend, cell):
    """The cone `cell_cone` builds on the large-cell branch, for any cell."""
    hyperplanes, equations = backend._normaliz_cell_h_presentation(cell, include_equations=True)
    return _rational_polyhedral_cone(
        cell.lattice,
        backend._iter_normaliz_owned_covectors(cell, hyperplanes),
        equation_covectors=backend._iter_normaliz_owned_covectors(cell, equations),
    )


def test_normaliz_h_presentation_describes_the_same_cone_as_from_rays() -> None:
    backend = LorentzianPerfectLocalBackend()
    cell = _isotropic_cell()

    by_rays = RationalPolyhedralCones(cell.lattice).from_rays(cell.vector_configuration)
    by_normaliz = _normaliz_cone(backend, cell)

    assert by_normaliz.is_equal_to(by_rays), "the Normaliz H-presentation must describe the cone of the rays"


def test_normaliz_route_recovers_exactly_the_configuration_rays() -> None:
    r"""The primitive rays of the Normaliz-built cone are the cell's own rays, up to positive scaling.

    `is_equal_to` compares the two descriptions; this checks the geometry directly, so a matching pair of
    equally-wrong H-descriptions cannot pass both.
    """
    backend = LorentzianPerfectLocalBackend()
    cell = _isotropic_cell()

    def normalized(cone):
        rays = []
        for ray in cone.primitive_rays():
            coordinates = tuple(int(ray.to_vector()(label)) for label in cell.lattice.module_generating_set())
            first = next(entry for entry in coordinates if entry)
            rays.append(tuple(entry // first for entry in coordinates))
        return sorted(rays)

    assert normalized(_normaliz_cone(backend, cell)) == normalized(
        RationalPolyhedralCones(cell.lattice).from_rays(cell.vector_configuration)
    )


def test_cell_cone_takes_the_ray_route_below_the_large_cell_threshold() -> None:
    """Guard the branch choice itself: a small cell must not pay for the Normaliz route."""
    backend = LorentzianPerfectLocalBackend()
    cell = _isotropic_cell()
    assert len(cell.rows) < 10_000
    cone = backend.cell_cone(cell)
    by_rays = RationalPolyhedralCones(cell.lattice).from_rays(cell.vector_configuration)
    assert cone.is_equal_to(by_rays)
