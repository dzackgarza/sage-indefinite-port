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
from sage.all import QQ, ZZ, matrix

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


def _coordinates(element, lattice_):
    """Exact coordinates of ``element`` in the lattice's declared label order.

    Owned ring elements do not coerce into Sage's `QQ` (`QQ(value)` raises "Integer Ring is not in Category of
    sets with partial maps"), so the conversion goes through the exact numerator/denominator pair that
    `groupoid_cache._gram_in_label_order` uses for the same reason.
    """
    coordinates = element.to_vector()
    return tuple(
        QQ(int(coordinates(label).numerator())) / QQ(int(coordinates(label).denominator()))
        for label in lattice_.module_generating_set()
    )


def _positive_ratio(source, target):
    r"""``λ > 0`` with ``source = λ target``, or None.

    Positive scaling, not "same line": ``v`` and ``-v`` are different rays of a pointed cone and must not be
    identified. Division is exact over QQ, so a rational ratio is never truncated.
    """
    ratio = None
    for left, right in zip(source, target, strict=True):
        if right == 0:
            if left != 0:
                return None
            continue
        candidate = left / right
        if ratio is None:
            ratio = candidate
        elif ratio != candidate:
            return None
    if ratio is None or ratio <= 0:
        return None
    return ratio


def test_normaliz_route_recovers_exactly_the_configuration_rays() -> None:
    r"""The Normaliz-built cone's rays are exactly the cell's configuration rays, up to **positive** scaling.

    `is_equal_to` compares two descriptions; this compares the geometry against the configuration that defines
    the cell, so two equally-wrong descriptions cannot both pass. The match is a bijection: every primitive ray
    of the cone is a positive multiple of exactly one configuration row, and every row is used once — so no ray
    is spurious and none is missing.
    """
    backend = LorentzianPerfectLocalBackend()
    cell = _isotropic_cell()
    required = tuple(tuple(ZZ(entry) for entry in row) for row in cell.rows)

    cone = _normaliz_cone(backend, cell)
    rays = tuple(_coordinates(ray, cell.lattice) for ray in cone.primitive_rays())

    matched: dict[int, int] = {}
    for ray in rays:
        hits = [index for index, row in enumerate(required) if _positive_ratio(ray, row) is not None]
        assert len(hits) == 1, (
            f"cone ray {ray} must be a positive multiple of exactly one configuration row, matched {hits}"
        )
        matched[hits[0]] = matched.get(hits[0], 0) + 1
    assert set(matched) == set(range(len(required))), "every configuration row must appear as a cone ray"
    assert all(count == 1 for count in matched.values()), "no configuration row may be covered twice"
    assert len(rays) == len(required)


def test_cell_cone_takes_the_ray_route_below_the_large_cell_threshold() -> None:
    """Guard the branch choice itself: a small cell must not pay for the Normaliz route."""
    backend = LorentzianPerfectLocalBackend()
    cell = _isotropic_cell()
    assert len(cell.rows) < 10_000
    cone = backend.cell_cone(cell)
    by_rays = RationalPolyhedralCones(cell.lattice).from_rays(cell.vector_configuration)
    assert cone.is_equal_to(by_rays)


def test_positive_ratio_is_sign_sensitive_and_exact() -> None:
    """Pin the matcher the ray test relies on, because getting this wrong is how the first version failed.

    A canonical "divide by the first nonzero coordinate" key identifies `v` with `-v` (both become
    `(1, ...)`) and floor-divides rational ratios; that is not positive-scaling equivalence on a pointed
    cone. These cases are the ones that separate the two notions.
    """
    assert _positive_ratio((QQ(2), QQ(4)), (QQ(1), QQ(2))) == 2
    assert _positive_ratio((QQ(2), QQ(3)), (QQ(4), QQ(6))) == QQ(1) / 2, "a rational ratio must survive exactly"
    assert _positive_ratio((QQ(-2), QQ(-4)), (QQ(1), QQ(2))) is None, "v and -v are different rays"
    assert _positive_ratio((QQ(-2), QQ(4)), (QQ(1), QQ(-2))) is None, "a consistent negative ratio is still not positive"
    assert _positive_ratio((QQ(1), QQ(1)), (QQ(1), QQ(2))) is None, "proportional entries must be consistent"
    assert _positive_ratio((QQ(0), QQ(0)), (QQ(1), QQ(2))) is None, "zero is not a ray"
    assert _positive_ratio((QQ(1), QQ(0)), (QQ(2), QQ(0))) == QQ(1) / 2
