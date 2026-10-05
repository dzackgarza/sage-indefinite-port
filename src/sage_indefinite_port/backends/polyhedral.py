"""Facet incidence and finite facet actions for Lorentzian cells."""

from __future__ import annotations

from functools import cache

from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods
from dzack_research.preamble.categories.polyhedral_cones import RationalPolyhedralCones
from sage.libs.gap.libgap import libgap

from sage_indefinite_port.backends.canonization import (
    CellConfiguration,
    _coordinate_matrix,
)
from sage_indefinite_port.groups.integral_structures import RationalMatrixGroup

type FacetIncidence = frozenset[int]


@cache
def configuration_cone(
    configuration: CellConfiguration,
) -> RationalPolyhedralCones.ParentMethods:
    """Return the exact preamble cone spanned by the configuration vectors."""
    return RationalPolyhedralCones(configuration.lattice).from_rays(configuration.vectors)


@cache
def configuration_facets(
    configuration: CellConfiguration,
) -> tuple[FacetIncidence, ...]:
    """Return facets as incidence subsets of the original vector family."""
    cone = configuration_cone(configuration)
    coordinates = _coordinate_matrix(configuration)
    facets = tuple(
        frozenset(position for position, value in enumerate(coordinates * inequality.A().column()) if value == 0) for inequality in cone._engine_polyhedron().inequalities()
    )
    if any(not facet for facet in facets):
        raise ArithmeticError("a perfect-cell facet has no incident configuration vectors")
    return tuple(sorted(facets, key=lambda facet: tuple(sorted(facet))))


def facet_orbits(
    configuration: CellConfiguration,
    stabilizer: RationalMatrixGroup,
) -> tuple[tuple[FacetIncidence, ...], ...]:
    """Return facet orbits under an exact cell stabilizer via libGAP."""
    facets = configuration_facets(configuration)
    if not facets:
        return ()
    facet_position = {facet: position for position, facet in enumerate(facets)}
    gap_generators = []
    for generator in stabilizer.generators():
        permutation = _configuration_permutation(configuration, generator)
        images = []
        for facet in facets:
            moved = frozenset(permutation[position] for position in facet)
            target = facet_position.get(moved)
            if target is None:
                raise ArithmeticError("a cell-stabilizer generator does not preserve the facet set")
            images.append(target + 1)
        gap_generators.append(libgap.PermList(images))
    group = libgap.Group(gap_generators)
    domain = libgap(list(range(1, len(facets) + 1)))
    return tuple(tuple(facets[int(position) - 1] for position in orbit) for orbit in libgap.Orbits(group, domain).sage())


@cache
def _configuration_permutation(
    configuration: CellConfiguration,
    generator: LatticeIsometryMethods,
) -> tuple[int, ...]:
    coordinates = _coordinate_matrix(configuration)
    position_by_row = {tuple(row): position for position, row in enumerate(coordinates.rows())}
    if len(position_by_row) != len(configuration.vectors):
        raise ArithmeticError("a cell configuration contains duplicate coordinate rows")
    action = configuration.lattice.Aut()._row_action_matrix(generator)
    images = []
    for row in (coordinates * action).rows():
        position = position_by_row.get(tuple(row))
        if position is None:
            raise ArithmeticError("a cell-stabilizer generator does not preserve the configuration vectors")
        images.append(position)
    if len(set(images)) != len(configuration.vectors):
        raise ArithmeticError("a cell-stabilizer generator does not act bijectively on the configuration")
    return tuple(images)
