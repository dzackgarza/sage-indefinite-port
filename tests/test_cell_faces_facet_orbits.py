"""Baseline of the cell facet-orbit computation, before the GAP calls are replaced.

`_CellFaces` currently computes its facet stabilizers and facet orbits with direct libGAP
calls (`libgap.Group`, `libgap.Stabilizer`, `libgap.Orbits` with `OnSets`), while issue
#15's boundary assigns cell stabilizers and facet orbits to research's maintained
operations. Replacing those calls has to preserve the combinatorics exactly, so this file
pins what the current path produces on a small frozen cell.

These assertions are about the *result*, not the implementation: any substitution that
keeps the facet set, the incidence classes, the cell group's order, the facet-orbit
decomposition and the flag count is free to change how they are computed.

The numbers are from `Result_Enumeration_0` (4 rays, 4 facets): group order 8, one facet
orbit of size 4, one maximal flag.
"""
import json
from pathlib import Path

import pytest

from sage.libs.gap.libgap import libgap

from sage_extensions.all import Lattices, ZZ

from sage_indefinite_port.indefinite.lorentzian_cells import (
    NONNEGATIVE_SQUARES,
    LorentzianPerfectLocalBackend,
    _CellFaces,
)

FIXTURES = Path(__file__).parent / "fixtures" / "lorentzian_perfect_domains.json"
CASE = "Result_Enumeration_0"


@pytest.fixture(scope="module")
def faces():
    with open(FIXTURES, encoding="utf-8") as handle:
        cases = json.load(handle)
    case = next(case for case in cases if case["id"] == CASE)
    backend = LorentzianPerfectLocalBackend()
    cell = backend.initial_cell(Lattices(ZZ)(case["gram"]), NONNEGATIVE_SQUARES)
    return _CellFaces(backend, cell)


def test_facet_set_and_incidences(faces):
    assert len(faces.rows) == 4
    assert len(faces.incidences) == 4
    assert len(set(faces.facet_incidences)) == 4, "distinct facets must have distinct ray incidences"
    assert faces.all_facet_positions == frozenset(range(4))


def test_cell_group_order(faces):
    assert int(libgap.Size(faces.group)) == 8


def test_facet_orbits_under_the_cell_group(faces):
    facets = sorted(tuple(sorted(index + 1 for index in facet)) for facet in faces.incidences)
    orbits = libgap.Orbits(faces.group, facets, libgap.OnSets)
    assert len(orbits) == 1
    assert sorted(len(orbit) for orbit in orbits) == [4]


def test_flag_count_under_the_cell_group(faces):
    flags = faces.flag_simplices(faces.full_face, faces.group)
    assert len(flags) == 1


def test_facet_stabilizer_is_a_subgroup_of_the_cell_group(faces):
    facet = sorted(index + 1 for index in sorted(faces.incidences)[0])
    stabilizer = libgap.Stabilizer(faces.group, facet, libgap.OnSets)
    assert int(libgap.Size(faces.group)) % int(libgap.Size(stabilizer)) == 0
    assert int(libgap.Size(stabilizer)) == 2
