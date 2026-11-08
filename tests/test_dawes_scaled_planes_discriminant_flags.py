"""The discriminant quadric of 2U(2)+A2 has nine points and six lines.

Explicit primitive null planes lift every maximal isotropic plane of its
2-primary discriminant form. Their flags give the coloured incidence graph
in Dawes's diagram for this lattice; black vertices denote rational null lines.
"""

from sage.all import Graph

from tests.acceptance.consumer import basis, lattice
from tests.fixtures.oracle_fixtures import load_dawes_buildings


def test_scaled_planes_discriminant_flags_match_the_published_vertex_colours():
    fixture = load_dawes_buildings()
    case = next(case for case in fixture["buildings"] if case["id"] == "dawes_2U2_A2_stable")
    ambient = lattice(case["gram"])
    e1, f1, e2, f2, a, b = basis(ambient)
    vectors = (e1, f1, e2, f2, e1 + e2, e1 + f2, f1 + e2, f1 + f2, e1 + f1 + e2 - f2)
    assert all(point.q() == 0 and point.is_primitive() and int(point.div()) == 2 for point in vectors)
    classes = tuple(ambient.divided_discriminant_class(point) for point in vectors)
    assert len(set(classes)) == 9
    assert all(int(point.additive_order()) == 2 for point in classes)
    assert len({frozenset((point, -point)) for point in classes}) == 9
    null_bases = (
        (e1, e2), (e1, f2), (f1, e2), (f1, f2),
        (e1 + e2, f1 - f2), (e1 + f2, f1 - e2),
    )
    plane_classes = []
    for first, second in null_bases:
        plane = ambient.primitive_isotropic_subobject(first, second)
        assert int(plane.module_rank()) == 2
        assert plane.is_primitive() and plane.is_totally_isotropic()
        plane_classes.append(frozenset(ambient.divided_discriminant_class(point) for point in (first, second, first + second)))
    expected_planes = {
        frozenset((first, second, first + second))
        for first in classes for second in classes if first != second and first + second in classes
    }
    assert set(plane_classes) == expected_planes
    assert len(expected_planes) == 6
    incidences = tuple((point, plane) for plane in plane_classes for point in plane)
    assert len(incidences) == 18
    assert all(sum(point in plane for plane in plane_classes) == 2 for point in classes)
    published = next(graph for graph in fixture["incidence_graphs"] if graph["building_id"] == case["id"])
    canonical = []
    for points, curves, edges in (
        (classes, plane_classes, incidences),
        (published["curves"], published["points"], tuple((curve, point) for point, curve in published["point_curve_incidences"])),
    ):
        point_vertices = [("p", point) for point in points]
        curve_vertices = [("c", curve) for curve in curves]
        graph = Graph([point_vertices + curve_vertices, [(("p", point), ("c", curve)) for point, curve in edges]], format="vertices_and_edges")
        canonical.append(graph.canonical_label(partition=[point_vertices, curve_vertices]))
    assert canonical[0] == canonical[1]
