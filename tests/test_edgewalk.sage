from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.integer_ring import ZZ

from sage_indefinite_port.indefinite.edgewalk import edgewalk_fundamental_domain


def _square(gram, row):
    v = vector(ZZ, row)
    return (v * gram * v.column())[0]


def test_allcock_edgewalk_on_hyperbolic_plane_has_one_wall() -> None:
    gram = matrix(ZZ, [[0, -1], [-1, 0]])
    record = edgewalk_fundamental_domain(gram)

    assert record["is_reflective"] is True
    assert record["simple_root_rows"] == ((1, -1),)
    assert len(record["vertices"]) == 1
    assert _square(gram, record["simple_root_rows"][0]) == 2
    assert record["isometry_generator_rows"] == ()


def test_allcock_edgewalk_on_u_plus_a1_recovers_the_reflective_triangle() -> None:
    gram = matrix(ZZ, [[0, -1, 0], [-1, 0, 0], [0, 0, 2]])
    record = edgewalk_fundamental_domain(gram)

    expected_roots = {(0, 0, -1), (1, -1, 0), (0, 1, 1)}
    expected_vertices = {(0, 1, 0), (1, 1, 0), (2, 2, 1)}

    assert record["is_reflective"] is True
    assert set(record["simple_root_rows"]) == expected_roots
    assert {generator for generator, _roots in record["vertices"]} == expected_vertices
    assert all(_square(gram, root) == 2 for root in record["simple_root_rows"])
    assert record["isometry_generator_rows"] == ()
