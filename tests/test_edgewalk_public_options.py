"""Public edgewalk option dispatch through the maintained lattice provider."""


def test_public_edgewalk_accepts_norm_termination_and_initial_vertex_options():
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    ambient = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [0, 1, 0], [1, 0, 0], [0, 0, 2],
    ]))
    default = ambient.allcock_edgewalk()
    assert default.is_reflective()
    assert len(default.simple_roots()) == 3

    k3 = ambient.allcock_edgewalk(option_norms="K3")
    assert k3.is_reflective()
    assert all(int(root.q()) == 2 for root in k3.simple_roots())

    exhaustive = ambient.allcock_edgewalk(early_termination=False)
    assert exhaustive.is_reflective()
    assert len(exhaustive.simple_roots()) == 3

    initial = default.vertices()[0]
    ray = initial.primitive_rays()[0]
    vector_start = ambient.allcock_edgewalk(initial_vertex=ray)
    assert vector_start.is_reflective()
    assert len(vector_start.simple_roots()) == 3

    explicit_start = ambient.allcock_edgewalk(
        initial_vertex=(ray, tuple(initial.wall_roots()))
    )
    assert explicit_start.is_reflective()
    assert len(explicit_start.simple_roots()) == 3


def test_public_edgewalk_rejects_unsupported_root_norm_selection():
    import pytest
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    ambient = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [0, 1, 0], [1, 0, 0], [0, 0, 2],
    ]))
    with pytest.raises(ValueError, match="root-norm option"):
        ambient.allcock_edgewalk(option_norms="unknown")


def test_public_binary_edgewalk_accepts_explicit_starting_wall():
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    ambient = HyperbolicLattices(ZZ)(Lattices(ZZ)([[1, 0], [0, -2]]))
    default = ambient.allcock_edgewalk()
    assert default.is_reflective()
    start = default.vertices()[0]
    point = start.primitive_rays()[0]
    roots = tuple(start.wall_roots())
    assert len(roots) == 1
    explicit = ambient.allcock_edgewalk(initial_vertex=(point, roots))
    assert explicit.is_reflective()
    assert len(explicit.simple_roots()) == 2


def test_isotropic_binary_edgewalk_records_source_nonreflectivity_reason():
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    ambient = HyperbolicLattices(ZZ)(Lattices(ZZ)("U"))
    report = ambient.allcock_edgewalk()
    assert not report.is_reflective()
    assert len(report.vertices()) == 0
    assert len(report.simple_roots()) == 0
    assert report.reason_non_reflective() == "list of roots is empty"


def test_public_nonreflective_edgewalk_can_complete_without_early_termination():
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    ambient = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [0, 1, 0], [1, 0, 0], [0, 0, 46],
    ]))
    early = ambient.allcock_edgewalk()
    complete = ambient.allcock_edgewalk(early_termination=False)
    assert not early.is_reflective() and not complete.is_reflective()
    assert early.reason_non_reflective()
    assert complete.reason_non_reflective() is None
    assert len(complete.vertices()) >= len(early.vertices()) == 3


def test_public_edgewalk_isometry_uses_integral_vertex_transporters():
    from sage.matrix.constructor import matrix
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    gram = [[0, 1, 0], [1, 0, 0], [0, 0, 2]]
    source = HyperbolicLattices(ZZ)(Lattices(ZZ)(gram))
    target = HyperbolicLattices(ZZ)(Lattices(ZZ)(
        gram, module_generators=("left", "right", "root")))
    witness = source.edgewalk_isometry(target)
    assert witness is not None
    assert witness.domain() is source and witness.codomain() is target
    assert all(target.b(witness(x), witness(y)) == source.b(x, y)
               for x in source.module_generators() for y in source.module_generators())

    nonisometric = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [0, 1, 0], [1, 0, 0], [0, 0, 4],
    ]))
    assert source.edgewalk_isometry(nonisometric) is None
    odd_same_determinant = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [1, 0, 0], [0, 1, 0], [0, 0, -2],
    ]))
    assert source.edgewalk_isometry(odd_same_determinant) is None

    change = matrix([[1, 1, 0], [0, 1, 0], [0, 0, 1]])
    transformed = change * matrix(gram) * change.transpose()
    conjugate = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [int(entry) for entry in row] for row in transformed.rows()
    ]))
    conjugate_witness = source.edgewalk_isometry(conjugate)
    assert conjugate_witness is not None
    assert all(conjugate.b(conjugate_witness(x), conjugate_witness(y)) == source.b(x, y)
               for x in source.module_generators() for y in source.module_generators())


def test_public_edgewalk_isometry_of_isotropic_binary_lattices():
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    first = HyperbolicLattices(ZZ)(Lattices(ZZ)("U"))
    second = HyperbolicLattices(ZZ)(Lattices(ZZ)(
        [[0, 1], [1, 0]], module_generators=("left", "right")))
    witness = first.edgewalk_isometry(second)
    assert witness is not None
    assert witness.domain() is first and witness.codomain() is second
    assert all(second.b(witness(x), witness(y)) == first.b(x, y)
               for x in first.module_generators() for y in first.module_generators())


def test_edgewalk_isometry_rejects_equal_determinant_distinct_discriminant_groups():
    from sage.all import ZZ as SageZZ, matrix
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    first_gram = [[0, 1, 0], [1, 0, 0], [0, 0, 8]]
    second_gram = [[0, 2, 0], [2, 0, 0], [0, 0, 2]]
    first = HyperbolicLattices(ZZ)(Lattices(ZZ)(first_gram))
    second = HyperbolicLattices(ZZ)(Lattices(ZZ)(second_gram))
    assert matrix(SageZZ, first_gram).det() == matrix(SageZZ, second_gram).det() == -8
    assert tuple(matrix(SageZZ, first_gram).elementary_divisors()) == (1, 1, 8)
    assert tuple(matrix(SageZZ, second_gram).elementary_divisors()) == (2, 2, 2)
    assert first.edgewalk_isometry(second) is None


def test_public_edgewalk_isometry_of_anisotropic_binary_presentations():
    from sage.matrix.constructor import matrix
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    gram = matrix([[1, 0], [0, -2]])
    change = matrix([[1, 1], [0, 1]])
    transformed = change * gram * change.transpose()
    first = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [int(entry) for entry in row] for row in gram.rows()]))
    second = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [int(entry) for entry in row] for row in transformed.rows()]))
    witness = first.edgewalk_isometry(second)
    assert witness is not None
    assert all(second.b(witness(x), witness(y)) == first.b(x, y)
               for x in first.module_generators() for y in first.module_generators())


def test_edgewalk_missing_starting_corner_uses_exact_lorentzian_equivalence(monkeypatch):
    import sage_indefinite_port.indefinite.edgewalk as edgewalk
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    gram = [[0, 1, 0], [1, 0, 0], [0, 0, 2]]
    source = HyperbolicLattices(ZZ)(Lattices(ZZ)(gram))
    target = HyperbolicLattices(ZZ)(Lattices(ZZ)(
        gram, module_generators=("a", "b", "c")))
    monkeypatch.setattr(edgewalk._Edgewalk, "initial_vertex", lambda self: None)
    witness = source.edgewalk_isometry(target)
    assert witness is not None
    assert all(target.b(witness(x), witness(y)) == source.b(x, y)
               for x in source.module_generators() for y in source.module_generators())


def test_public_edgewalk_handles_opposite_signature_convention():
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    gram = [[0, -1, 0], [-1, 0, 0], [0, 0, -2]]
    source = HyperbolicLattices(ZZ)(Lattices(ZZ)(gram))
    target = HyperbolicLattices(ZZ)(Lattices(ZZ)(
        gram, module_generators=("a", "b", "c")))
    assert tuple(int(x) for x in source.signature_pair()) == (1, 2)
    report = source.allcock_edgewalk(option_norms="K3")
    assert report.is_reflective()
    assert len(report.simple_roots()) == 3
    assert all(int(root.q()) == -2 for root in report.simple_roots())
    start = report.vertices()[0]
    explicit = source.allcock_edgewalk(
        initial_vertex=(start.primitive_rays()[0], tuple(start.wall_roots()))
    )
    assert explicit.is_reflective()
    assert len(explicit.simple_roots()) == 3
    witness = source.edgewalk_isometry(target)
    assert witness is not None
    assert all(target.b(witness(x), witness(y)) == source.b(x, y)
               for x in source.module_generators() for y in source.module_generators())


def test_edge_root_intersection_rejects_spacelike_nonvertices():
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ
    from sage_indefinite_port.indefinite.edgewalk import _Edgewalk

    ambient = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [0, 1, 0], [1, 0, 0], [0, 0, 46],
    ]))
    walk = _Edgewalk(ambient)
    vertex = walk.initial_vertex()
    assert vertex is not None
    roots, discarded = next(
        (roots, discarded)
        for roots, discarded in walk.edge_directions(vertex)
        if len(roots) == 1 and roots[0] == 23 * ambient.module_generators()[0] - ambient.module_generators()[2]
    )
    alpha = -ambient.module_generators()[0] - ambient.module_generators()[1]
    assert alpha.q() == 2
    assert walk.vertex_from_root(vertex.generator, roots, discarded, alpha) is None


def test_anisotropic_edge_without_finite_extension_is_a_nonreflective_boundary(monkeypatch):
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ
    from sage_indefinite_port.indefinite.edgewalk import _Edgewalk

    ambient = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [1, 0, 0], [0, 1, 0], [0, 0, -3],
    ]))
    walk = _Edgewalk(ambient)
    vertex = walk.initial_vertex()
    assert vertex is not None
    roots, discarded = walk.edge_directions(vertex)[0]
    assert not walk.edge_plane(roots, discarded, vertex.generator, walk.norms[0]).plane.is_isotropic()
    monkeypatch.setattr(walk, "finite_extensions", lambda _roots: ())
    assert walk.edge_step(vertex, roots, discarded) is None
    report = walk.walk()
    assert not report.is_reflective()
    assert "anisotropic chamber edge" in report.reason_non_reflective()
