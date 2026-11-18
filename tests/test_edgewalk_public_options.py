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


def test_public_nonreflective_edgewalk_can_complete_without_early_termination():
    from sage_extensions.all import HyperbolicLattices, Lattices, ZZ

    ambient = HyperbolicLattices(ZZ)(Lattices(ZZ)([
        [0, 1, 0], [1, 0, 0], [0, 0, 46],
    ]))
    early = ambient.allcock_edgewalk()
    complete = ambient.allcock_edgewalk(early_termination=False)
    assert not early.is_reflective() and not complete.is_reflective()
    assert len(complete.vertices()) >= len(early.vertices()) == 3


def test_public_edgewalk_isometry_uses_integral_vertex_transporters():
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
