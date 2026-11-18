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
