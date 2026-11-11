r"""Options of Allcock's edgewalk through the public ``HyperbolicLattices`` API (issue #15).

Each case is the first example that issue #15 names for the option:

- the full walk (``early_termination=False``, source EarlyTerminationIfNotReflective=false)
  on U+<46>: every edge of every listed vertex ends at a vertex equivalent to a listed vertex;
- ``option_norms="K3"`` on U+<4>: (0,0,1) is a norm-4 wall in "all" mode, and K3 mode
  keeps only norm-2 walls;
- a prescribed initial vertex (source FileVertex and FileVertexRoots) on U+<2>;
- ``reason_non_reflective`` on U+<46>: the finiteness test stops the walk;
- ``edgewalk_isometry`` from U+<2> to a unimodular basis change of it.

The last case covers ``edge_directions`` at an ideal vertex whose root cone is not
simplicial (pinned edgewalk.h:525-539 GetAdjacencyDirection), so that a facet omits
more than one wall.
"""

from sage_extensions.all import ZZ as PreambleZZ
from sage_extensions.all import HyperbolicLattices
from sage_extensions.categories.polyhedral_cones import RationalPolyhedralCones
from sage.all import ZZ, matrix

from sage_indefinite_port.indefinite.edgewalk import _Edgewalk, _Vertex
from tests.acceptance.consumer import coordinates, lattice


def _hyperbolic(gram):
    return HyperbolicLattices(PreambleZZ)(lattice(gram))


def _u_plus(value):
    return [[int(0), int(1), int(0)], [int(1), int(0), int(0)], [int(0), int(0), int(value)]]


def _root_gram_invariant(lattice_, roots):
    """The sorted diagonal and the sorted off-diagonal entries of the roots' Gram matrix."""
    diagonal = sorted(int(root.q()) for root in roots)
    off_diagonal = sorted(int(lattice_.b(left, right)) for index, left in enumerate(roots) for right in roots[index + 1 :])
    return diagonal, off_diagonal


def _listed_vertex(walker, cone):
    (generator,) = tuple(cone.primitive_rays())
    return _Vertex(generator, tuple(cone.wall_roots()))


def test_full_walk_on_u_plus_46_closes_every_edge():
    """Every edge neighbour of a listed vertex is O(L)-equivalent to a listed vertex."""
    lattice_ = _hyperbolic(_u_plus(46))
    report = lattice_.allcock_edgewalk(early_termination=False)

    assert not report.is_reflective()
    walker = _Edgewalk(lattice_, early_termination=False)
    listed = tuple(_listed_vertex(walker, cone) for cone in report.vertices())
    assert len(listed) >= 3
    listed_configurations = tuple(walker.vertex_configuration(vertex) for vertex in listed)
    for vertex in listed:
        for roots, discarded in walker.edge_directions(vertex):
            neighbour = walker.edge_step(vertex, roots, discarded)
            neighbour_configuration = walker.vertex_configuration(neighbour)
            assert any(
                lattice_.O().transporter(known, neighbour_configuration) is not None for known in listed_configurations
            ), f"the edge {roots} of {vertex.generator} ends at {neighbour.generator}, which is not equivalent to a listed vertex"


def test_k3_norms_drop_the_norm_four_wall_of_u_plus_4():
    """(0,0,1) has norm 4 and is a wall in "all" mode; K3 mode keeps only norm-2 walls."""
    lattice_ = _hyperbolic(_u_plus(4))
    all_roots = tuple(coordinates(root) for root in lattice_.allcock_edgewalk().simple_roots())
    k3_report = lattice_.allcock_edgewalk(option_norms="K3")
    k3_roots = tuple(coordinates(root) for root in k3_report.simple_roots())

    assert (0, 0, 1) in all_roots
    assert (0, 0, 1) not in k3_roots
    assert all(int(root.q()) == 2 for root in k3_report.simple_roots())


def test_prescribed_initial_vertex_on_u_plus_2():
    """A walk from a listed vertex, with or without its roots, finds the same chamber."""
    lattice_ = _hyperbolic(_u_plus(2))
    report = lattice_.allcock_edgewalk()
    assert report.is_reflective()
    expected = _root_gram_invariant(lattice_, tuple(report.simple_roots()))

    for cone in report.vertices():
        (generator,) = tuple(cone.primitive_rays())
        for initial in (generator, (generator, tuple(cone.wall_roots()))):
            prescribed = lattice_.allcock_edgewalk(initial_vertex=initial)
            assert prescribed.is_reflective()
            assert _root_gram_invariant(lattice_, tuple(prescribed.simple_roots())) == expected
            starts = tuple(tuple(vertex.primitive_rays())[0] for vertex in prescribed.vertices())
            assert coordinates(starts[0]) == coordinates(generator)


def test_u_plus_46_reports_the_finiteness_reason():
    """The early-terminating walk on U+<46> stops at the finiteness test."""
    report = _hyperbolic(_u_plus(46)).allcock_edgewalk()

    assert not report.is_reflective()
    assert report.reason_non_reflective() == "Not reflective by the finiteness test"


def test_edgewalk_isometry_of_u_plus_2_under_a_basis_change():
    """The edgewalk route returns an isometry U+<2> -> P (U+<2>) P^T, and it is one."""
    gram = matrix(ZZ, _u_plus(2))
    change = matrix(ZZ, [[int(1), int(1), int(0)], [int(0), int(1), int(1)], [int(0), int(0), int(1)]])
    assert change.det() == 1
    changed = change * gram * change.transpose()
    first = _hyperbolic(gram.rows())
    second = _hyperbolic([[int(entry) for entry in row] for row in changed.rows()])

    witness = first.edgewalk_isometry(second)

    assert witness is not None
    basis = tuple(first.module_generators())
    for left in basis:
        for right in basis:
            assert second.b(witness(left), witness(right)) == first.b(left, right)
    images = matrix(ZZ, [list(coordinates(witness(generator))) for generator in basis])
    assert images.det() in (1, -1)


def test_ideal_facets_omitting_several_walls_give_one_neighbour():
    r"""At the cusp of U+<2>+<2> under K3 norms every omitted wall orients the edge alike.

    k^perp/k is the even binary lattice of determinant 4, which is A1+A1, so the cusp's
    walls form two affine A1 components: four roots in the rank-three space k^perp. Each
    facet of their cone contains one root of each component and omits the other two.
    """
    lattice_ = _hyperbolic(
        [
            [int(0), int(1), int(0), int(0)],
            [int(1), int(0), int(0), int(0)],
            [int(0), int(0), int(2), int(0)],
            [int(0), int(0), int(0), int(2)],
        ]
    )
    walker = _Edgewalk(lattice_, option_norms="K3")
    vertex = walker.initial_isotropic_vertex()
    assert vertex is not None
    assert int(vertex.generator.q()) == 0
    assert len(vertex.roots) == 4

    cone = RationalPolyhedralCones(lattice_).from_rays(vertex.roots)
    several = 0
    for covector in cone.facet_covectors():
        included = tuple(root for root in vertex.roots if cone.evaluate_covector(covector, root) == 0)
        omitted = tuple(root for root in vertex.roots if cone.evaluate_covector(covector, root) != 0)
        several += len(omitted) > 1
        neighbours = {coordinates(walker.edge_step(vertex, included, discarded).generator) for discarded in omitted}
        assert len(neighbours) == 1, f"the facet {included} has neighbours {neighbours} for different omitted walls"
    assert several == 4
