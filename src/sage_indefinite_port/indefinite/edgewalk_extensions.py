"""Normed Dynkin extensions of a root configuration, for the Allcock edgewalk (ComputePossibleExtensions)."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations

from dzack_research.preamble.all import Lattices
from dzack_research.preamble.tensors.tensor import _engine_component_matrix
from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector
from sage.rings.infinity import Infinity, PlusInfinity
from sage.rings.integer import Integer
from sage.rings.rational import Rational
from sage.rings.rational_field import QQ

from .edgewalk_scalars import _rational

type CoxeterLabel = Integer | PlusInfinity

_FOUR_COS_SQUARED: dict[CoxeterLabel, Rational] = {
    Integer(2): QQ(0),
    Integer(3): QQ(1),
    Integer(4): QQ(2),
    Integer(6): QQ(3),
    Infinity: QQ(4),
}
_FINITE_LABELS: tuple[CoxeterLabel, ...] = (
    Integer(2),
    Integer(3),
    Integer(4),
    Integer(6),
)


@dataclass(frozen=True)
class NormedDynkinExtension:
    r"""A root a of norm ``norm`` with prescribed Coxeter labels against the roots r_i.

    ``u_component`` is the projection of a onto the span of the r_i, in coordinates of L tensor Q,
    and ``residual_norm`` is q(a) - q(u_component), the norm left for the component orthogonal to them.
    """

    u_component: tuple[Rational, ...]
    residual_norm: Rational
    norm: Rational


def _norm_ratio_allows(label: CoxeterLabel, old_norm: Rational, new_norm: Rational) -> bool:
    if label == 3:
        return new_norm == old_norm
    if label == 4:
        return 2 * new_norm == old_norm or new_norm == 2 * old_norm
    if label == 6:
        return 3 * new_norm == old_norm or new_norm == 3 * old_norm
    return True


_INFINITE_LABEL = -1
r"""The Coxeter label of parallel walls in the integer rows passed below (``practical_infinity``)."""


def _recognize_irreducible_diagram(rows: tuple[tuple[int, ...], ...]) -> tuple[str, int] | None:
    r"""The type and rank of a connected spherical or Euclidean Coxeter diagram, or ``None``.

    # source: coxeter_dynkin.h:305-640 RecognizeIrreducibleSphericalEuclideanDiagram

    ``rows`` is the Coxeter matrix of the component (2 for commuting walls).
    A connected graph has exactly one cycle when it has as many edges as
    vertices and more than one when it has more (GRAPH_FindAllCycles).
    """
    size = len(rows)
    if size == 1:
        return "A", 1
    adjacency = tuple(tuple(j for j in range(size) if j != i and rows[i][j] != 2) for i in range(size))
    degrees = tuple(len(neighbours) for neighbours in adjacency)
    isolated_adjacent = tuple(rows[i][adjacency[i][-1]] if adjacency[i] else None for i in range(size))
    multiplicity: dict[int, int] = {}
    higher_edges = 0
    edges = 0
    for i in range(size):
        for j in adjacency[i]:
            if j > i:
                edges += 1
                multiplicity[rows[i][j]] = multiplicity.get(rows[i][j], 0) + 1
                if rows[i][j] != 3:
                    higher_edges += 1
    degree_one = [i for i in range(size) if degrees[i] == 1]
    degree_two = [i for i in range(size) if degrees[i] == 2]
    degree_three = [i for i in range(size) if degrees[i] == 3]
    degree_four = [i for i in range(size) if degrees[i] == 4]
    if any(degree > 4 for degree in degrees):
        return None
    if degree_four:
        if size != 5 or len(degree_four) != 1 or len(degree_one) != 4 or degree_two or degree_three:
            return None
        centre = degree_four[0]
        if any(rows[i][centre] != 3 for i in range(size) if i != centre):
            return None
        return "tildeD", 4
    if edges > size:
        return None
    if edges == size:
        if len(degree_two) != size or higher_edges:
            return None
        return "tildeA", size - 1
    if len(degree_one) == 2 and len(degree_two) == size - 2 and not higher_edges:
        return "A", size
    if len(degree_three) > 2:
        return None
    if len(degree_three) == 2:
        if higher_edges:
            return None
        for centre in degree_three:
            if sum(1 for neighbour in adjacency[centre] if degrees[neighbour] == 1) != 2:
                return None
        return "tildeD", size - 1
    if not degree_three:
        if multiplicity.get(4, 0) == 2:
            if higher_edges != 2:
                return None
            if any(isolated_adjacent[vertex] != 4 for vertex in degree_one):
                return None
            return "tildeC", size - 1
        if multiplicity.get(4, 0) == 1:
            if size == 2:
                return "B", 2
            if higher_edges != 1:
                return None
            single = sum(1 for vertex in degree_one if isolated_adjacent[vertex] == 3)
            four = sum(1 for vertex in degree_one if isolated_adjacent[vertex] == 4)
            if single == 2:
                if size == 4:
                    return "F", 4
                if size == 5:
                    return "tildeF", 4
            if four == 1 and single == 1:
                return "B", size
            return None
        if multiplicity.get(5, 0) == 1:
            if size == 2:
                return "H", 2
            if size > 5:
                return None
            single = sum(1 for vertex in degree_one if isolated_adjacent[vertex] == 3)
            five = sum(1 for vertex in degree_one if isolated_adjacent[vertex] == 5)
            if single == 1 and five == 1 and size in (3, 4):
                return "H", size
            return None
        if multiplicity.get(6, 0) == 1:
            if higher_edges != 1:
                return None
            if size == 2:
                return "G", 2
            if size == 3:
                return "tildeG", 2
            return None
        if size == 2:
            return ("tildeI", 1) if rows[0][1] == _INFINITE_LABEL else ("I", 2)
        return None
    (centre,) = degree_three
    if multiplicity.get(4, 0) == 1:
        leaves = [neighbour for neighbour in adjacency[centre] if degrees[neighbour] == 1]
        simple_leaves = [neighbour for neighbour in leaves if isolated_adjacent[neighbour] == 3]
        if len(leaves) == 3:
            return ("tildeB", 3) if len(simple_leaves) == 2 else None
        if len(leaves) != 2 or len(simple_leaves) != 2:
            return None
        if any(isolated_adjacent[vertex] == 4 for vertex in degree_one):
            return "tildeB", size - 1
        return None
    if higher_edges:
        return None
    lengths: dict[int, int] = {}
    for neighbour in adjacency[centre]:
        previous, current, length = centre, neighbour, 1
        while degrees[current] != 1:
            following = next(vertex for vertex in adjacency[current] if vertex != previous)
            previous, current = current, following
            length += 1
        lengths[length] = lengths.get(length, 0) + 1
    if lengths.get(1, 0) == 3:
        return "D", 4
    if lengths.get(1, 0) == 2:
        return "D", size
    if lengths.get(1, 0) == 1 and lengths.get(2, 0) == 2:
        return "E", 6
    if lengths.get(1, 0) == 1 and lengths.get(2, 0) == 1 and lengths.get(3, 0) == 1:
        return "E", 7
    if lengths.get(1, 0) == 1 and lengths.get(2, 0) == 1 and lengths.get(4, 0) == 1:
        return "E", 8
    if lengths.get(2, 0) == 3:
        return "tildeE", 6
    if lengths.get(1, 0) == 1 and lengths.get(3, 0) == 2:
        return "tildeE", 7
    if lengths.get(1, 0) == 1 and lengths.get(2, 0) == 1 and lengths.get(5, 0) == 1:
        return "tildeE", 8
    return None


def _diagram_extension_labels(coxeter_rows: list[list[Integer]], spherical: bool) -> tuple[tuple[CoxeterLabel, ...], ...]:
    return _cached_diagram_extension_labels(tuple(tuple(int(entry) for entry in row) for row in coxeter_rows), spherical)


@lru_cache(maxsize=4096)
def _cached_diagram_extension_labels(coxeter_rows: tuple[tuple[int, ...], ...], spherical: bool) -> tuple[tuple[CoxeterLabel, ...], ...]:
    r"""Pinned FindDiagramExtensions_Efficient, coxeter_dynkin.h:710–1671.

    Enumerate the spherical and affine Dynkin extension families in source
    order, rather than an arbitrary bounded bond search. The preamble
    Coxeter-diagram predicate identifies the existing irreducible components;
    admissibility of each extension comes from the pinned classification.
    """
    n = len(coxeter_rows)
    two, three, four, six = tuple(Integer(i) for i in (2, 3, 4, 6))
    basic: tuple[CoxeterLabel, ...] = (two,) * n
    candidates: list[tuple[CoxeterLabel, ...]] = [basic]

    def add(*bonds: tuple[int, CoxeterLabel]) -> None:
        labels = list(basic)
        for index, value in bonds:
            labels[index] = value
        candidates.append(tuple(labels))

    def one(index: int, value: CoxeterLabel = three) -> None:
        add((index, value))

    def pair(first: int, second: int, value_first: CoxeterLabel = three, value_second: CoxeterLabel = three) -> None:
        if first == second:
            raise ArithmeticError("a Dynkin extension pair must use different vertices")
        add((first, value_first), (second, value_second))

    def triple(first: int, second: int, third: int) -> None:
        if len({first, second, third}) != 3:
            raise ArithmeticError("a Dynkin extension triple must use distinct vertices")
        add((first, three), (second, three), (third, three))

    adjacency = tuple(tuple(j for j in range(n) if j != i and coxeter_rows[i][j] != two) for i in range(n))
    degrees = tuple(len(neighbors) for neighbors in adjacency)
    isolated = tuple(i for i in range(n) if degrees[i] == 0)
    components: list[tuple[int, ...]] = []
    component_index = [-1] * n
    component_sizes = [0] * n
    for first in range(n):
        if component_index[first] != -1:
            continue
        seen = {first}
        fringe = [first]
        while fringe:
            for neighbor in adjacency[fringe.pop()]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    fringe.append(neighbor)
        component = tuple(sorted(seen))
        index = len(components)
        components.append(component)
        for vertex in component:
            component_index[vertex] = index
            component_sizes[vertex] = len(component)

    # coxeter_dynkin.h:765–1034, type-dependent eligible attachments.
    # A parabolic component remains disjoint from admissible new components.
    a2: list[tuple[int, ...]] = []
    a3: list[tuple[int, ...]] = []
    an: list[tuple[int, ...]] = []
    vert_a2: list[int] = []
    vert_g2: list[int] = []
    middle_a3: list[int] = []
    ends_a3: list[int] = []
    ends_a4: list[int] = []
    ends_a5: list[int] = []
    ends_f4: list[int] = []
    center_d4: list[int] = []
    expand_b: list[int] = []
    expand_minus_one_b: list[int] = []
    nonexpand_b: list[int] = []
    expand_d: list[int] = []
    expand_minus_one_d: list[int] = []
    nonexpand_d: list[int] = []
    minus_one_a: list[int] = []
    minus_two_a: list[int] = []
    minus_three_a: list[int] = []
    dist2_e6: list[int] = []
    dist1_e6: list[int] = []
    dist3_e7: list[int] = []
    dist2_e7: list[int] = []
    dist4_e8: list[int] = []

    for component in components:
        size = len(component)
        if size == 1:
            continue
        # coxeter_dynkin.h:820-827: every component must be spherical or Euclidean.
        recognized = _recognize_irreducible_diagram(tuple(tuple(int(coxeter_rows[i][j]) for j in component) for i in component))
        if recognized is None:
            raise ValueError("The recognition did not get something Euclidean or Spherical")
        kind, _rank = recognized
        if kind == "A":
            if size == 2:
                a2.append(component)
                an.append(component)
                vert_a2.extend(component)
            else:
                ends = tuple(i for i in component if degrees[i] == 1)
                an.append(ends)
                if size == 3:
                    a3.append(ends)
                    ends_a3.extend(ends)
                    middle_a3.extend(i for i in component if degrees[i] == 2)
                if size == 4:
                    ends_a4.extend(ends)
                if size == 5:
                    ends_a5.extend(ends)
                predecessors = []
                for position, end in enumerate(ends):
                    if position == 0 or size > 3:
                        predecessor = adjacency[end][0]
                        predecessors.append(predecessor)
                        minus_one_a.append(predecessor)
                if size > 4:
                    prior = []
                    first = True
                    for predecessor in predecessors:
                        for neighbor in adjacency[predecessor]:
                            if neighbor not in ends:
                                if first or size > 5:
                                    prior.append(neighbor)
                                    minus_two_a.append(neighbor)
                                first = False
                    if size == 7:
                        first = True
                        for predecessor in prior:
                            for neighbor in adjacency[predecessor]:
                                if neighbor not in predecessors:
                                    if first:
                                        minus_three_a.append(neighbor)
                                    first = False
        elif kind == "B":
            if size == 2:
                expand_b.extend(component)
            else:
                for vertex in component:
                    if degrees[vertex] == 1:
                        bond = coxeter_rows[vertex][adjacency[vertex][0]]
                        if bond == three:
                            expand_b.append(vertex)
                            expand_minus_one_b.append(adjacency[vertex][0])
                        if bond == four:
                            nonexpand_b.append(vertex)
        elif kind == "D":
            if size == 4:
                expand_d.extend(vertex for vertex in component if degrees[vertex] == 1)
                center_d4.extend(vertex for vertex in component if degrees[vertex] == 3)
            else:
                center = next(vertex for vertex in component if degrees[vertex] == 3)
                for vertex in component:
                    if degrees[vertex] == 1:
                        if adjacency[vertex][0] == center:
                            nonexpand_d.append(vertex)
                        else:
                            expand_d.append(vertex)
                            expand_minus_one_d.append(adjacency[vertex][0])
        elif kind == "G" and size == 2:
            vert_g2.extend(component)
        elif kind == "F" and size == 4:
            ends_f4.extend(vertex for vertex in component if degrees[vertex] == 1)
        elif kind == "E":
            center = next(vertex for vertex in component if degrees[vertex] == 3)
            for neighbor in adjacency[center]:
                previous, current, length = center, neighbor, 1
                while degrees[current] != 1:
                    following = next(i for i in adjacency[current] if i != previous)
                    previous, current = current, following
                    length += 1
                if size == 6:
                    if length == 2:
                        dist2_e6.append(current)
                    if length == 1:
                        dist1_e6.append(current)
                elif size == 7:
                    if length == 3:
                        dist3_e7.append(current)
                    if length == 2:
                        dist2_e7.append(current)
                elif size == 8 and length == 4:
                    dist4_e8.append(current)

    # coxeter_dynkin.h:1036–1197: spherical one-bond extensions.
    for ends in an:
        for vertex in ends:
            one(vertex)
    for vertex in isolated:
        one(vertex)
    for ends in an:
        for vertex in ends:
            one(vertex, four)
    for vertex in isolated:
        one(vertex, four)
    for vertex in (*expand_b, *expand_d, *minus_one_a):
        one(vertex)
    for family, size in ((minus_two_a, 5), (nonexpand_d, 5), (minus_two_a, 6), (nonexpand_d, 6)):
        for vertex in family:
            if component_sizes[vertex] == size:
                one(vertex)
    for vertex in dist2_e6:
        one(vertex)
    for family, size in ((minus_two_a, 7), (nonexpand_d, 7)):
        for vertex in family:
            if component_sizes[vertex] == size:
                one(vertex)
    for vertex in dist3_e7:
        one(vertex)
    for vertex in isolated:
        one(vertex, six)
    for vertex in nonexpand_b:
        if component_sizes[vertex] == 3:
            one(vertex)

    if not spherical:
        # Affine single-bond cases, retaining the pinned insertion order.
        for vertex in expand_d:
            one(vertex, four)
        for vertex in expand_minus_one_b:
            one(vertex)
        for vertex in middle_a3:
            one(vertex, four)
        for vertex in expand_b:
            one(vertex, four)
        for vertex in (*expand_minus_one_d, *center_d4, *dist1_e6, *minus_three_a, *dist2_e7):
            one(vertex)
        for family in (nonexpand_d, minus_two_a):
            for vertex in family:
                if component_sizes[vertex] == 8:
                    one(vertex)
        for vertex in (*dist4_e8, *vert_g2):
            one(vertex)
        for vertex in vert_a2:
            one(vertex, six)
        for vertex in ends_f4:
            one(vertex)
        for vertex in nonexpand_b:
            if component_sizes[vertex] == 4:
                one(vertex)
        for vertex in isolated:
            one(vertex, Infinity)

    # coxeter_dynkin.h:1199–1453: two-bond families.
    for left_group, right_group in combinations(an, 2):
        for left in left_group:
            for right in right_group:
                pair(left, right)
    for left in isolated:
        for group in an:
            for right in group:
                pair(left, right)
    for left, right in combinations(isolated, 2):
        pair(left, right)
    for left in isolated:
        for right in isolated:
            if left != right:
                pair(left, right, four, three)
    for left in isolated:
        for group in an:
            for right in group:
                pair(left, right, four, three)
    for left in expand_b:
        for group in an:
            for right in group:
                pair(left, right)
    for left in expand_b:
        for right in isolated:
            pair(left, right)
    for left in middle_a3:
        for group in an:
            for right in group:
                if component_index[left] != component_index[right]:
                    pair(left, right)
    for left in middle_a3:
        for right in isolated:
            pair(left, right)
    for left in expand_d:
        for group in an:
            for right in group:
                pair(left, right)
    for left in expand_d:
        for right in isolated:
            pair(left, right)
    # E6 from A4+A1, E7 from A5+A1, and the A4+A2 variant.
    for left in isolated:
        for right in minus_one_a:
            if component_sizes[right] == 4:
                pair(left, right)
    for left in isolated:
        for right in minus_one_a:
            if component_sizes[right] == 5:
                pair(left, right)
    for left in vert_a2:
        for right in minus_one_a:
            if component_sizes[right] == 4:
                pair(left, right)
    for left in isolated:
        for right in nonexpand_d:
            if component_sizes[right] == 5:
                pair(left, right)
    for left in dist2_e6:
        for right in isolated:
            pair(left, right)
    for left in isolated:
        for right in minus_one_a:
            if component_sizes[right] == 6:
                pair(left, right)
    for left in minus_one_a:
        if component_sizes[left] == 4:
            for right in ends_a3:
                pair(left, right)
    for left in nonexpand_d:
        if component_sizes[left] == 5:
            for right in vert_a2:
                pair(left, right)
    for left in vert_a2:
        for right in isolated:
            pair(left, right, four, three)

    if not spherical:
        # Affine two-bond cases from pinned coxeter_dynkin.h.
        for ends in an:
            pair(ends[0], ends[1])
        for left in expand_b:
            for right in middle_a3:
                pair(left, right)
        for left in isolated:
            for right in middle_a3:
                pair(left, right, four, three)
        for left in expand_b:
            for right in expand_d:
                pair(left, right)
        for left in isolated:
            for right in expand_d:
                pair(left, right, four, three)
        for left, right in combinations(expand_b, 2):
            if component_index[left] != component_index[right]:
                pair(left, right)
        for left in isolated:
            for right in expand_b:
                pair(left, right, four, three)
        for left, right in combinations(isolated, 2):
            pair(left, right, four, four)
        for left in isolated:
            for right in isolated:
                if component_index[left] != component_index[right]:
                    pair(left, right, three, six)
        for left, right in combinations(middle_a3, 2):
            pair(left, right)
        for left in middle_a3:
            for right in expand_d:
                pair(left, right)
        for left, right in combinations(expand_d, 2):
            if component_index[left] != component_index[right]:
                pair(left, right)
        for left in minus_two_a:
            if component_sizes[left] == 5:
                for right in isolated:
                    pair(left, right)
        for left in minus_one_a:
            if component_sizes[left] == 5:
                for right in vert_a2:
                    pair(left, right)
        for left in nonexpand_d:
            if component_sizes[left] == 6:
                for right in isolated:
                    pair(left, right)
        for left in minus_one_a:
            if component_sizes[left] == 7:
                for right in isolated:
                    pair(left, right)
        for left in minus_one_a:
            if component_sizes[left] == 4:
                for right in ends_a4:
                    if component_index[left] != component_index[right]:
                        pair(left, right)
        for left in nonexpand_d:
            if component_sizes[left] == 5:
                for right in ends_a3:
                    pair(left, right)
        for left in dist2_e6:
            for right in vert_a2:
                pair(left, right)
        for left in dist3_e7:
            for right in isolated:
                pair(left, right)
        for left in ends_a3:
            for right in isolated:
                pair(left, right, four, three)
        for left in vert_a2:
            for right in vert_a2:
                if component_index[left] != component_index[right]:
                    pair(left, right, four, three)
        for left in nonexpand_b:
            if component_sizes[left] == 3:
                for right in isolated:
                    pair(left, right)

    # coxeter_dynkin.h:1454–1671: exceptional triple and quadruple bonds.
    for first_group, second_group in combinations(a2, 2):
        for first in first_group:
            for second in second_group:
                for third in isolated:
                    triple(first, second, third)
    for first in vert_a2:
        for second in ends_a3:
            for third in isolated:
                triple(first, second, third)
    for first in vert_a2:
        for second in ends_a4:
            for third in isolated:
                triple(first, second, third)
    for first, second, third in combinations(isolated, 3):
        triple(first, second, third)
    for first, second in combinations(isolated, 2):
        for group in an:
            for third in group:
                triple(first, second, third)
    if not spherical:
        for first_group, second_group, third_group in combinations(a2, 3):
            for first in first_group:
                for second in second_group:
                    for third in third_group:
                        triple(first, second, third)
        for first_group, second_group in combinations(a3, 2):
            for first in first_group:
                for second in second_group:
                    for third in isolated:
                        triple(first, second, third)
        for first in vert_a2:
            for second in ends_a5:
                for third in isolated:
                    triple(first, second, third)
        for vertices in combinations(isolated, 3):
            for singled in vertices:
                add(*((vertex, four if vertex == singled else three) for vertex in vertices))
        for first, second in combinations(isolated, 2):
            for third in expand_b:
                triple(first, second, third)
        for first, second in combinations(isolated, 2):
            for third in expand_d:
                triple(first, second, third)
        for first, second in combinations(isolated, 2):
            for third in middle_a3:
                triple(first, second, third)
        for four_vertices in combinations(isolated, 4):
            add(*((vertex, three) for vertex in four_vertices))

    # The source inserts these classified families directly.  Its CheckDiagram
    # test is conditional on SANITY_CHECK_EFFICIENT_ENUMERATION, not part of
    # FindDiagramExtensions_Efficient's execution path: reconstructing an
    # owned Coxeter diagram for every candidate changes the algorithm and
    # dominates large edgewalks.
    return tuple(candidates)


def _extensions(
    lattice: Lattices.ParentMethods,
    roots: Sequence[Lattices.ElementMethods],
    norms: Sequence[int],
    labels: tuple[CoxeterLabel, ...],
) -> tuple[NormedDynkinExtension, ...]:
    if not roots:
        return tuple(NormedDynkinExtension((QQ.zero(),) * int(lattice.rank()), _rational(norm), _rational(norm)) for norm in norms)
    coordinates = matrix([[_rational(root.to_vector()(label)) for label in lattice.module_generating_set()] for root in roots])
    gram = _engine_component_matrix(lattice.gram_tensor()).change_ring(QQ)
    pairings = coordinates * gram * coordinates.T
    root_norms = tuple(pairings[i, i] for i in range(len(roots)))
    if not pairings.is_positive_definite():
        raise ValueError(f"the roots {roots} span a form that is not positive definite: they are not the walls through a point of hyperbolic space")
    inverse = pairings.inverse()
    if any(pairings[i, j] > 0 for i in range(len(roots)) for j in range(i)):
        raise ValueError("Dynkin extensions require simple roots with nonpositive mutual pairings")
    coxeter_rows = [[Integer(1) if i == j else Integer(2) for j in range(len(roots))] for i in range(len(roots))]
    for i in range(len(roots)):
        for j in range(i):
            old_pair = pairings[i, j]
            squares = 4 * old_pair * old_pair / (root_norms[i] * root_norms[j])
            label = next((edge for edge, value in _FOUR_COS_SQUARED.items() if squares == value), None)
            if label is None:
                raise ValueError("the supplied roots have noncrystallographic Coxeter pairings")
            coxeter_rows[i][j] = coxeter_rows[j][i] = Integer(_INFINITE_LABEL) if isinstance(label, PlusInfinity) else label

    # ComputePossibleExtensions, coxeter_dynkin.h:2210–2336: each diagram is
    # visited before its compatible norms, in the supplied norm order.
    out: list[NormedDynkinExtension] = []
    for diagram in _diagram_extension_labels(coxeter_rows, spherical=labels == _FINITE_LABELS):
        for new_norm in (_rational(norm) for norm in norms):
            if any(not _norm_ratio_allows(label, old_norm, new_norm) for label, old_norm in zip(diagram, root_norms, strict=True)):
                continue
            scalar_products: list[Rational] = []
            for label, old_norm in zip(diagram, root_norms, strict=True):
                square = _FOUR_COS_SQUARED[label] * old_norm * new_norm / 4
                if not square.is_square():
                    break
                scalar_products.append(-square.sqrt())
            if len(scalar_products) != len(roots):
                continue
            weights = inverse * vector(scalar_products)
            projected_square = _rational(vector(scalar_products).dot_product(weights))
            out.append(NormedDynkinExtension(tuple(weights * coordinates), new_norm - projected_square, new_norm))
    return tuple(out)


def finite_type_extensions(
    lattice: Lattices.ParentMethods,
    roots: Sequence[Lattices.ElementMethods],
    norms: Sequence[int],
) -> tuple[NormedDynkinExtension, ...]:
    r"""Extensions by a root at angles pi/2, pi/3, pi/4, pi/6 with a nonzero orthogonal component."""
    return tuple(extension for extension in _extensions(lattice, roots, norms, _FINITE_LABELS) if extension.residual_norm != 0)


def parabolic_extensions(
    lattice: Lattices.ParentMethods,
    roots: Sequence[Lattices.ElementMethods],
    norms: Sequence[int],
) -> tuple[NormedDynkinExtension, ...]:
    r"""Extensions by a root lying in the span of the roots (residual norm 0), at any crystallographic angle including parallel walls."""
    return tuple(extension for extension in _extensions(lattice, roots, norms, (*_FINITE_LABELS, Infinity)) if extension.residual_norm == 0)
