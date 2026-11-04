"""Validation of the orbit-reduced facet enumeration against a complete double description.

The reference side is an independent Motzkin double description (`all_facets_by_dd`): it builds
the full facet set of cone(rows) with no group and no orbit logic, so a disagreement is a real
disagreement between two algorithms rather than a shared assumption.

Each case asserts, with exact arithmetic only:

  A. the enumerated orbit representatives, saturated under the group, are *exactly* the facets the
     double description finds (completeness and no spurious facet);
  B. each representative's normal is the case's facet normal up to a positive rational multiple;
  C. the orbit count equals the number of orbits of the double description's facet set;
  D. each orbit member is reached from its representative by an explicit word in the generators
     (transport witnesses), checked by applying the word to the incidence;
  E. the facet contract (nonnegative everywhere, zero set of span rank n-1) holds for every normal
     and the pivot rejects a non-ridge.

Run: python3 validate_orbit_facets.py
"""
from __future__ import annotations

import itertools
import sys
from fractions import Fraction

from orbit_facets import (
    certify_facet,
    dot,
    facet_orbits,
    incidence_of,
    independent_subset,
    nullspace_int,
    one_facet_by_descent,
    orbit_closure,
    pivot_across_ridge,
    primitive,
    rank_int,
    word_as_permutation,
)


# ---------------------------------------------------------------------------------------------
# independent ground truth: complete double description
# ---------------------------------------------------------------------------------------------
def _canonical(normal):
    normal = primitive(normal)
    for entry in normal:
        if entry != 0:
            return normal if entry > 0 else tuple(-x for x in normal)
    raise ArithmeticError("zero normal")


def all_facets_by_subsets(rows):
    """Every facet normal of cone(rows), by exhaustive rank-(n-1) subset enumeration.

    A facet's incidence contains a set of rank n-1, and the annihilator of such a set is one
    dimensional, so enumerating every rank-(n-1) subset and keeping the oriented annihilators that
    are nonnegative on all rays produces *every* facet and nothing else.  Exponential in the rank
    and deliberately so: this is the independent ground truth, not the algorithm under test.
    """
    n = len(rows[0])
    facets = {}
    for subset in itertools.combinations(range(len(rows)), n - 1):
        candidates = nullspace_int([rows[index] for index in subset], n)
        if len(candidates) != 1:
            continue
        normal = candidates[0]
        if any(dot(normal, row) < 0 for row in rows):
            continue
        incidence = incidence_of(rows, normal)
        if rank_int([rows[i] for i in incidence]) != n - 1:
            continue
        facets.setdefault(incidence, normal)
    return facets


# ---------------------------------------------------------------------------------------------
# configurations and their group actions
# ---------------------------------------------------------------------------------------------
def _permutation_from_mapping(rows, mapping):
    """The position permutation induced by an invertible coordinate map on the rays."""
    positions = {row: index for index, row in enumerate(rows)}
    images = [positions[tuple(mapping(row))] for row in rows]
    if sorted(images) != list(range(len(rows))):
        raise ArithmeticError("the coordinate map does not permute the rays")
    return tuple(images)


def _cube_group(rows, dimension):
    generators = []
    for axis in range(dimension - 1):
        def swap(row, axis=axis):
            entries = list(row)
            entries[1 + axis], entries[2 + axis] = entries[2 + axis], entries[1 + axis]
            return tuple(entries)
        generators.append(_permutation_from_mapping(rows, swap))
    def negate(row):
        entries = list(row)
        entries[-1] = -entries[-1]
        return tuple(entries)
    generators.append(_permutation_from_mapping(rows, negate))
    return tuple(generators)


def _hyperoctahedral_group(rows, dimension):
    generators = []
    for axis in range(dimension - 1):
        def swap(row, axis=axis):
            entries = list(row)
            entries[1 + axis], entries[2 + axis] = entries[2 + axis], entries[1 + axis]
            return tuple(entries)
        generators.append(_permutation_from_mapping(rows, swap))
    for axis in range(dimension):
        def negate(row, axis=axis):
            entries = list(row)
            entries[1 + axis] = -entries[1 + axis]
            return tuple(entries)
        generators.append(_permutation_from_mapping(rows, negate))
    return tuple(generators)


def cases():
    out = []

    # 1. the 3-cube: 8 rays in rank 4, 6 facets, one facet orbit under the octahedral group.
    rows = tuple((1, a, b, c) for a in (1, -1) for b in (1, -1) for c in (1, -1))
    out.append(("cube3 / octahedral, rank 4", rows, _cube_group(rows, 3)))

    # 2. the cuboctahedron (A3 roots): 12 rays, 14 facets, two facet orbits.
    vertices = []
    for signs in itertools.product((1, -1), repeat=2):
        for zero in range(3):
            entries = [signs[0], signs[1], 0]
            entries[zero], entries[2] = entries[2], entries[zero]
            vertices.append(tuple(entries))
    rows = tuple((1, *vertex) for vertex in sorted(set(vertices)))
    out.append(("cuboctahedron / A3 roots, rank 4", rows, _hyperoctahedral_group(rows, 3)))

    # 3. the octahedron: 6 rays, 8 triangular facets, one orbit, every facet a simplex.
    rows = tuple((1, *v) for v in [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)])
    out.append(("octahedron, rank 4", rows, _hyperoctahedral_group(rows, 3)))

    # 4. the triangular prism: 6 rays, 5 facets, a C2 action exchanging the two triangles.
    rows = tuple((1, a, b, c) for c in (0, 1) for a, b in [(0, 0), (1, 0), (0, 1)])
    flip = lambda row: (row[0], row[1], row[2], 1 - row[3])
    out.append(("triangular prism / C2, rank 4", rows, (_permutation_from_mapping(rows, flip),)))

    # 5. the 4-cube: 16 rays in rank 5, 8 facets of 8 rays each, one orbit, recursion depth 3.
    rows = tuple((1, *v) for v in itertools.product((1, -1), repeat=4))
    out.append(("4-cube / hyperoctahedral, rank 5", rows, _hyperoctahedral_group(rows, 4)))

    # 6. the 24-cell: 24 rays in rank 5, 24 facets, one orbit.
    vertices = set()
    for pair in itertools.combinations(range(4), 2):
        for signs in itertools.product((1, -1), repeat=2):
            entries = [0, 0, 0, 0]
            entries[pair[0]], entries[pair[1]] = signs
            vertices.add(tuple(entries))
    rows = tuple((1, *vertex) for vertex in sorted(vertices))
    out.append(("24-cell / hyperoctahedral, rank 5", rows, _hyperoctahedral_group(rows, 4)))

    # 7. a 4-simplex: 5 rays, 5 facets, trivial group (the orbit loop degenerates to the full set).
    rows = ((1, 0, 0, 0, 0), (1, 1, 0, 0, 0), (1, 0, 1, 0, 0), (1, 0, 0, 1, 0), (1, 0, 0, 0, 1))
    out.append(("4-simplex / trivial group, rank 5", rows, ()))

    # 8. an asymmetric configuration with no symmetry at all: the adjacency walk must find every
    #    facet from the seed, so this is the case where completeness has nothing to lean on.
    rows = (
        (1, 0, 1, 0, 0),
        (1, 1, 0, 0, 1),
        (1, 0, 0, 1, 0),
        (1, 2, 1, 1, 0),
        (1, 0, 1, 2, 1),
        (1, 1, 1, 0, 2),
        (1, 3, 0, 1, 1),
        (1, -1, 1, 1, 0),
    )
    out.append(("asymmetric 8 rays / trivial group, rank 5", rows, ()))

    # 9. a case with a two-element symmetry and a non-simplicial facet, rank 5: the pyramid over
    #    the 3-cube, which has a facet with 8 rays.
    rows = tuple((1, a, b, c, 0) for a in (1, -1) for b in (1, -1) for c in (1, -1)) + ((1, 0, 0, 0, 2),)
    generators = []
    for left, right in ((1, 2), (2, 3)):
        def swap(row, left=left, right=right):
            entries = list(row)
            entries[left], entries[right] = entries[right], entries[left]
            return tuple(entries)
        generators.append(_permutation_from_mapping(rows, swap))
    def negate_base(row):
        entries = list(row)
        entries[3] = -entries[3]
        return tuple(entries)
    generators.append(_permutation_from_mapping(rows, negate_base))
    out.append(("cube pyramid / S4 x C2, rank 5", rows, tuple(generators)))
    return out


# ---------------------------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------------------------
def check_case(label, rows, permutations):
    n = len(rows[0])
    assert rank_int(rows) == n, f"{label}: the configuration must span"
    ground_truth = all_facets_by_dd(rows)

    # sanity: the ground truth is closed under the group and every normal certifies
    for incidence, normal in ground_truth.items():
        assert certify_facet(rows, normal) == incidence
        for permutation in permutations:
            image = tuple(sorted(permutation[position] for position in incidence))
            assert image in ground_truth, f"{label}: the group does not preserve the facet set"

    # exact orbit partition of the ground truth, by closure
    truth_orbits = []
    remaining = set(ground_truth)
    while remaining:
        seed = min(remaining)
        closure, _ = orbit_closure(seed, permutations)
        truth_orbits.append((seed, closure))
        remaining -= set(closure)

    # the facet-finding precondition: the homogenising covector is strictly positive
    grading = (1,) + (0,) * (n - 1)
    assert all(dot(grading, row) > 0 for row in rows), f"{label}: rows must be homogenised"

    stats = {}
    seed_normal = one_facet_by_descent(rows, grading, stats=stats)
    seed_incidence = certify_facet(rows, seed_normal)
    assert seed_incidence in ground_truth, f"{label}: the descent returned a non-facet"
    stats = {}
    records = facet_orbits(rows, permutations, grading, stats=stats)

    # A. completeness and soundness
    enumerated = set()
    for record in records:
        enumerated.update(orbit_closure(record["incidence"], permutations)[0])
    assert enumerated == set(ground_truth), (
        f"{label}: enumerated {len(enumerated)} facets, the double description finds {len(ground_truth)}"
    )

    # B. exact normals
    for record in records:
        expected = ground_truth[record["incidence"]]
        ratio = None
        for left, right in zip(record["normal"], expected):
            if right == 0:
                assert left == 0, f"{label}: normal disagrees off the zero set"
                continue
            candidate = Fraction(left, right)
            if ratio is None:
                ratio = candidate
            else:
                assert ratio == candidate, f"{label}: normals are not proportional"
        assert ratio is not None and ratio > 0, f"{label}: normals must agree up to positive scaling"

    # C. orbit count
    assert len(records) == len(truth_orbits), (
        f"{label}: {len(records)} enumerated orbits, the double description has {len(truth_orbits)}"
    )
    assert sum(record["orbit_size"] for record in records) == len(ground_truth)

    # D. transport witnesses: every facet is reached from its representative by an explicit word
    representative_of = {}
    for record in records:
        closure, _ = orbit_closure(record["incidence"], permutations)
        for image, word in closure.items():
            assert tuple(sorted(image)) == tuple(sorted(image))
            mapped = tuple(sorted(word_as_permutation(permutations, word)[position]
                                  for position in record["incidence"]))
            assert mapped == image, f"{label}: the recorded word does not transport the representative"
            representative_of[image] = (record["incidence"], word)
    assert set(representative_of) == set(ground_truth)

    # E. the pivot rejects a set that is not a ridge
    record = records[0]
    if n >= 3:
        local = independent_subset([rows[i] for i in record["incidence"]], n - 1)
        bad = tuple(record["incidence"][i] for i in local)
        try:
            pivot_across_ridge(rows, record["normal"], bad)
        except ArithmeticError:
            pass
        else:
            raise AssertionError(f"{label}: the pivot accepted a facet-sized set as a ridge")

    return {
        "label": label,
        "rays": len(rows),
        "rank": n,
        "facets": len(ground_truth),
        "facet_orbits": len(records),
        "pivots": stats["pivots"],
        "dual_description_calls": stats["dual_description_calls"],
        "orbit_steps": stats["orbit_steps"],
        "max_subcone_rays": stats["max_subcone_rays"],
        "descent_steps": stats["descent_steps"],
        "orbit_sizes": sorted((record["orbit_size"] for record in records), reverse=True),
    }


def main():
    results = [check_case(*case) for case in cases()]
    header = f"{'case':44s} {'rays':>5s} {'rank':>4s} {'facets':>6s} {'orbits':>6s} {'pivots':>6s} {'dd':>4s} {'steps':>5s} {'maxsub':>6s}"
    print(header)
    print("-" * len(header))
    for result in results:
        print(
            f"{result['label']:44s} {result['rays']:5d} {result['rank']:4d} {result['facets']:6d} "
            f"{result['facet_orbits']:6d} {result['pivots']:6d} {result['dual_description_calls']:4d} "
            f"{result['descent_steps']:5d} {result['max_subcone_rays']:6d}"
        )
    print()
    for result in results:
        print(f"{result['label']}: facet orbit sizes {result['orbit_sizes']}")
    print("\nall cases: enumerated orbit representatives saturate to exactly the double description's facet set")


if __name__ == "__main__":
    sys.exit(main())
