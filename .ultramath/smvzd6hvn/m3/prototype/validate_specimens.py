"""Validate the orbit-wise dual description against exact full facet enumerations.

Specimens are perfect-cell-shaped: a Lorentzian lattice of Gram
``I_k (+) <-1>`` (rank ``d = k + 1``) and the ray set
``{(p, 1) : p in Z^k, |p|^2 = R}``, i.e. the lattice vectors of norm ``R - 1``
lying in the affine hyperplane ``x_d = 1``.  The stabiliser generators are the
integral isometries induced by signed coordinate permutations of the first
``k`` coordinates; they permute the ray set exactly.

For each specimen the script checks, against PPL's complete exact
H-description of the same cone:

  * every facet found by PPL lies in exactly one enumerated orbit;
  * the union of the orbits is exactly the PPL facet set (no extras, none missing);
  * every representative normal is an exact supporting covector whose zero set
    is the representative incidence and whose incidence has rank d-1;
  * every orbit member's Schreier word carries the representative onto it, and
    the transported representative normal is, up to a positive rational scale,
    the member's own primitive normal; the transported isometry is an integral
    isometry of the Gram matrix with determinant +-1.
"""

from __future__ import annotations

import time
from itertools import product

from sage.all import QQ, ZZ, matrix, vector

from orbit_addm import (
    Stats,
    apply_word,
    facet_normal,
    facet_orbits,
    full_facet_incidences,
    incidence,
    pair,
    rank_of,
    same_ray,
)


def shell_points(k, square):
    """All ``p in Z^k`` with ``|p|^2 = square``, sorted."""
    points = []

    def recurse(prefix, remaining, accumulated):
        if len(prefix) == k:
            if accumulated == square:
                points.append(tuple(prefix))
            return
        left = k - len(prefix) - 1
        for value in range(-square, square + 1):
            rest = accumulated + value * value
            if rest > square:
                continue
            if rest + left * square < square:
                continue
            prefix.append(value)
            recurse(prefix, remaining, rest)
            prefix.pop()

    recurse([], k, 0)
    return sorted(set(points))


def gram_of(k):
    gram = matrix(ZZ, k + 1, k + 1, lambda i, j: 0)
    for index in range(k):
        gram[index, index] = 1
    gram[k, k] = -1
    return gram


def generators_of(rows, k, *, sign_flips=True):
    """Position permutations induced by signed coordinate permutations."""
    dimension = k + 1
    positions = {row: index for index, row in enumerate(rows)}
    actions = []
    for index in range(k - 1):
        action = matrix(ZZ, dimension, dimension, lambda i, j: 0)
        for label in range(dimension):
            action[label, label] = 1
        action[index, index] = 0
        action[index + 1, index + 1] = 0
        action[index + 1, index] = 1
        action[index, index + 1] = 1
        actions.append(action)
    if sign_flips:
        for index in range(k):
            action = matrix(ZZ, dimension, dimension, lambda i, j: 0)
            for label in range(dimension):
                action[label, label] = 1
            action[index, index] = -1
            actions.append(action)
    permutations = []
    for action in actions:
        images = []
        for row in rows:
            image = tuple(int(x) for x in (vector(ZZ, row) * action))
            images.append(positions[image])
        permutations.append(tuple(images))
    return permutations, actions


def basis_positions(rows):
    return tuple(int(index) for index in matrix(QQ, rows).transpose().pivots())


def transport_normal(rows, generator_permutations, word, normal):
    """The row action recovered from a Schreier word, and the transported covector."""
    positions = basis_positions(rows)
    source_basis = matrix(QQ, [rows[index] for index in positions])
    image_positions = list(positions)
    for generator_index in word:
        permutation = generator_permutations[generator_index]
        image_positions = [permutation[position] for position in image_positions]
    action = source_basis.inverse() * matrix(QQ, [rows[index] for index in image_positions])
    # Rows act by x -> x . action, so a covector u acts as u . (action^{-1})^T:
    # that is the covector of the same name in the transported pairings.
    transported = vector(QQ, normal) * action.inverse().transpose()
    return action, tuple(transported)


def run_specimen(name, k, square, *, sign_flips=True, verbose=True):
    points = shell_points(k, square)
    rows = tuple(tuple(list(point) + [1]) for point in points)
    dimension = k + 1
    assert len(set(rows)) == len(rows)
    assert rank_of(rows) == dimension, (name, rank_of(rows), dimension)
    height = tuple([0] * k + [1])
    assert all(pair(height, row) == 1 for row in rows)
    gram = gram_of(k)
    generators, actions = generators_of(rows, k, sign_flips=sign_flips)

    started = time.perf_counter()
    stats = Stats()
    orbits = facet_orbits(rows, generators, dimension, height=height, stats=stats)
    elapsed = time.perf_counter() - started

    reference_started = time.perf_counter()
    reference = full_facet_incidences(rows)
    reference_elapsed = time.perf_counter() - reference_started
    reference_set = set(reference)
    assert len(reference_set) == len(reference), "PPL reported a repeated incidence"

    members: dict = {}
    for index, record in enumerate(orbits):
        for member in record["members"]:
            assert member not in members, "an enumerated facet appears in two orbits"
            members[member] = index
    assert set(members) == reference_set, (
        f"orbit union != facet set: missing {len(reference_set - set(members))}, "
        f"extra {len(set(members) - reference_set)}"
    )

    # exact normals of the representatives
    for index, record in enumerate(orbits):
        normal = record["normal"]
        values = [pair(normal, row) for row in rows]
        assert all(value >= 0 for value in values)
        assert incidence(rows, normal) == record["face"]
        assert rank_of([rows[position] for position in sorted(record["face"])]) == dimension - 1
        assert facet_normal(rows, record["face"], dimension) == normal, (
            "the representative normal is not the canonical one"
        )

    # exact transport: word -> isometry -> covector, for every member
    for index, record in enumerate(orbits):
        representative = record["face"]
        normal = record["normal"]
        for member, word in record["members"].items():
            assert apply_word(representative, generators, word) == member
            action, transported = transport_normal(rows, generators, word, normal)
            assert all(entry.denominator() == 1 for entry in action.list()), "nonintegral action"
            action = action.change_ring(ZZ)
            assert abs(int(action.det())) == 1
            assert action * gram * action.transpose() == gram, "the word is not a lattice isometry"
            assert all(entry.denominator() == 1 for entry in transported)
            transported = tuple(int(entry) for entry in transported)
            assert incidence(rows, transported) == member, "the transported normal does not cut out the member"
            member_normal = facet_normal(rows, member, dimension)
            assert same_ray(transported, member_normal), "transported normal is not the member's normal"

    orbit_sizes = sorted(len(record["members"]) for record in orbits)
    if verbose:
        print(
            f"{name:<28} rays={len(rows):>5} dim={dimension} facets={len(reference_set):>4} "
            f"orbits={len(orbits):>3} sizes={orbit_sizes if len(orbit_sizes) <= 12 else str(orbit_sizes[:6]) + '...'} "
            f"gens={len(generators)} | {stats} | adm={elapsed:.3f}s full={reference_elapsed:.3f}s"
        )
    return {
        "name": name,
        "rays": len(rows),
        "dimension": dimension,
        "facets": len(reference_set),
        "orbits": len(orbits),
        "orbit_sizes": orbit_sizes,
        "stats": stats,
        "elapsed": elapsed,
        "reference_elapsed": reference_elapsed,
    }


SPECIMENS = [
    ("Z2 shell1 (isotropic)", 2, 1, True),
    ("Z2 shell2", 2, 2, True),
    ("Z2 shell4", 2, 4, True),
    ("Z3 shell1 (isotropic)", 3, 1, True),
    ("Z3 shell2", 3, 2, True),
    ("Z3 shell4", 3, 4, True),
    ("Z4 shell1 (isotropic)", 4, 1, True),
    ("Z4 shell2", 4, 2, True),
    ("Z4 shell4", 4, 4, True),
    ("Z4 shell4 no flips", 4, 4, False),
    ("Z5 shell2", 5, 2, True),
    ("Z5 shell4", 5, 4, True),
    ("Z5 shell6", 5, 6, True),
    ("Z6 shell2", 6, 2, True),
    ("Z6 shell4", 6, 4, True),
]


def main():
    results = []
    for name, k, square, flips in SPECIMENS:
        results.append(run_specimen(name, k, square, sign_flips=flips))
    # trivial group: the enumeration must produce every facet as its own orbit
    for name, k, square in [("Z3 shell2 trivial", 3, 2), ("Z4 shell2 trivial", 4, 2)]:
        results.append(run_specimen(name, k, square, sign_flips=False, verbose=True))
        rows = tuple(tuple(list(point) + [1]) for point in shell_points(k, square))
        orbits = facet_orbits(rows, (), k + 1, height=tuple([0] * k + [1]))
        assert len(orbits) == len(full_facet_incidences(rows)), "trivial group must give one orbit per facet"
        assert all(len(record["members"]) == 1 for record in orbits)
    print("\nall specimens agree with the exact full dual description")


if __name__ == "__main__":
    main()
