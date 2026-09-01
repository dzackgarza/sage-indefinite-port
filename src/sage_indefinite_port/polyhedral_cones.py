"""
TASK-03-3: Exact Rational Polyhedral Cone Operations.

Wraps Sage polyhedron operations (which delegate to cddlib/PPL/Normaliz)
for exact rational cone computations: extreme rays, facet enumeration,
incidence, and face stabilizers.

Port of polyhedral cone operations from polyhedral_common.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import ZZ, QQ, matrix, vector, Polyhedron


def cone_from_rays(ray_matrix):
    """
    Construct a rational polyhedral cone from a matrix of generating rays.

    INPUT:
    - ray_matrix: nested list of rationals/integers (rows are rays)

    OUTPUT:
    - Sage Polyhedron object
    """
    rays = [vector(QQ, row) for row in ray_matrix]
    return Polyhedron(rays=rays, base_ring=QQ)


def extreme_rays(cone):
    """
    Compute extreme rays of a polyhedral cone.

    INPUT:
    - cone: Sage Polyhedron

    OUTPUT:
    - list of vectors (extreme rays)
    """
    return cone.rays()


def facet_inequalities(cone):
    """
    Compute facet-defining inequalities of a cone.

    INPUT:
    - cone: Sage Polyhedron

    OUTPUT:
    - list of (normal_vector, rhs) pairs
    """
    H = cone.Hrepresentation()
    ineqs = []
    for h in H:
        if h.is_inequality():
            coeffs = vector(ZZ, h.A())
            ineqs.append((coeffs, ZZ(h.b())))
    return ineqs


def cone_incidence(cone):
    """
    Compute ray-facet incidence matrix.

    INPUT:
    - cone: Sage Polyhedron

    OUTPUT:
    - nested list (1 if ray on facet, 0 otherwise)
    """
    inc = cone.incidence_matrix()
    n_rays = len(cone.rays())
    n_facets = len(cone.Hrepresentation())
    return [[int(inc[i, j]) for j in range(n_facets)] for i in range(n_rays)]


def face_from_rays(cone, ray_indices):
    """
    Compute the face of a cone spanned by given rays.

    INPUT:
    - cone: Sage Polyhedron
    - ray_indices: list of ray indices

    OUTPUT:
    - Sage Polyhedron (the face)
    """
    rays = cone.rays()
    face_rays = [rays[i] for i in ray_indices]
    return Polyhedron(rays=face_rays, base_ring=QQ)


def cone_intersection(cone1, cone2):
    """
    Compute intersection of two cones.

    INPUT:
    - cone1, cone2: Sage Polyhedra

    OUTPUT:
    - Sage Polyhedron
    """
    return cone1.intersection(cone2)


def cone_stabilizer_generators(cone, isometry_gens):
    """
    Compute generators of the stabilizer of a cone under a group of isometries.

    INPUT:
    - cone: Sage Polyhedron
    - isometry_gens: list of Sage matrices (isometries)

    OUTPUT:
    - list of isometry generators preserving the cone
    """
    ray_set = set(tuple(r) for r in cone.rays())
    stabilizer = []

    for M in isometry_gens:
        preserves = True
        for r in cone.rays():
            img = tuple(M * r)
            # Check if image is in cone
            if not cone.contains(Polyhedron(vertices=[img])):
                preserves = False
                break
        if preserves:
            stabilizer.append(M)

    return stabilizer
