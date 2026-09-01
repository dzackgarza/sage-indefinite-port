"""
TASK-03-2: Pairing Configuration Graphs and Sage Bliss/Nauty Canonical Labeling.

Encode finite vector/ray configurations into colored graphs, compute
canonical labelings via Sage Bliss/Nauty, and lift permutations to
integral lattice isometries.

Port of graph canonization from polyhedral_common.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

from sage.all import ZZ, QQ, matrix, Graph
from sage.graphs.graph import Graph as SageGraph


def pairing_configuration_graph(gram):
    """
    Build a colored graph encoding the pairing structure of a lattice.

    Each vertex represents a basis vector. Edge colors encode the
    value of b(e_i, e_j) = gram[i][j]. Vertex colors encode q(e_i) = gram[i][i].

    INPUT:
    - gram: nested list of integers (Gram matrix)

    OUTPUT:
    - Sage colored graph with vertex/edge labels
    """
    n = len(gram)
    G = SageGraph(multiedges=False)

    for i in range(n):
        G.add_vertex(i)

    # Color vertices by diagonal entry (norm)
    for i in range(n):
        G.set_vertex(i, gram[i][i])

    # Add edges with pairing values as colors
    for i in range(n):
        for j in range(i + 1, n):
            val = gram[i][j]
            if val != 0:
                G.add_edge(i, j, val)

    return G


def canonical_labeling(gram, backend="bliss"):
    """
    Compute canonical labeling of pairing configuration graph.

    INPUT:
    - gram: nested list of integers (Gram matrix)
    - backend: "bliss" or "nauty"

    OUTPUT:
    - permutation (list mapping old labels to new labels), canonical graph
    """
    G = pairing_configuration_graph(gram)

    if backend == "bliss":
        try:
            perm_dict, canon_g = G.canonical_label(
                certificate=True, edge_labels=True
            )
        except Exception:
            perm_dict, canon_g = G.canonical_label(certificate=True)
    else:
        perm_dict, canon_g = G.canonical_label(
            certificate=True, edge_labels=True
        )

    perm = [perm_dict[i] if i in perm_dict else i for i in range(G.order())]
    return perm, canon_g


def lift_automorphism(gram, perm):
    """
    Lift a graph automorphism (permutation of vertices) to an
    integral lattice isometry matrix.

    The permutation perm[i] = j means vertex i maps to vertex j.
    The matrix P has P[i, perm[i]] = 1, P[i, j] = 0 for j != perm[i].

    For this to be a lattice isometry, we need P^T Q P = Q.
    Returns None if the lift is not an isometry.

    INPUT:
    - gram: nested list of integers (Gram matrix)
    - perm: list mapping old vertex labels to new labels

    OUTPUT:
    - Sage integer matrix (the isometry), or None
    """
    n = len(gram)
    Q = matrix(ZZ, gram)
    P = matrix(ZZ, n, n)

    for i in range(n):
        P[i, perm[i]] = 1

    if P.transpose() * Q * P == Q:
        return P
    else:
        return None


def graph_automorphism_isometries(gram, backend="bliss"):
    """
    Compute all lattice isometries by lifting graph automorphisms.

    INPUT:
    - gram: nested list of integers (Gram matrix)
    - backend: "bliss" or "nauty"

    OUTPUT:
    - list of Sage integer matrices (verified isometries)
    """
    G = pairing_configuration_graph(gram)
    automs = G.automorphisms_group(
        edge_labels=True, representation="partition"
    )

    n = len(gram)
    isometries = []
    seen = set()

    for aut_partition in automs:
        # Convert partition to permutation
        perm = [0] * n
        for block in aut_partition:
            for i, v in enumerate(block):
                perm[v] = block[0] if i > 0 else block[(i + 1) % len(block)]

        # Try direct permutation
        M = lift_automorphism(gram, perm)
        if M is not None:
            key = tuple(M[i, j] for i in range(n) for j in range(n))
            if key not in seen:
                seen.add(key)
                isometries.append(M)

    return isometries
