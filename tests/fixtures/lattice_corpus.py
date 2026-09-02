"""
Test lattice corpus for differential verification against upstream polyhedral_common.

Each entry provides:
- name: human-readable identifier
- gram: Gram matrix as nested list (symmetric, integral)
- signature: (p, q, z) = (positive, negative, zero eigenvalues)
- det: determinant
- is_even: whether diagonal entries are all even
- isometry_class: known class invariant
- isotropic_vector_count: number of primitive isotropic orbit reps (if known)
- automorphism_group_order: |O(Q,Z)| for finite cases, None for infinite
- notes: source reference
"""

from typing import TypedDict


class LatticeEntry(TypedDict):
    name: str
    gram: list[list[int]]
    signature: tuple[int, int, int]
    det: int
    is_even: bool
    isometry_class: str
    isotropic_vector_count: int | None
    automorphism_group_order: int | None
    notes: str


CORPUS: list[LatticeEntry] = [
    # === Hyperbolic plane U ===
    {
        "name": "hyperbolic_plane_U",
        "gram": [[0, 1], [1, 0]],
        "signature": (1, 1, 0),
        "det": -1,
        "is_even": True,
        "isometry_class": "U",
        "isotropic_vector_count": 1,
        "automorphism_group_order": None,
        "notes": "The standard hyperbolic plane. Infinite orthogonal group.",
    },
    # === Hyperbolic plane U(2) ===
    {
        "name": "hyperbolic_plane_U2",
        "gram": [[0, 2], [2, 0]],
        "signature": (1, 1, 0),
        "det": -4,
        "is_even": True,
        "isometry_class": "U(2)",
        "isotropic_vector_count": 1,
        "automorphism_group_order": None,
        "notes": "Scaled hyperbolic plane.",
    },
    # === SPLAG 51a: det -128, two spinor genera ===
    {
        "name": "SPLAG_51a",
        "gram": [[-1, 1, 0], [1, 63, 0], [0, 0, 2]],
        "signature": (2, 1, 0),
        "det": -128,
        "is_even": False,
        "isometry_class": "51a",
        "isotropic_vector_count": None,
        "automorphism_group_order": None,
        "notes": "SPLAG Table 15.1, genus 51a (forms with same discriminant and local invariants, different spinor genus).",
    },
    # === SPLAG 51b: det -128, other genus ===
    {
        "name": "SPLAG_51b",
        "gram": [[-9, 1, 0], [1, 7, 0], [0, 0, 2]],
        "signature": (2, 1, 0),
        "det": -128,
        "is_even": False,
        "isometry_class": "51b",
        "isotropic_vector_count": None,
        "automorphism_group_order": None,
        "notes": "SPLAG Table 15.1, genus 51b. In same genus as 51a, not in same spinor genus.",
    },
    # === Lorentzian unimodular II_{1,9} ===
    {
        "name": "II_1_9",
        "gram": [
            # U + E8(-1)
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, -2, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, -2, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, -2, 1, 0, 0, 0, 1],
            [0, 0, 0, 0, 1, -2, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 1, -2, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 1, -2, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, -2, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, -2],
        ],
        "signature": (1, 9, 0),
        "det": -1,
        "is_even": True,
        "isometry_class": "II_{1,9}",
        "isotropic_vector_count": 1,
        "automorphism_group_order": None,
        "notes": "Unimodular Lorentzian lattice of rank 10. Automorphism group is W(E10).",
    },
    # === Nikulin Enriques lattice N = U + U(2) + E8(-2) ===
    {
        "name": "Enriques_N",
        "gram": [
            # U: 2x2
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            # U(2): 2x2
            [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            # E8(-2): 8x8
            [0, 0, 0, 0, -4, 2, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 2, -4, 2, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 2, -4, 2, 0, 0, 0, 2],
            [0, 0, 0, 0, 0, 0, 2, -4, 2, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 2, -4, 2, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 2, -4, 2, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -4, 0],
            [0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, -4],
        ],
        "signature": (2, 10, 0),
        "det": -(2**10),
        "is_even": True,
        "isometry_class": "U + U(2) + E8(-2)",
        "isotropic_vector_count": 2,
        "automorphism_group_order": None,
        "notes": "Enriques surface second cohomology lattice (Nikulin 1984, Sterk 1985). Two cusp orbits.",
    },
    # === Definite: A2 ===
    {
        "name": "A2",
        "gram": [[2, -1], [-1, 2]],
        "signature": (2, 0, 0),
        "det": 3,
        "is_even": True,
        "isometry_class": "A2",
        "isotropic_vector_count": 0,
        "automorphism_group_order": 12,
        "notes": "Positive definite A2 root lattice. Automorphism group is W(A2) x {+-1} = D6 x Z2 = 12.",
    },
    # === Definite: D4 ===
    {
        "name": "D4",
        "gram": [
            [2, -1, 0, 0],
            [-1, 2, -1, -1],
            [0, -1, 2, 0],
            [0, -1, 0, 2],
        ],
        "signature": (4, 0, 0),
        "det": 4,
        "is_even": True,
        "isometry_class": "D4",
        "isotropic_vector_count": 0,
        "automorphism_group_order": 1152,
        "notes": "Positive definite D4 root lattice. Triality: Aut(D4) = W(D4) : S3, order 1152.",
    },
    # === Definite: E8 ===
    {
        "name": "E8",
        "gram": [
            [2, -1, 0, 0, 0, 0, 0, 0],
            [-1, 2, -1, 0, 0, 0, 0, 0],
            [0, -1, 2, -1, 0, 0, 0, -1],
            [0, 0, -1, 2, -1, 0, 0, 0],
            [0, 0, 0, -1, 2, -1, 0, 0],
            [0, 0, 0, 0, -1, 2, -1, 0],
            [0, 0, 0, 0, 0, -1, 2, 0],
            [0, 0, -1, 0, 0, 0, 0, 2],
        ],
        "signature": (8, 0, 0),
        "det": 1,
        "is_even": True,
        "isometry_class": "E8",
        "isotropic_vector_count": 0,
        "automorphism_group_order": 696729600,
        "notes": "Positive definite E8 lattice. Reference for automorphism group size.",
    },
]


def get_gram(name: str) -> list[list[int]]:
    """Return the Gram matrix for a named lattice."""
    for entry in CORPUS:
        if entry["name"] == name:
            return entry["gram"]
    raise KeyError(f"Unknown lattice: {name}")


def get_all_indefinite() -> list[LatticeEntry]:
    """Return only indefinite lattices (signature with both p>0 and q>0)."""
    return [e for e in CORPUS if e["signature"][0] > 0 and e["signature"][1] > 0]


def get_all_definite() -> list[LatticeEntry]:
    """Return only definite lattices."""
    return [e for e in CORPUS if e["signature"][1] == 0 or e["signature"][0] == 0]
