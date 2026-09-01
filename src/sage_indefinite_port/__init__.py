"""sage_indefinite_port: Pure SageMath indefinite lattice algorithms.

Port of polyhedral_common/src_indefinite algorithms to exact SageMath
primitives with zero temporary-file or binary-executable dependencies.

Modules:
- lorentzian_perfect: Attack scheme, perfect domain traversal (Phase 04-1)
- approximate_models: 2U-Eichler models (Phase 04-2)
- higher_witt: Full O(L) assembly and equivalence (Phase 04-3)
- subgroup_constructors: Arithmetic subgroup types (Phase 04-4)
- isotropic_orbits: Parabolic recursion and orbit decomposition (Phase 05)
- equivariant: Centralizers and equivariant lattices (Phase 06-1)
- semantic_wrapper: Coordinate-free Sage API (Phase 06-2)
- validation: Literature milestone tests (Phase 06-3)
- graph_canonization: Bliss/Nauty canonical labeling (Phase 03-2)
- polyhedral_cones: Exact cone operations (Phase 03-3)
- definite_leaf: Definite lattice adapter (Phase 03-4)
- pre_sieves: Invariant pre-sieves (Phase 02-4)
- rational_group_integralization: Matrix group integralization (Phase 03-1)
"""

__version__ = "0.1.0"

# Public API
from .semantic_wrapper import (
    indefinite_lattice,
    orthogonal_generators,
    test_isometry,
    isotropic_orbits,
    IndefiniteLattice,
    OrthogonalGroup,
    IsotropicSubspace,
)
from .higher_witt import (
    indefinite_form_automorphism_group,
    indefinite_form_test_equivalence,
)
from .lorentzian_perfect import (
    get_attack_scheme,
    is_lorentzian,
    lorentzian_generators_autom,
    lorentzian_full_automorphism_group,
)
from .pre_sieves import (
    LatticeInvariant,
    VectorInvariant,
    lattice_pre_sieve,
    vector_pre_sieve,
)
from .equivariant import (
    EquivariantLattice,
    with_isometry,
    centralizer_group,
)
