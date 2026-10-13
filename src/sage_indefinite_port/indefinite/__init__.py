"""Indefinite-lattice reduction and lifting algorithms."""

from sage_indefinite_port.indefinite.eichler import (
    EichlerEnvelope,
    EichlerOrbitCover,
    InfiniteLocusError,
    OrbitCover,
    OrbitCoverModel,
    TwoHyperbolicPlaneDecomposition,
    eichler_transvection,
    find_hyperbolic_pair,
    square_divisors,
)
from sage_indefinite_port.indefinite.isotropic_lifts import (
    CodimensionOneIsotropicExtension,
    CodimensionOneIsotropicExtensionResult,
    IntegralParameterCoset,
    IsometryExtensionTorsor,
    MatrixEquationSolution,
    NoIntegralExtensionError,
    PointwisePerpendicularKernel,
    pointwise_perpendicular_kernel,
    solve_isotropic_extension_equation,
)
from sage_indefinite_port.indefinite.isotropic_reductions import FlagType
from sage_indefinite_port.indefinite.presentations import ReducedLatticePresentation, presentation_bucket_key
from sage_indefinite_port.indefinite.vector_sections import (
    IsotropicVectorSection,
    NonIsotropicVectorSection,
    VectorOrthogonalSection,
    orthogonal_section,
)

__all__ = [
    "CodimensionOneIsotropicExtension",
    "EichlerEnvelope",
    "EichlerOrbitCover",
    "eichler_transvection",
    "find_hyperbolic_pair",
    "InfiniteLocusError",
    "FlagType",
    "CodimensionOneIsotropicExtensionResult",
    "IntegralParameterCoset",
    "IsometryExtensionTorsor",
    "IsotropicVectorSection",
    "MatrixEquationSolution",
    "NoIntegralExtensionError",
    "NonIsotropicVectorSection",
    "OrbitCover",
    "OrbitCoverModel",
    "PointwisePerpendicularKernel",
    "ReducedLatticePresentation",
    "VectorOrthogonalSection",
    "orthogonal_section",
    "pointwise_perpendicular_kernel",
    "solve_isotropic_extension_equation",
    "square_divisors",
    "TwoHyperbolicPlaneDecomposition",
    "presentation_bucket_key",
]
