"""Typed exception hierarchy for indefinite lattice algorithms."""

from __future__ import annotations


class IndefinitePortError(Exception):
    """Base exception for sage-indefinite-port."""


class NotPrimitiveError(IndefinitePortError, ValueError):
    """Raised when a vector or sublattice is required to be primitive but is not."""


class NotIsotropicError(IndefinitePortError, ValueError):
    """Raised when an object is required to be isotropic but is not."""


class NotEvenError(IndefinitePortError, ValueError):
    """Raised when a form or lattice is required to be even but is not."""


class InfiniteLocusError(IndefinitePortError, ValueError):
    """Raised when a requested vector locus is infinite (e.g. non-primitive norm 0)."""


class NoIntegralExtensionError(IndefinitePortError, ArithmeticError):
    """Raised when a rational isometry has an empty integral extension locus."""


class UnsupportedRankError(IndefinitePortError, ValueError):
    """Raised when an isotropic rank or lattice dimension exceeds current capabilities."""


class NoLocalMarkTheoremError(IndefinitePortError, RuntimeError):
    """Raised when Lorentzian marks cannot be resolved locally."""
