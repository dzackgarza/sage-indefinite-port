"""Exact scalar conversion for the Allcock edgewalk algorithms."""

from sage.rings.rational import Rational
from sage.rings.rational_field import QQ as SageQQ


def _rational(entry) -> Rational:
    """Convert exact rational coordinates without cross-category Sage coercion."""
    numerator = getattr(entry, "numerator", entry)
    denominator = getattr(entry, "denominator", 1)
    numerator = numerator() if callable(numerator) else numerator
    denominator = denominator() if callable(denominator) else denominator
    return SageQQ(int(numerator)) / SageQQ(int(denominator))


def _integer(entry) -> int:
    r"""An integral coordinate; a rational entry with a denominator is a defect, not a value to truncate."""
    integer = int(entry)
    if entry != integer:
        raise ArithmeticError(f"the coordinate {entry} of a lattice vector is not integral")
    return integer
