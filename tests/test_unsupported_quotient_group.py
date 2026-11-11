r"""A null line whose stabilizer needs a quotient group without a backend.

L = U + K with K = diag(2, 2, -14, -14). The form x^2 + y^2 - 7 z^2 - 7 w^2 is
anisotropic over Q: a sum of two squares has even 7-adic valuation, and 7 times
a nonzero sum of two squares has odd valuation. So K has signature (2, 2) and
rational Witt index zero. K is not isometric to the supported quaternion model
2(x^2 + y^2 - 3 z^2 - 3 w^2), whose determinant is 144, not 784.

The line Ze of U exists, and its stabilizer is built from O(e^perp/e) = O(K).
The pinned approximate model needs two hyperbolic planes and stops for K at
IndefiniteForms.g:830 ("Something more to write"). The stabilizer request must
fail with the precise unsupported-algorithm error, not an incidental one.
"""
import pytest

from tests.acceptance.consumer import lattice


def test_line_stabilizer_with_an_anisotropic_rank_four_quotient_is_unsupported():
    gram = [[0, 1, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0], [0, 0, 2, 0, 0, 0], [0, 0, 0, 2, 0, 0], [0, 0, 0, 0, -14, 0], [0, 0, 0, 0, 0, -14]]
    ambient = lattice(gram)
    e = next(iter(ambient.module_generators()))
    line = ambient.primitive_isotropic_subobject(e)

    with pytest.raises(NotImplementedError):
        tuple(ambient.O().isotropic_stabilizer_generators(line))
