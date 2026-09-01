"""
Monkey-patch py_polyhedral.binaries to replace binary shell-outs
with pure SageMath implementations.

This is the integration seam: import this module before using
any preamble lattice API on indefinite lattices.
"""

import sys
import os
sys.path.insert(0, os.path.expanduser("~/research/src"))

import py_polyhedral.binaries as _binaries


def _patch():
    """Replace all indefinite binary functions with pure SageMath versions."""
    from sage_indefinite_port.indefinite_algorithms import (
        indefinite_form_automorphism_group,
        indefinite_form_test_equivalence,
        indefinite_form_test_equivalence_vector,
        indefinite_form_get_orbit_representative,
        indefinite_form_isotropic_k_plane,
        indefinite_form_isotropic_k_flag,
        indefinite_form_stabilizer_vector,
        indefinite_form_stabilizer_isotropic_subspace,
    )

    _binaries.indefinite_form_automorphism_group = indefinite_form_automorphism_group
    _binaries.indefinite_form_test_equivalence = indefinite_form_test_equivalence
    _binaries.indefinite_form_test_equivalence_vector = indefinite_form_test_equivalence_vector
    _binaries.indefinite_form_get_orbit_representative = indefinite_form_get_orbit_representative
    _binaries.indefinite_form_isotropic_k_plane = indefinite_form_isotropic_k_plane
    _binaries.indefinite_form_isotropic_k_flag = indefinite_form_isotropic_k_flag
    _binaries.indefinite_form_stabilizer_vector = indefinite_form_stabilizer_vector
    _binaries.indefinite_form_stabilizer_isotropic_subspace = indefinite_form_stabilizer_isotropic_subspace


_patch()
