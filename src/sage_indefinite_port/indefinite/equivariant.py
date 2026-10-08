r"""Centralizers of isometries of indefinite lattices (map section 15).

The research preamble's centralizer subgroup ``G.centralizer(f)`` of ``G = L.O()`` takes its
generators from here for an indefinite lattice ``L``. They come from the eigensublattices
``L_+`` and ``L_-`` of an involution and the stabilizer of the gluing datum
``L/(L_+ + L_-)`` in ``O(L_+) x O(L_-)``, extended to ``L``. The stabilizer itself is a
preamble operation.
"""

from dzack_research.preamble.categories.lattice_morphisms import LatticeIsometryMethods, LatticeIsometryMor
from sage.structure.parent import Parent

from sage_indefinite_port.readiness import unfinished


@unfinished(23)
def centralizer_generators(group: LatticeIsometryMor, isometry: LatticeIsometryMethods) -> Parent:
    r"""Generators of the centralizer of ``isometry`` in ``G``."""


__all__ = ["centralizer_generators"]
