"""Minimal path for the slow generator enumeration in Enriques case 1.

Measured phase: `OwnedGroups.ParentMethods.group_generators`, 21 calls totalling
177.52 s with a worst single call of 118.67 s, on a run whose other three phases are
negligible (subgroup admission 1 call / 0.03 s, 957 cardinality calls / 0.01 s,
21 discriminant morphisms / 13.75 s).

The exact parent, from source rather than from the test's outer line:

  lattice_centralizers.py:852  coinvariant_extension_subgroup
      stabilizer = extension.invariant.O().stabilizer(invariant_polarization)
  lattice_centralizers.py:985  the retained generating morphism
      invariant_image = invariant_subgroup.image_under(...)      <- needs generators
      ... image_under -> group_generators -> selected_group_resolution
      -> realize -> lattice_morphisms.vector_stabilizer_generators      (:3664)
      -> _marking_stabilizer_generators                                (:3691)
      -> lattices._marked_automorphism_generators                     (:1874)
      -> port marked_lattices.automorphism_group_generators            (:109)
      -> port recursive.pointed_automorphism_generators                (:591)
      -> recursive._lattice_action (:457) -> orthogonal_group_generators (:333)
      -> lattices._core_automorphism_generators                        (:1849)
      -> lattice_morphisms._computed_group_generators (:3575) / _engine_group (:3456)
      -> lattice_engines.definite_isometry_group                       (:1146) -> Julia

So the parent is the stabilizer of the case-1 polarization inside O(M), M = U(2) + E8(-2).
Its retained presentation is a marked-lattice stabilizer subgroup, realised through the
port's pointed-automorphism route; `eichler.OrbitCoverModel` is not involved, which is
consistent with the separate measurement of zero OrbitCoverModel calls.

Repeat requests are expected to be free: `selected_resolution` (objects.py:383) replaces a
`_PendingResolution` with the realised resolution in the per-owner registry, and
`group_generators` itself is not cached. This script measures that directly instead of
assuming it: realise once, then request again.

Run after the sole bounded probe has finished. Not run yet.
"""
import time

from sage.all import QQ

from sage_extensions.catalogue import Involutions

from tests.acceptance.consumer import element, lattice
from tests.fixtures.oracle_fixtures import (
    load_e10_fundamental_domain,
    load_enriques_87_polarizations,
)

case = next(row for row in load_enriques_87_polarizations() if row["case"] == 1)
print("case 1 polarization:", case["face_polarization"], "order:", case["group_order"])

involution = Involutions.I_En
invariant = involution.primitive_extension().invariant
transport = lattice(load_e10_fundamental_domain()["gram"]).twist(QQ(2)).Isom(invariant).an_element()
polarization = invariant.inclusion()(transport(element(transport.domain(), case["face_polarization"])))

stabilizer = invariant.O().stabilizer(polarization)
print("parent:", stabilizer)

started = time.monotonic()
first = stabilizer.group_generators()
first_seconds = time.monotonic() - started

started = time.monotonic()
second = stabilizer.group_generators()
second_seconds = time.monotonic() - started

print(f"initial realisation : {first_seconds:.2f}s")
print(f"repeat request      : {second_seconds:.4f}s")
print(f"ratio               : {first_seconds / max(second_seconds, 1e-9):.0f}x")
print(f"same object returned: {first is second}")
