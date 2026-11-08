"""Attribute the non-returning case_3 stage-D call by sub-stage, on the published input.

Same input as the retained repro (.tmp/repro_marked_stabilizer.py at sage-extensions 0dad201c): the direct-gram
N from tests/fixtures/unpolarized_enriques.json and the published line representative I_1_1 (Z e1). It prints
before/after each sub-stage of IsotropicInduction._stabilizer_generators so the non-returned call is named by
markers, and faulthandler dumps the live stack so its frame is evidence, not inference.
"""
import faulthandler
import json
import os
import time

EVERY = int(os.environ.get("STACK_EVERY", "20"))
faulthandler.dump_traceback_later(EVERY, repeat=True)
START = time.monotonic()


def stage(label):
    print(f"[t={time.monotonic() - START:8.1f}s] {label}", flush=True)


stage("imports")
from tests.acceptance.consumer import element, lattice
from sage_indefinite_port.indefinite.recursive import EichlerRecursion

FIXTURE = json.load(open("/work/tests/fixtures/unpolarized_enriques.json"))
N = lattice(FIXTURE["lattice"]["gram"])
v = element(N, FIXTURE["component_preserving_group"]["line_representatives"][0]["basis"][0])
stage(f"N built: rank {N.rank()}, signature {N.signature_pair()}, det {FIXTURE['lattice']['det']}; v = {FIXTURE['component_preserving_group']['line_representatives'][0]['id']}")

REC = EichlerRecursion()

t0 = time.monotonic()
reduction = N.primitive_isotropic_subobject(v).isotropic_reduction()
stage(f"isotropic_reduction done in {time.monotonic() - t0:.1f}s")

quotient = reduction.quotient_lattice()
stage(f"quotient K: rank {quotient.rank()}, signature {quotient.signature_pair()}, det {quotient.discriminant()}")

from sage_indefinite_port.invariants import AttackProfile
profile = AttackProfile.from_lattice(quotient)
stage(f"AttackProfile: positive_index {profile.positive_index}, negative_index {getattr(profile, 'negative_index', '?')}")

t0 = time.monotonic()
generators = REC.orthogonal_group_generators(quotient)
stage(f"orthogonal_group_generators(K) done in {time.monotonic() - t0:.1f}s: {len(tuple(generators))} generators")

t0 = time.monotonic()
kernel = tuple(reduction.pointwise_perpendicular_kernel().group_generators())
stage(f"pointwise_perpendicular_kernel done in {time.monotonic() - t0:.1f}s: {len(kernel)} generators")

stage("full isotropic_stabilizer_action")
t0 = time.monotonic()
action = REC.isotropic_stabilizer_action((v,), ((-v,),))
stage(f"isotropic_stabilizer_action done in {time.monotonic() - t0:.1f}s")

t0 = time.monotonic()
stab = action.stabilizer().group_generators()
stage(f"action.stabilizer() done in {time.monotonic() - t0:.1f}s: {len(tuple(stab))} generators")
print("ATTRIB_OK", flush=True)
