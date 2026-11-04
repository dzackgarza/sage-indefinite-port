"""Stage-D attribution on provider e609d032 (hash fix 46b2c392), same published input.

Stages are run once each in this order, so the reduction and O(K) are warmed before the action; the stabilizer
call is the same call either way. Prints provenance so the imported provider tree is a fact, not an assumption.
"""
import faulthandler
import hashlib
import json
import os
import time

faulthandler.dump_traceback_later(int(os.environ.get("STACK_EVERY", "30")), repeat=True)
START = time.monotonic()


def stage(label):
    print(f"[t={time.monotonic() - START:8.1f}s] {label}", flush=True)


stage("imports")
import sage_extensions
from tests.acceptance.consumer import element, lattice, is_isometry, rows_of
from sage_indefinite_port.indefinite.recursive import EichlerRecursion

MODULE = sage_extensions.__file__
SOURCE = MODULE.rsplit("/sage_extensions/", 1)[0] + "/sage_extensions/categories/lattice_morphisms.py"
DIGEST = hashlib.sha256(open(SOURCE, "rb").read()).hexdigest()[:16]
stage(f"provider package: {MODULE}")
stage(f"lattice_morphisms.py sha256[:16] = {DIGEST}")

FIXTURE = json.load(open("/work/tests/fixtures/unpolarized_enriques.json"))
GRAM = FIXTURE["lattice"]["gram"]
N = lattice(GRAM)
REP = FIXTURE["component_preserving_group"]["line_representatives"][0]
v = element(N, REP["basis"][0])
stage(f"N: rank {N.rank()}, signature {N.signature_pair()}, det {FIXTURE['lattice']['det']}; line {REP['id']}")

REC = EichlerRecursion()
t0 = time.monotonic()
reduction = N.primitive_isotropic_subobject(v).isotropic_reduction()
stage(f"isotropic_reduction {time.monotonic() - t0:.1f}s")
quotient = reduction.quotient_lattice()
stage(f"K: rank {quotient.rank()}, signature {quotient.signature_pair()}, det {quotient.discriminant()}")

t0 = time.monotonic()
group_generators = REC.orthogonal_group_generators(quotient)
stage(f"orthogonal_group_generators(K) {time.monotonic() - t0:.1f}s: {len(tuple(group_generators))} gens")

t0 = time.monotonic()
action = REC.isotropic_stabilizer_action((v,), ((-v,),))
stage(f"isotropic_stabilizer_action {time.monotonic() - t0:.1f}s")

stage("action.stabilizer().group_generators()  <- the previously non-returning call")
t0 = time.monotonic()
stabilizer_generators = tuple(action.stabilizer().group_generators())
stage(f"action.stabilizer() {time.monotonic() - t0:.1f}s: {len(stabilizer_generators)} generators")

# the integral-generator obligation: genuine isometries of N, fixing the line setwise
integral = tuple(N.integral_isometry(generator) for generator in stabilizer_generators)
signs = []
for generator in integral:
    image = generator(v)
    if image == v:
        signs.append(1)
    elif image == -v:
        signs.append(-1)
    else:
        signs.append(None)
stage(f"line action per generator (1 = fixed, -1 = negated): {signs}")
stage(f"all preserve the form: {all(is_isometry(rows_of(generator), GRAM) for generator in integral)}")
stage(f"all fix the line setwise: {all(sign in (1, -1) for sign in signs)}")
print("ATTRIB_HASHFIX_OK", flush=True)
