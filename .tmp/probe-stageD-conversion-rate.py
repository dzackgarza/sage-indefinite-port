"""Measured processed count, evaluation and conversion rates on the family letters on N (I_1_1).

Prior phase isolated: atom A and dispatch stages are timed separately, then the FIRST BLOCK letters of
the resolution's generating family are read one at a time - each letter's word image and each image's
integral restriction are timed independently - and counted. Not the capped rerun: a bounded
rate-measurement probe; the tuple of the whole family is never materialized, so no 1383-1898s prior
phase is rebuilt. Rates are printed live every 50 letters; the derived full-family figure reports its
formula explicitly (processed count * measured rate), not as a measured total.
"""
import os
import time

START = time.monotonic()
STACK_EVERY = int(os.environ.get("STACK_EVERY", "120"))
import faulthandler

faulthandler.dump_traceback_later(STACK_EVERY, repeat=True)


def stage(label):
    print(f"[t={time.monotonic() - START:8.1f}s] {label}", flush=True)


stage("imports")
from tests.acceptance.consumer import element, lattice
from sage_indefinite_port.indefinite.recursive import EichlerRecursion
import json

FIXTURE = json.load(open("/work/tests/fixtures/unpolarized_enriques.json"))
GRAM = FIXTURE["lattice"]["gram"]
N = lattice(GRAM)
REP = FIXTURE["component_preserving_group"]["line_representatives"][0]
v = element(N, REP["basis"][0])
stage(f"N, line {REP['id']}")

REC = EichlerRecursion()
t0 = time.monotonic()
action = REC.isotropic_stabilizer_action((v,), ((-v,),))
stage(f"isotropic_stabilizer_action {time.monotonic() - t0:.1f}s")

t0 = time.monotonic()
rep = action.finite_representation()
stage(f"atomA finite_representation {time.monotonic() - t0:.1f}s (orbit {len(rep._submodules)})")

t0 = time.monotonic()
rho = rep.permutation_representation()
fstab = rep.finite_stabilizer()
stage(f"rho + finite_stabilizer {time.monotonic() - t0:.1f}s")

t0 = time.monotonic()
stab = action.stabilizer()
stage(f"stab construct {time.monotonic() - t0:.1f}s")

res = stab.selected_group_resolution()
labels = list(res.generating_set())
stage(f"family labels {len(labels)}")

BLOCK = int(os.environ.get("RATE_BLOCK", "300"))
eval_total = 0.0
conv_total = 0.0
seen = 0
stage(f"rate loop over the first {BLOCK} letters")
for label in labels[:BLOCK]:
    t1 = time.monotonic()
    g = res.generator(res.generating_set().retract(label) if False else label)
    eval_total += time.monotonic() - t1
    t1 = time.monotonic()
    g_int = N.integral_isometry(g)
    conv_total += time.monotonic() - t1
    seen += 1
    if seen % 50 == 0:
        stage(
            f"processed {seen}/{BLOCK}: eval {1000 * eval_total / seen:.1f} ms/lr, "
            f"conv {1000 * conv_total / seen:.1f} ms/lr"
        )

stage(
    f"MEASURED processed={seen} eval_total={eval_total:.2f}s ({1000 * eval_total / seen:.2f} ms/lr) "
    f"conv_total={conv_total:.2f}s ({1000 * conv_total / seen:.2f} ms/lr)"
)
total_letters = len(labels)
stage(
    "DERIVED full-family figures (formula: labels x measured rate; NOT a measured total): "
    f"group_generators_est={total_letters * eval_total / seen:.0f}s, conversions_est={total_letters * conv_total / seen:.0f}s"
)
stage("RATE_OK")
