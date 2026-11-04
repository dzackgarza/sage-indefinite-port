"""Bounded cost instrumentation of the post-repair stage-D stabilizer on N (I_1_1), PEAKRSS-guarded.

Separates the two cost atoms after the tensor-hash repair made the orbit honest:
  atom A: the finite representation's orbit BFS (finite_representation() -> from_finite_orbits)
  atom B: the Schreier generating family of the point stabilizer (stabilizer -> orbit_schreier_morphism)
A poller thread prints the orbit point count every 15 s so the growth rate and the per-atom
split are visible even at the cap. Not a stage-D rerun: capped diagnostics only.
"""
import os
import threading
import time

import faulthandler

START = time.monotonic()
STACK_EVERY = int(os.environ.get("STACK_EVERY", "120"))
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
stage(f"N from unpolarized_enriques.json: rank {N.rank()}, det {FIXTURE['lattice']['det']}; line {REP['id']}")

REC = EichlerRecursion()
stage("build isotropic_stabilizer_action (warmed like the stage-D probes)")
action = REC.isotropic_stabilizer_action((v,), ((-v,),))
stage("action built")

holder = {"rep": None}
stop = [False]


def poller():
    while not stop[0]:
        rep = holder.get("rep")
        if rep is not None:
            try:
                print(
                    f"[t={time.monotonic() - START:8.1f}s] POLL orbit points={len(rep._submodules)} positions={len(rep._positions)}",
                    flush=True,
                )
            except Exception as exc:
                print(f"POLL error: {exc!r}", flush=True)
        time.sleep(15)


threading.Thread(target=poller, daemon=True).start()

t0 = time.monotonic()
holder["rep"] = action.finite_representation()
t1 = time.monotonic()
rep = holder["rep"]
stage(f"ATOM A finite_representation() {t1 - t0:.1f}s: orbit points={len(rep._submodules)}")

t0 = time.monotonic()
stab = action.stabilizer()
t1 = time.monotonic()
sgens = tuple(stab.group_generators())
stage(f"ATOM B stabilizer() + group_generators() {t1 - t0:.1f}s: {len(sgens)} generators")

stop[0] = True
stage("COST_ATOMS_OK")
