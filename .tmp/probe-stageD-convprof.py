"""cProfile over the measured conversion feedback on N (I_1_1, first BLOCK letters).

Same bounded specimen as the rate probe (atom A + dispatch, then per-letter evaluation and
integral restriction), profiling the whole loop, so the cumulative cost inside one conversion is
named by its own frames and never by attribution. No family tuple is materialized.
"""
import cProfile
import os
import pstats
import time

PROFILE_BLOCK = int(os.environ.get("RATE_BLOCK", "60"))


def body():
    t0 = time.monotonic()
    from tests.acceptance.consumer import element, lattice
    from sage_indefinite_port.indefinite.recursive import EichlerRecursion
    import json

    FIXTURE = json.load(open("/work/tests/fixtures/unpolarized_enriques.json"))
    N = lattice(FIXTURE["lattice"]["gram"])
    REP = FIXTURE["component_preserving_group"]["line_representatives"][0]
    v = element(N, REP["basis"][0])
    REC = EichlerRecursion()
    action = REC.isotropic_stabilizer_action((v,), ((-v,),))
    print(f"[t={time.monotonic() - t0:.1f}] action", flush=True)
    rep = action.finite_representation()
    print(f"[t={time.monotonic() - t0:.1f}] atomA", flush=True)
    rep.permutation_representation()
    rep.finite_stabilizer()
    stab = action.stabilizer()
    print(f"[t={time.monotonic() - t0:.1f}] stab construct", flush=True)
    res = stab.selected_group_resolution()
    labels = list(res.generating_set())
    print(f"[t={time.monotonic() - t0:.1f}] labels {len(labels)}", flush=True)
    eval_total = 0.0
    conv_total = 0.0
    for label in labels[:PROFILE_BLOCK]:
        t1 = time.monotonic()
        g = res.generator(label)
        eval_total += time.monotonic() - t1
        t1 = time.monotonic()
        N.integral_isometry(g)
        conv_total += time.monotonic() - t1
    print(
        f"[t={time.monotonic() - t0:.1f}] MEASURED n={PROFILE_BLOCK} "
        f"eval={1000 * eval_total / PROFILE_BLOCK:.1f}ms conv={1000 * conv_total / PROFILE_BLOCK:.1f}ms",
        flush=True,
    )


prof = cProfile.Profile()
prof.enable()
body()
prof.disable()
stats = pstats.Stats(prof)
stats.sort_stats("cumulative").print_stats(24)
print("CONVPROF_OK", flush=True)
