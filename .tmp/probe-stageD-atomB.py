"""Disaggregate atom B of the post-repair stage-D stabilizer on N (I_1_1), PEAKRSS-guarded.

After ATOM A (finite representation: 154.8 s, orbit 1024) finished, the stabilizer stages did not
end at cap. Time each stage of lattice_stabilizer separately so the wall has an owner:
  B1 permutation_representation() (rho: G -> Sym(orbit))
  B2 finite_stabilizer()          (point stabilizer in the permutation codomain)
  B3 orbit_schreier_morphism      (the words tuple: |orbit| x |letters|)
  B4 preimage_subgroup + group_generators() (materialization of the Schreier family)
For B3, first time 200 (point, generator) iterations of the exact same comprehension and
extrapolate before paying the whole tuple.
"""
import os
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
stage(f"N from unpolarized_enriques.json, line {REP['id']}")

REC = EichlerRecursion()
action = REC.isotropic_stabilizer_action((v,), ((-v,),))
t0 = time.monotonic()
rep = action.finite_representation()
stage(f"ATOM A finite_representation() {time.monotonic() - t0:.1f}s: points={len(rep._submodules)}")

t0 = time.monotonic()
rho = rep.permutation_representation()
stage(f"B1 permutation_representation {time.monotonic() - t0:.1f}s (codomain {rho.codomain()})")

t0 = time.monotonic()
fstab = rep.finite_stabilizer()
stage(f"B2 finite_stabilizer {time.monotonic() - t0:.1f}s")

# B3 micro-measurement: the exact comprehension body of _orbit_schreier_morphism, 200 iterations
gset = rep.g_set()
orbit, words, augmentation = gset._orbit_of(rep._selected_point())
ranking = orbit.ranking_map()
from sage_extensions.categories.sets.set_categories import Sets
word_section = Sets().Mor(orbit, augmentation.domain())(lambda image: words[int(ranking(image))])
free = augmentation.domain()
letters = tuple(free.group_generators())
points = tuple(orbit)
stage(f"B3 dims: points={len(points)} letters={len(letters)} -> {len(points) * len(letters)} words")
aug_letters = tuple(augmentation(letter) for letter in letters)
stage("B3 letter images cached")
acts = {}
t0 = time.monotonic()
rows = 0
TEST_N = 200
for point in points:
    for s, letter_img in zip(letters, aug_letters):
        val = (~word_section(orbit.act(letter_img, point))) * s * word_section(point)
        rows += 1
        if rows % 50 == 0:
            elapsed = time.monotonic() - t0
            stage(f"B3 rows {rows}: {elapsed:.2f}s so far ({1000 * elapsed / rows:.1f} ms/row)")
        if rows >= TEST_N:
            break
    if rows >= TEST_N:
        break
elapsed = time.monotonic() - t0
total_rows = len(points) * len(letters)
stage(f"B3 estimate: {1000 * elapsed / rows:.1f} ms/row -> full words tuple ~{total_rows * elapsed / rows:.0f}s")

# per-operation split within a row, on a fresh sample
import collections
splits = collections.defaultdict(float)
s = letters[0]
g = aug_letters[0]
point = points[1]
for _ in range(30):
    t1 = time.monotonic(); a1 = orbit.act(g, point); splits["act"] += time.monotonic() - t1
    t1 = time.monotonic(); w1 = word_section(a1); splits["section(act_p)"] += time.monotonic() - t1
    t1 = time.monotonic(); w2 = word_section(point); splits["section(p)"] += time.monotonic() - t1
    t1 = time.monotonic(); v = (~w1) * s * w2; splits["word ops"] += time.monotonic() - t1
for name, total in splits.items():
    stage(f"B3 row split {name}: {1000 * total / 30:.1f} ms avg")

# B3 full (pay it once)
t0 = time.monotonic()
sch = gset.orbit_schreier_morphism(rep._selected_point())
stage(f"B3 orbit_schreier_morphism FULL {time.monotonic() - t0:.1f}s")

t0 = time.monotonic()
stabsub = rho.preimage_subgroup(fstab, generating_morphism=sch)
stage(f"B4 preimage_subgroup construct {time.monotonic() - t0:.1f}s")

t0 = time.monotonic()
gens = tuple(stabsub.group_generators())
stage(f"B4 group_generators() {time.monotonic() - t0:.1f}s: {len(gens)} generators")

stage("ATOMB_OK")
