# Enriques reproduction setup for the `sage-extensions` owner — measured, not reconstructed

Everything below was checked against what is on this VM right now. Where a value is a pin it was verified by
comparing content, not by trusting a directory name.

**Status of the heavy slot at the time of writing: free.** No port process is running, so the whole machine is
the library's for its probe. This document is a reproduction handoff, not a request to launch anything here.

## 0. Why the library's import fails, and the one-line fix

The library bound the archived checkout `~/port-a4d4aa4` at **`73ef6713`** ("docs: use the shared runtime in local
reproduction commands"). That commit is an **ancestor** of the port's migration to the shared library:

```
$ git merge-base --is-ancestor 73ef6713 61ac45a0 ; echo $?     # 0 — it is pre-migration
$ git log -1 --format='%h %ci %s' 61ac45a0
61ac45a0 2026-10-31 23:26:04 +0000 feat: migrate the port to sage_extensions
```

and its source imports the pre-migration package:

```
$ grep -n dzack_research ~/port-a4d4aa4/src/sage_indefinite_port/indefinite/recursive.py
23:from dzack_research.preamble.all import ZZ, Lattices
24:from dzack_research.preamble.categories.group.groups import OwnedGroups
```

That is exactly the `dzack_research.preamble.categories` import that fails, and it fails at module import, before
any provider stage can run. The current port has no such import at all:

```
$ grep -rn dzack_research --include=*.py --include=*.sage src/ tests/ | wc -l
0
$ grep -n "from sage_extensions.all import" src/sage_indefinite_port/indefinite/recursive.py
25:from sage_extensions.all import GObjects, Lattices, Modules, Sets, ZZ
```

**Matching working consumer: any port revision at or after `61ac45a0`.** The concrete one every measurement below
was taken on is **`1051bf1e`** (this repository's `main`). Binding `73ef6713` cannot work with the current
provider under any settings — it is a different library.

## 1. Pins, verified

| thing | value | how it was verified |
| --- | --- | --- |
| port | this repository, `1051bf1e` | `git rev-parse HEAD`; every run below used it |
| provider staging | `~/sage-extensions-42d7361c` | **byte-identical** to commit `42d7361cfda3d8505f365b26d9aa687447d6df9d` of `~/sage-image/rootfs/repos/sage-extensions`: full-tree `diff -rq` excluding caches is empty, and `pyproject.toml` compares equal |
| provider archive source | `~/sage-image/rootfs/repos/sage-extensions` | `git cat-file -t 42d7361c…` → `commit`; `git archive 42d7361c… | tar -t | wc -l` → 1096 entries |
| pytest plugins | `~/ai-review-ci/tool-artifacts/pytest_plugins` | present (bound to `/ci/tool-artifacts/pytest_plugins`) |
| fixture | `tests/fixtures/enriques_87_polarizations.json` (63 041 B) | in the port repository |
| acceptance node | `tests/acceptance/test_boundaries.sage::test_enriques_modular_groups_match_dutour_sikiric_hulek[case_1]` | collected and run |

The provider staging has **no `.git`** (it is an `git archive` extraction). Do not try to `git -C
~/sage-extensions-42d7361c rev-parse`; it returns nothing. Use the pin above, or re-extract with
`git -C ~/sage-image/rootfs/repos/sage-extensions archive 42d7361c… | tar -x -C <dest>`.

## 2. Runner and the original process-tree RSS guard

The runner is PRoot over a pinned rootfs; `run.sh` execs `env -i`, so **host environment variables do not reach
Python** — anything the process needs is set inside the `/bin/sh -c` string.

```
run.sh     ~/sage-image/run.sh        (721 B)  — mounts + `env -i` from env.oci
peakrss.py ~/sage-image/peakrss.py    (4576 B) — process-TREE RSS sampling and the guard
node_run.py ~/sage-image/node_run.py  (2382 B) — one acceptance node, cap+guard, keeps the transcript
stackdump.py ~/sage-image/stackdump.py (778 B) — 60 s faulthandler dumps (a pytest plugin)
```

The guard is the **tree** RSS, not a single process: `peakrss.py` sums the process tree and kills at
`PEAKRSS_LIMIT_MIB`, and reports `peak tree RSS` plus a per-executable breakdown. The Enriques runs used
`PEAKRSS_LIMIT_MIB=6144` with `timeout -s KILL` as the time cap (2 400 s for the staged probe, 1 800 s for the
acceptance node) — the same limits as the acceptance family, unchanged.

## 3. The proven commands

**(a) Staged probe** — this is the run that produced the stage table (4.3 s invariant, 19.5 s transport, 249 s
`coinvariant_extension_subgroup`, then > 2 127 s unfinished inside the discriminant bundle). It wrote
`~/sage-image/enriques-case1.log` on the acceptance run and its own transcript for the probe; the probe itself is
`~/sage-image/probe_enriques_stages.py` (3 999 B) and prints before every stage, so a cap still leaves the ordered
attribution.

```bash
cd $HOME/sage-image
REPO=/workspace/nucleus/work/projects/p-165a6922b36b40789ddc90db34bb23c7/repo
CI=$HOME/ai-review-ci
PROV=$HOME/sage-extensions-42d7361c
PEAKRSS_LABEL=enrstages PEAKRSS_LIMIT_MIB=6144 \
python3 peakrss.py -- timeout -s KILL 2400 ./run.sh -w /work \
  -b "$REPO:/work" \
  -b "$CI:/ci" \
  -b "$PROV:/repos/sage-extensions" \
  -b "$HOME/sage-image/probe_enriques_stages.py:/probe.py" \
  -- /bin/sh -c 'cd /work && PYTHONPATH=/ci/tool-artifacts/pytest_plugins:/work/src:/work python -u -s /probe.py'
```

Note `PYTHONPATH` is set **inside** the `sh -c` string: `/ci/tool-artifacts/pytest_plugins:/work/src:/work`. The
provider is imported from `/repos/sage-extensions` by its install, not from `PYTHONPATH`, which is why the mount
to that exact path matters.

**(b) The untouched acceptance node** — the original assertions, 1 800 s cap, 6 144 MiB guard, transcript kept:

```bash
cd $HOME/sage-image
SAGE_EXTENSIONS_TREE=$HOME/sage-extensions-42d7361c python3 node_run.py \
  "tests/acceptance/test_boundaries.sage::test_enriques_modular_groups_match_dutour_sikiric_hulek[case_1]" \
  1800 6144 "$HOME/sage-image/enriques-case1.log" -s -p stackdump \
  --bind "$HOME/sage-image/stackdump.py:/probe/stackdump.py"
```

Measured: **CAPPED**, killed at 1 800 s, peak 4 182 MiB, 30 stack dumps landing in 26 distinct innermost frames.

**(c) The short morphism probe** — `~/sage-image/probe_enriques_morphism_cost.py` (3 014 B), same invocation form
with a 900 s cap; it measured 219 s to build the coinvariant subgroup and **> 681 s unfinished** inside
`group.group_generators()`.

## 4. Artifacts, where they are and what they hold

| artifact | location | contents |
| --- | --- | --- |
| `enriques-case1.log` | `~/sage-image/enriques-case1.log` (60 615 B) | the acceptance node's full transcript: peakrss guard lines, 30 stack dumps, the CAPPED verdict |
| staged stage table | `.tmp/handoff-owner-enriques-oqn.md` §2, §2b, §2c | measured per-stage times and what was retracted |
| probe | `~/sage-image/probe_enriques_stages.py` | reproduces `_gamma_h` and stage 1 of the assertion, printing before each stage |
| probe | `~/sage-image/probe_enriques_morphism_cost.py` | isolates `group_generators()` from the discriminant morphisms |
| fixture | `tests/fixtures/enriques_87_polarizations.json` | the 87 cases; `case_1` has published order 696 729 600 and orbits 5 / 9 / 18 |
| acceptance | `tests/acceptance/test_boundaries.sage:219-230` | the unchanged node |

## 5. What the specimen and the acceptance are, unchanged

`case_1` of `test_enriques_modular_groups_match_dutour_sikiric_hulek`, with the published
`group_order = 696 729 600` and orbit counts `5 / 9 / 18` for lines, planes and flags. Both assertions stay as
they are: `int(image.cardinality()) == case["group_order"]` and `_counts(group) == (line, plane, flag)`. Nothing
in this handoff changes the assertion, the cap or the guard, and no port-side patch to the provider is proposed.
