"""Plan the Acceptance jobs from the derived per-test budgets.

The nodes of a corpus are packed into jobs by their costs, largest first, each node into the least-loaded job that
still has room (longest-processing-time scheduling, Graham 1969). A job holds at
most JOB_BUDGET seconds of budget, so the number of jobs follows the work in
the corpus. A node whose own budget exceeds JOB_BUDGET runs alone (#79).

The cost of a tier1 or tier2 node is its recorded price ("prices" in
ci/test_budgets.json: the least recorded passed call time, the "Acceptance
tiers" table of references/complexity/triage.md), so the tier job timeouts
come from the recorded prices. Every tier1 node has a price. An unpriced tier2
node and every full node cost their derived budget.

    python3 ci/plan_acceptance.py matrix tier1|tier2|full|probe|probe_long
    python3 ci/plan_acceptance.py nodes tier1|tier2|full|probe|probe_long JOB

tier1 is ci/acceptance_tier1.txt (the smallest independent node per pinned
branch), tier2 is ci/acceptance_tier2.txt (the published-result reproductions
of #64), full is every acceptance node (references/complexity/triage.md), and
probe is ci/acceptance_probe.txt, the seeded pricing sample of every family
(up to five nodes per test function, the ListReflect and IsotropicCases
functions split by rank). A probe node costs PROBE_KILL_SECONDS of
ci/test_budgets.py, its flat kill limit, and a probe job holds at most
PROBE_JOB_BUDGET seconds.

probe_long is ci/acceptance_probe_long.txt, the unfinished nodes of probe run
44000728665. Each receives at least PROBE_LONG_KILL_SECONDS for completion
pricing, retaining the existing allowance when its recorded price requires
more time. Another timeout remains a lower bound, never a completed measurement.
"""

import ast
import heapq
import json
import math
import sys
from pathlib import Path

CI = Path(__file__).parent
JOB_BUDGET = 5 * 3600  # GitHub ends a job at 6 hours; the rest is setup.
SETUP_MINUTES = 20  # checkout, Sage and Julia provisioning, provider clone
JOB_LIMIT_MINUTES = 360
PROBE_JOB_BUDGET = 30 * 60


def probe_kill_seconds(corpus, node, prices):
    # ci/test_budgets.py imports cysignals, absent on the plan runner; read its constant.
    tree = ast.parse((CI / "test_budgets.py").read_text())
    constant = {"probe": "PROBE_KILL_SECONDS", "probe_long": "PROBE_LONG_KILL_SECONDS"}[corpus]
    constants = {
        assignment.targets[0].id: assignment.value.value
        for assignment in tree.body
        if isinstance(assignment, ast.Assign) and isinstance(assignment.value, ast.Constant)
    }
    limit = constants[constant]
    if corpus == "probe_long" and node in prices:
        return max(limit, constants["PRICE_FACTOR"] * prices[node] + constants["WARMUP_SECONDS"])
    return limit


def corpus_nodes(corpus):
    record = json.loads((CI / "test_budgets.json").read_text())
    budgets = record["budgets"]
    if corpus in ("probe", "probe_long"):
        nodes = (CI / f"acceptance_{corpus}.txt").read_text().splitlines()
        missing = [node for node in nodes if node not in budgets]
        assert not missing, f"nodes without a budget: {missing}"
        return {node: probe_kill_seconds(corpus, node, record["prices"]) for node in nodes}
    if corpus in ("tier1", "tier2"):
        nodes = (CI / f"acceptance_{corpus}.txt").read_text().splitlines()
    elif corpus == "full":
        nodes = sorted(node for node in budgets if node.startswith("tests/acceptance/"))
    else:
        raise SystemExit(f"unknown corpus {corpus!r}; expected tier1, tier2, full, probe or probe_long")
    missing = [node for node in nodes if node not in budgets]
    assert not missing, f"nodes without a budget: {missing}"
    if corpus == "full":
        return {node: budgets[node] for node in nodes}
    prices = json.loads((CI / "test_budgets.json").read_text())["prices"]
    unpriced = [node for node in nodes if node not in prices]
    assert corpus == "tier2" or not unpriced, f"tier1 nodes without a recorded price: {unpriced}"
    return {node: prices.get(node, budgets[node]) for node in nodes}


def plan(corpus):
    budgets = corpus_nodes(corpus)
    job_budget = PROBE_JOB_BUDGET if corpus in ("probe", "probe_long") else JOB_BUDGET
    jobs = []  # (load, index, nodes)
    heap = []
    for node in sorted(budgets, key=lambda n: (-budgets[n], n)):
        cost = budgets[node]
        if heap and heap[0][0] + cost <= job_budget:
            load, index = heapq.heappop(heap)
            jobs[index][0] += cost
            jobs[index][1].append(node)
            heapq.heappush(heap, (jobs[index][0], index))
            continue
        jobs.append([cost, [node]])
        heapq.heappush(heap, (cost, len(jobs) - 1))
    return jobs


def main():
    command, corpus = sys.argv[1], sys.argv[2]
    jobs = plan(corpus)
    if command == "matrix":
        include = [
            {
                "family": "acceptance",
                "shard": index,
                "shards": len(jobs),
                "selection": corpus,
                "timeout": min(JOB_LIMIT_MINUTES, math.ceil(load / 60) + SETUP_MINUTES),
            }
            for index, (load, _nodes) in enumerate(jobs)
        ]
        print("matrix=" + json.dumps({"include": include}))
    elif command == "nodes":
        print("\n".join(jobs[int(sys.argv[3])][1]))
    else:
        raise SystemExit(f"unknown command {command!r}; expected matrix or nodes")


if __name__ == "__main__":
    main()
