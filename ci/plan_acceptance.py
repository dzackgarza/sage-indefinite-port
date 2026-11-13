"""Plan the Acceptance jobs from the derived per-test budgets.

The nodes of a corpus are packed into jobs by their budgets in
ci/test_budgets.json, largest first, each node into the least-loaded job that
still has room (longest-processing-time scheduling, Graham 1969). A job holds at
most JOB_BUDGET seconds of budget, so the number of jobs follows the work in
the corpus. A node whose own budget exceeds JOB_BUDGET runs alone (#79).

    python3 ci/plan_acceptance.py matrix sample|full   # workflow matrix
    python3 ci/plan_acceptance.py nodes sample|full JOB  # the node ids of one job
"""

import heapq
import json
import math
import sys
from pathlib import Path

CI = Path(__file__).parent
JOB_BUDGET = 5 * 3600  # GitHub ends a job at 6 hours; the rest is setup.
SETUP_MINUTES = 20  # checkout, Sage and Julia provisioning, provider clone
JOB_LIMIT_MINUTES = 360


def corpus_nodes(corpus):
    budgets = json.loads((CI / "test_budgets.json").read_text())["budgets"]
    if corpus == "sample":
        nodes = (CI / "acceptance_sample.txt").read_text().splitlines()
    elif corpus == "full":
        nodes = sorted(node for node in budgets if node.startswith("tests/acceptance/"))
    else:
        raise SystemExit(f"unknown corpus {corpus!r}; expected sample or full")
    missing = [node for node in nodes if node not in budgets]
    assert not missing, f"nodes without a budget: {missing}"
    return {node: budgets[node] for node in nodes}


def plan(corpus):
    budgets = corpus_nodes(corpus)
    jobs = []  # (load, index, nodes)
    heap = []
    for node in sorted(budgets, key=lambda n: (-budgets[n], n)):
        cost = budgets[node]
        if heap and heap[0][0] + cost <= JOB_BUDGET:
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
