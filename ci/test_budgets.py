"""Per-test time budgets for pytest-timeout under Sage, with a durable per-test record.

``ci/test_budgets.json`` maps every test node id to its estimated maximum CI
seconds. A parametrized case has its own entry: two lattices of one family can
differ in cost by orders of magnitude. It is versioned with the code, so a raised budget is visible in git
history beside the change that slowed the test. A collected test without an
entry stops the session at collection.

Each finished test prints one ``BUDGET`` line (node, outcome, duration, budget)
to the log at once; ``--report-log`` keeps the same record as JSON lines.

Sage patch: when pytest-timeout's SIGALRM arrives inside a cysignals ``sig_on``
block (GAP, PARI, FLINT), cysignals raises ``AlarmInterrupt``, a
``KeyboardInterrupt`` subclass, and pytest ends the whole session on it. The
call wrapper turns that interrupt into pytest-timeout's per-test failure,
keeping the stack, so the next test still runs.
"""

import json
import os
from pathlib import Path

import pytest
from cysignals.signals import AlarmInterrupt

_RECORD = json.loads((Path(__file__).parent / "test_budgets.json").read_text())
# A model budget is the time a correct implementation needs, a target, not a
# kill limit. A node with a recorded price (least measured passing time) is
# killed at PRICE_FACTOR times that price plus WARMUP_SECONDS. PRICE_FACTOR covers
# the gap between the local measuring machine and the CI runner; WARMUP_SECONDS
# covers runner jitter and the first GAP or PARI call in a process, which
# dominate a sub-second node.
PRICE_FACTOR = 3
WARMUP_SECONDS = 10
# The probe corpus is a pricing run, not acceptance: every node is killed at
# one flat limit, and a killed node is priced as at least that many seconds.
PROBE_KILL_SECONDS = 60
# Completion pricing for nodes unfinished in Acceptance run 44000728665.
# This is a retry limit, not a measured price or a correctness budget.
PROBE_LONG_KILL_SECONDS = 600
_BUDGETS = _RECORD["budgets"] | {
    node: PRICE_FACTOR * price + WARMUP_SECONDS for node, price in _RECORD["prices"].items()
}
if os.environ["ACCEPTANCE_CORPUS"] == "probe":
    _BUDGETS = {node: PROBE_KILL_SECONDS for node in _BUDGETS}
elif os.environ["ACCEPTANCE_CORPUS"] == "probe_long":
    _BUDGETS = {node: PROBE_LONG_KILL_SECONDS for node in _BUDGETS}


def pytest_collection_modifyitems(items):
    missing = [item.nodeid for item in items if item.nodeid not in _BUDGETS]
    if missing:
        raise pytest.UsageError("tests without a budget in ci/test_budgets.json:\n" + "\n".join(missing))
    for item in items:
        budget = _BUDGETS[item.nodeid]
        item.add_marker(pytest.mark.timeout(budget))
        item.user_properties.append(("budget_seconds", budget))


@pytest.hookimpl(wrapper=True)
def pytest_runtest_call(item):
    try:
        return (yield)
    except AlarmInterrupt as interrupt:
        budget = _BUDGETS[item.nodeid]
        raise pytest.fail.Exception(f"Timeout (>{budget}s) inside a cysignals block", pytrace=True) from interrupt


def pytest_runtest_logreport(report):
    if report.when == "call" or report.outcome != "passed":
        print(f"\nBUDGET\t{report.nodeid}\t{report.when}\t{report.outcome}\t{report.duration:.1f}s\t{_BUDGETS[report.nodeid]}s", flush=True)
