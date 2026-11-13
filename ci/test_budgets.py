"""Per-test time budgets for pytest-timeout under Sage, with a durable per-test record.

``ci/test_budgets.json`` maps every test node id to its estimated maximum CI
seconds. It is versioned with the code, so a raised budget is visible in git
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
from pathlib import Path

import pytest
from cysignals.signals import AlarmInterrupt

_BUDGETS = json.loads((Path(__file__).parent / "test_budgets.json").read_text())["budgets"]


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
