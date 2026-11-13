"""Per-test time budgets for pytest-timeout under Sage.

Budgets live in ``ci/test_budgets.json`` (one source of truth, keyed by node id;
its ``rule`` field states how each value was set). Collection gives every test
``pytest.mark.timeout(budget)``; pytest-timeout enforces it with SIGALRM.

Sage patch: when SIGALRM arrives inside a cysignals ``sig_on`` block (GAP,
PARI, FLINT), cysignals raises ``AlarmInterrupt``, a ``KeyboardInterrupt``
subclass, and pytest ends the whole session on it. The call wrapper turns that
interrupt into pytest-timeout's ordinary per-test failure, keeping the stack, so
the next test still runs.
"""

import json
from pathlib import Path

import pytest
from cysignals.signals import AlarmInterrupt

_TABLE = json.loads((Path(__file__).parent / "test_budgets.json").read_text())


def pytest_collection_modifyitems(items):
    for item in items:
        budget = _TABLE["budgets"].get(item.nodeid, _TABLE["default_seconds"])
        item.add_marker(pytest.mark.timeout(budget))
        item.user_properties.append(("budget_seconds", budget))


@pytest.hookimpl(wrapper=True)
def pytest_runtest_call(item):
    try:
        return (yield)
    except AlarmInterrupt as interrupt:
        budget = item.get_closest_marker("timeout").args[0]
        raise pytest.fail.Exception(f"Timeout (>{budget}s) inside a cysignals block", pytrace=True) from interrupt
