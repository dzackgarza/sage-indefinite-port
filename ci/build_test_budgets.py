"""Write ci/test_budgets.json from the cost models in ci/test_budget_models.

Each model module states, per collected test node, the time a correct
implementation needs on a GitHub ubuntu-latest runner, derived from the
complexity theory in references/complexity/. Budgets are not measured runtimes.
"""

import importlib.util
import json
import re
from pathlib import Path

MODELS = Path(__file__).parent / "test_budget_models"
# E owns the cubic-division Voronoi witnesses that F also priced; G re-derives part of E's scope.
ORDER = "ABCDFEG"


def _load(group):
    spec = importlib.util.spec_from_file_location(f"budget_model_{group}", MODELS / f"{group}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.estimates()


def main():
    budgets, basis, owner = {}, {}, {}
    for group in ORDER:
        for node, (seconds, reason) in _load(group).items():
            if node in budgets:
                assert (group, owner[node]) in {("E", "F"), ("G", "E")}, (node, group, owner[node])
            assert seconds > 0, (node, seconds)
            budgets[node], basis[node], owner[node] = float(f"{seconds:.3g}"), f"[{group}] {reason}", group
    # The model ledger contains historical test nodes. Only tests that still
    # exist in the independently sourced acceptance corpus may be scheduled.
    current_functions = {
        f"tests/acceptance/{path.name}": set(
            re.findall(r"(?m)^def (test_\w+)\(", path.read_text())
        )
        for path in (MODELS.parent.parent / "tests/acceptance").glob("test_*.sage")
    }
    nodes = [
        node for node in (MODELS / "nodes.txt").read_text().splitlines()
        if "::" in node
        and node.split("::", 1)[0] in current_functions
        and node.split("::", 1)[1].split("[", 1)[0] in current_functions[node.split("::", 1)[0]]
    ]
    missing = [node for node in nodes if node not in budgets]
    assert not missing, missing[:20]
    path = Path(__file__).parent / "test_budgets.json"
    out = {
        # Recorded prices are measurements, not derived budgets: they are kept.
        "prices": {
            node: price
            for node, price in json.loads(path.read_text())["prices"].items()
            if node in nodes
        },
        "source": "ci/build_test_budgets.py over ci/test_budget_models, from references/complexity",
        "budgets": {node: budgets[node] for node in nodes},
        "basis": {node: basis[node] for node in nodes},
    }
    path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n")
    print(len(nodes), "budgets")


if __name__ == "__main__":
    main()
