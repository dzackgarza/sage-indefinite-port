"""Write ci/test_budgets.json from the cost models in ci/test_budget_models.

Each model module states, per collected test node, the time a correct
implementation needs on a GitHub ubuntu-latest runner, derived from the
complexity theory in references/complexity/. Budgets are not measured runtimes.
"""

import importlib.util
import json
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
            budgets[node], basis[node], owner[node] = round(seconds, 3), f"[{group}] {reason}", group
    nodes = (MODELS / "nodes.txt").read_text().splitlines()
    missing = [node for node in nodes if node not in budgets]
    assert not missing, missing[:20]
    out = {
        "source": "ci/build_test_budgets.py over ci/test_budget_models, from references/complexity",
        "budgets": {node: budgets[node] for node in nodes},
        "basis": {node: basis[node] for node in nodes},
    }
    (Path(__file__).parent / "test_budgets.json").write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
    print(len(nodes), "budgets")


if __name__ == "__main__":
    main()
