"""Run exactly one explicitly selected mathematical regression per edit.

This plugin rejects an entire module or a parametrized function expanded to
multiple items before running tests.  The independent mathematical argument
belongs with the focused regression, not in this execution adapter.
"""

import pytest


class ExactlyOne:
    def pytest_collection_modifyitems(self, items):
        if len(items) != 1:
            raise pytest.UsageError(
                "Select exactly one mathematical case; qualify parametrized nodes with [case]"
            )


def run(node: str) -> int:
    return pytest.main(
        ["-q", "--maxfail=1", "--disable-warnings", node],
        plugins=(ExactlyOne(),),
    )
