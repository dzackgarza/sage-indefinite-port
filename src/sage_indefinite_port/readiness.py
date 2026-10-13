"""Entry points whose work unit is still open refuse every call immediately.

A port entry point the research preamble calls is marked ``@unfinished(N)`` while its work
unit, issue N, is open. The marked function raises on its first line, so a consumer and an
acceptance case reach the missing capability in no time instead of after the expensive
prerequisites that precede it. The work unit's last commit removes the marker; from then on
its acceptance cases run in full.
"""

from collections.abc import Callable
from functools import wraps


class UnfinishedCapability(NotImplementedError):
    """A port entry point was called while its owning work unit is open."""


UNFINISHED: dict[str, int] = {}
"""``module:qualname`` of every marked entry point, to the issue that owns it."""


def unfinished[**P, R](issue: int) -> Callable[[Callable[P, R]], Callable[P, R]]:
    def mark(function: Callable[P, R]) -> Callable[P, R]:
        @wraps(function)
        def refuse(*args: P.args, **kwargs: P.kwargs) -> R:
            raise UnfinishedCapability(f"{function.__module__}.{function.__qualname__} is unfinished: work unit #{issue} owns it")

        UNFINISHED[f"{function.__module__}:{function.__qualname__}"] = issue
        return refuse

    return mark
