"""SageMath indefinite lattice and orthogonal group port."""

__version__ = "0.1.0"


def version() -> str:
    """Return package version."""
    return __version__


def hello_world() -> str:
    """Return greeting string for sage-indefinite-port package."""
    return "Hello from sage-indefinite-port!"


__all__ = ["__version__", "hello_world", "version"]
