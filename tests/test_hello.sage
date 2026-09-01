"""SageMath test suite for sage_indefinite_port."""

from sage_indefinite_port import hello_world, version


def test_hello_world() -> None:
    """Verify hello_world returns greeting string."""
    assert hello_world() == "Hello from sage-indefinite-port!"


def test_version() -> None:
    """Verify package version is set."""
    assert version() == "0.1.0"
