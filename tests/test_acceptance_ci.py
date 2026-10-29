"""Regression tests for acceptance artifact and node completeness checks."""

import contextlib
import importlib.util
import io
import json
import xml.etree.ElementTree as ET
from pathlib import Path


def _load_ci_module(name: str):
    path = Path(__file__).resolve().parents[1] / "ci" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_acceptance_summary_fails_when_a_matrix_shard_has_no_artifact(tmp_path) -> None:
    summary = _load_ci_module("acceptance_summary")
    matrix = json.dumps({"include": [{"family": "reflective", "shard": 0, "shards": 1}]})
    output = io.StringIO()

    with contextlib.redirect_stdout(output):
        result = summary.main(tmp_path / "missing-results", matrix)

    assert result == 1
    assert "reflective-0 (shard artifact missing)" in output.getvalue()


def test_junit_node_identity_is_not_replaced_by_an_equal_result_count() -> None:
    batches = _load_ci_module("pytest_batches")
    node = (
        "tests/acceptance/test_perfect_domains.sage::"
        "test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_19]"
    )
    expected = (
        "tests.acceptance.test_perfect_domains.sage",
        "test_perfect_domain_orbit_counts_match_polyhedral_common[Result_Enumeration_19]",
    )
    assert batches.junit_key(node) == expected
    assert batches.junit_nodes_match([node], [{"classname": expected[0], "name": expected[1]}])
    assert not batches.junit_nodes_match([node], [{"classname": expected[0], "name": "different-test"}])
    assert not batches.junit_nodes_match([node], [{"classname": expected[0], "name": expected[1]}] * 2)


def test_acceptance_summary_accepts_a_present_shard_with_a_junit_report(tmp_path) -> None:
    summary = _load_ci_module("acceptance_summary")
    matrix = json.dumps({"include": [{"family": "perfect_domains", "shard": 0, "shards": 1}]})
    shard = tmp_path / "perfect_domains-0"
    shard.mkdir()
    suites = ET.Element("testsuites")
    suite = ET.SubElement(suites, "testsuite")
    ET.SubElement(suite, "testcase", classname="tests.acceptance.test_perfect_domains.sage", name="test_case")
    ET.ElementTree(suites).write(shard / "junit.xml")
    (shard / "commit").write_text("port-revision")
    (shard / "research_commit").write_text("provider-revision")
    output = io.StringIO()

    with contextlib.redirect_stdout(output):
        result = summary.main(tmp_path, matrix)

    assert result == 0, output.getvalue()
    assert "| perfect_domains | 1 | 0 | 0 | 0 |" in output.getvalue()
