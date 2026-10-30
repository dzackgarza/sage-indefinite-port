"""Regression tests for acceptance artifact and node completeness checks."""

import contextlib
import importlib.util
import io
import json
import subprocess
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


def test_batch_runner_continues_after_malformed_batch_junit(tmp_path, monkeypatch, capsys) -> None:
    batches = _load_ci_module("pytest_batches")
    first = "tests/acceptance/test_reflective.sage::test_first"
    second = "tests/acceptance/test_reflective.sage::test_second"
    calls = 0

    def fake_run(command, **kwargs):
        nonlocal calls
        calls += 1
        if "--collect-only" in command:
            return subprocess.CompletedProcess(command, 0, stdout=f"{first}\n{second}\n", stderr="")
        report = Path(next(arg.split("=", 1)[1] for arg in command if arg.startswith("--junitxml=")))
        if first in command:
            report.write_text("<testsuites><testcase>")
        else:
            suite = ET.Element("testsuites")
            child = ET.SubElement(suite, "testsuite")
            ET.SubElement(child, "testcase", classname="tests.acceptance.test_reflective.sage", name="test_second")
            ET.ElementTree(suite).write(report)
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr(batches.subprocess, "run", fake_run)
    monkeypatch.setattr(
        "sys.argv",
        ["pytest_batches.py", "--sage", "sage", "--out", str(tmp_path), "--batch-size", "1", "--", "tests/acceptance/test_reflective.sage"],
    )

    result = batches.main()

    assert result == 1
    assert calls == 3
    assert "produced malformed JUnit XML" in capsys.readouterr().err
    aggregate = ET.parse(tmp_path / "junit.xml").getroot()
    assert [case.get("name") for case in aggregate.iter("testcase")] == ["test_second"]


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


def test_acceptance_summary_rejects_mixed_port_or_provider_revisions(tmp_path) -> None:
    summary = _load_ci_module("acceptance_summary")
    matrix = json.dumps(
        {
            "include": [
                {"family": "perfect_domains", "shard": 0, "shards": 2},
                {"family": "perfect_domains", "shard": 1, "shards": 2},
            ]
        }
    )
    for shard_id, port, provider in ((0, "port-a", "research-a"), (1, "port-b", "research-b")):
        shard = tmp_path / f"perfect_domains-{shard_id}"
        shard.mkdir()
        suite = ET.Element("testsuites")
        ET.SubElement(suite, "testsuite")
        ET.ElementTree(suite).write(shard / "junit.xml")
        (shard / "commit").write_text(port)
        (shard / "research_commit").write_text(provider)
    output = io.StringIO()

    with contextlib.redirect_stdout(output):
        result = summary.main(tmp_path, matrix)

    assert result == 1
    assert "port revisions are not uniquely identified" in output.getvalue()
    assert "research revisions are not uniquely identified" in output.getvalue()


def test_acceptance_summary_rejects_missing_revision_metadata(tmp_path) -> None:
    summary = _load_ci_module("acceptance_summary")
    matrix = json.dumps({"include": [{"family": "perfect_domains", "shard": 0, "shards": 1}]})
    shard = tmp_path / "perfect_domains-0"
    shard.mkdir()
    suite = ET.Element("testsuites")
    ET.SubElement(suite, "testsuite")
    ET.ElementTree(suite).write(shard / "junit.xml")
    (shard / "commit").write_text("port-revision")
    output = io.StringIO()

    with contextlib.redirect_stdout(output):
        result = summary.main(tmp_path, matrix)

    assert result == 1
    assert "perfect_domains-0 (research_commit missing or empty)" in output.getvalue()

    (shard / "research_commit").write_text("   \n")
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        result = summary.main(tmp_path, matrix)
    assert result == 1
    assert "perfect_domains-0 (research_commit missing or empty)" in output.getvalue()


def test_acceptance_summary_rejects_a_failed_family_matrix(tmp_path) -> None:
    summary = _load_ci_module("acceptance_summary")
    matrix = json.dumps({"include": [{"family": "perfect_domains", "shard": 0, "shards": 1}]})
    shard = tmp_path / "perfect_domains-0"
    shard.mkdir()
    suites = ET.Element("testsuites")
    suite = ET.SubElement(suites, "testsuite")
    ET.SubElement(suite, "testcase", classname="sample", name="passed")
    ET.ElementTree(suites).write(shard / "junit.xml")
    (shard / "commit").write_text("port-revision")
    (shard / "research_commit").write_text("provider-revision")
    output = io.StringIO()

    with contextlib.redirect_stdout(output):
        result = summary.main(tmp_path, matrix, family_result="failure")

    assert result == 1
    assert "family matrix job result is failure, not success" in output.getvalue()


def test_acceptance_summary_rejects_an_empty_junit_report(tmp_path) -> None:
    summary = _load_ci_module("acceptance_summary")
    matrix = json.dumps({"include": [{"family": "perfect_domains", "shard": 0, "shards": 1}]})
    shard = tmp_path / "perfect_domains-0"
    shard.mkdir()
    ET.ElementTree(ET.Element("testsuites")).write(shard / "junit.xml")
    (shard / "commit").write_text("port-revision")
    (shard / "research_commit").write_text("provider-revision")
    output = io.StringIO()

    with contextlib.redirect_stdout(output):
        result = summary.main(tmp_path, matrix)

    assert result == 1
    assert "perfect_domains-0 (JUnit report has no testcases)" in output.getvalue()


def test_acceptance_summary_reports_malformed_junit_without_crashing(tmp_path) -> None:
    summary = _load_ci_module("acceptance_summary")
    matrix = json.dumps({"include": [{"family": "perfect_domains", "shard": 0, "shards": 1}]})
    shard = tmp_path / "perfect_domains-0"
    shard.mkdir()
    (shard / "junit.xml").write_text("<testsuites><testcase>")
    (shard / "commit").write_text("port-revision")
    (shard / "research_commit").write_text("provider-revision")
    output = io.StringIO()

    with contextlib.redirect_stdout(output):
        result = summary.main(tmp_path, matrix)

    assert result == 1
    assert "perfect_domains-0 (malformed JUnit report:" in output.getvalue()


def test_acceptance_summary_reports_a_failed_plan_without_valid_matrix_json(tmp_path) -> None:
    summary = _load_ci_module("acceptance_summary")
    output = io.StringIO()

    with contextlib.redirect_stdout(output):
        result = summary.main(tmp_path / "no-results", matrix_json="", family_result="skipped", plan_result="failure")

    assert result == 1
    assert "acceptance matrix planning job result is failure, not success" in output.getvalue()
    assert "acceptance matrix JSON is missing or invalid" in output.getvalue()
    assert "family matrix job result is skipped, not success" in output.getvalue()


def test_acceptance_summary_rejects_empty_or_malformed_matrix_includes(tmp_path) -> None:
    summary = _load_ci_module("acceptance_summary")
    for matrix, message in (
        ({"include": []}, "acceptance matrix include list is empty"),
        ({"include": [{"family": "perfect_domains", "shard": True}]}, "acceptance matrix include entries are malformed"),
        ({"include": ["perfect_domains-0"]}, "acceptance matrix include entries are malformed"),
        ({"include": [{"family": "perfect_domains", "shard": 0, "shards": 2}]}, "acceptance matrix shard indices are incomplete or duplicated for perfect_domains"),
        ({"include": [{"family": "perfect_domains", "shard": 0, "shards": 1}, {"family": "perfect_domains", "shard": 0, "shards": 1}]}, "acceptance matrix shard indices are incomplete or duplicated for perfect_domains"),
    ):
        output = io.StringIO()
        with contextlib.redirect_stdout(output):
            result = summary.main(tmp_path / "no-results", json.dumps(matrix))
        assert result == 1
        assert message in output.getvalue()
