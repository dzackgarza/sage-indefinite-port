"""Collate the acceptance shards' JUnit reports into one failure list.

Usage: acceptance_summary.py RESULTS_DIR [MATRIX_JSON [FAMILY_RESULT [PLAN_RESULT]]]

RESULTS_DIR holds one directory per shard (<family>-<shard>, from `just
ci-acceptance`), each with junit.xml, commit and research_commit. If MATRIX_JSON
is supplied, it is the workflow plan's `include` object and every planned shard
must have an artifact. The family and plan results must both be ``success``.
Writes a Markdown report to stdout and exits 1 when planning/execution failed,
any shard lacks an artifact/report, revisions are inconsistent, or a case failed.
"""

import json
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


def main(
    results_dir: Path,
    matrix_json: str | None = None,
    family_result: str = "success",
    plan_result: str = "success",
) -> int:
    shard_dirs = sorted(path for path in results_dir.iterdir() if path.is_dir()) if results_dir.exists() else []
    if not shard_dirs:
        if matrix_json is None:
            raise SystemExit(f"no shard results in {results_dir}")
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    failures: dict[str, list[str]] = defaultdict(list)
    missing: list[str] = []
    if plan_result != "success":
        missing.append(f"acceptance matrix planning job result is {plan_result}, not success")
    if matrix_json is not None:
        try:
            matrix = json.loads(matrix_json)
        except json.JSONDecodeError:
            matrix = None
            missing.append("acceptance matrix JSON is missing or invalid")
        if matrix is not None and (not isinstance(matrix, dict) or not isinstance(matrix.get("include"), list)):
            matrix = None
            missing.append("acceptance matrix JSON has no include list")
        if family_result != "success":
            missing.append(f"family matrix job result is {family_result}, not success")
        if matrix is not None:
            expected = {f"{item['family']}-{item['shard']}" for item in matrix["include"]}
            present = {shard.name for shard in shard_dirs}
            missing.extend(f"{name} (shard artifact missing)" for name in sorted(expected - present))
            unexpected = sorted(present - expected)
            if unexpected:
                missing.extend(f"{name} (unexpected shard artifact)" for name in unexpected)
    commits = {(shard / "commit").read_text().strip() for shard in shard_dirs if (shard / "commit").exists()}
    research = {(shard / "research_commit").read_text().strip() for shard in shard_dirs if (shard / "research_commit").exists()}
    for shard in shard_dirs:
        family = shard.name.rsplit("-", 1)[0]
        for revision_file in ("commit", "research_commit"):
            revision_path = shard / revision_file
            if not revision_path.exists() or not revision_path.read_text().strip():
                missing.append(f"{shard.name} ({revision_file} missing or empty)")
        report = shard / "junit.xml"
        if not report.exists():
            missing.append(shard.name)
            continue
        cases = tuple(ET.parse(report).getroot().iter("testcase"))
        if not cases:
            missing.append(f"{shard.name} (JUnit report has no testcases)")
        for case in cases:
            node = f"{case.get('classname')}::{case.get('name')}"
            outcome = next((child.tag for child in case if child.tag in ("failure", "error", "skipped")), "passed")
            counts[family][outcome] += 1
            if outcome != "passed":
                detail = next(child for child in case if child.tag == outcome).get("message", "")
                failures[family].append(f"- `{node}` ({outcome}): {detail.splitlines()[0] if detail else ''}")
    if shard_dirs and len(commits) != 1:
        missing.append(f"port revisions are not uniquely identified: {sorted(commits)!r}")
    if shard_dirs and len(research) != 1:
        missing.append(f"research revisions are not uniquely identified: {sorted(research)!r}")
    print(f"# Acceptance at port {', '.join(sorted(commits))}, research {', '.join(sorted(research))}\n")
    print("| Family | passed | failure | error | skipped |\n| --- | --- | --- | --- | --- |")
    for family in sorted(counts):
        row = counts[family]
        print(f"| {family} | {row['passed']} | {row['failure']} | {row['error']} | {row['skipped']} |")
    if missing:
        print("\n## Missing, unexpected, or report-less shards\n")
        print("\n".join(f"- {name}" for name in missing))
    for family in sorted(failures):
        print(f"\n## {family}\n")
        print("\n".join(failures[family]))
    return 1 if missing or failures else 0


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3, 4, 5):
        raise SystemExit("usage: acceptance_summary.py RESULTS_DIR [MATRIX_JSON [FAMILY_RESULT [PLAN_RESULT]]]")
    sys.exit(
        main(
            Path(sys.argv[1]),
            sys.argv[2] if len(sys.argv) >= 3 else None,
            sys.argv[3] if len(sys.argv) >= 4 else "success",
            sys.argv[4] if len(sys.argv) >= 5 else "success",
        )
    )
