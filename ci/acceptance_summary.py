"""Collate the acceptance shards' JUnit reports into one failure list.

Usage: acceptance_summary.py RESULTS_DIR

RESULTS_DIR holds one directory per shard (<family>-<shard>, from `just
ci-acceptance`), each with junit.xml, commit and research_commit. Writes a
Markdown report to stdout and exits 1 when any shard lacks a report or any case
did not pass.
"""

import json
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


def main(results_dir: Path, matrix_json: str | None = None) -> int:
    shard_dirs = sorted(path for path in results_dir.iterdir() if path.is_dir()) if results_dir.exists() else []
    if not shard_dirs:
        if matrix_json is None:
            raise SystemExit(f"no shard results in {results_dir}")
    counts: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
    failures: dict[str, list[str]] = defaultdict(list)
    missing: list[str] = []
    if matrix_json is not None:
        matrix = json.loads(matrix_json)
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
        report = shard / "junit.xml"
        if not report.exists():
            missing.append(shard.name)
            continue
        for case in ET.parse(report).getroot().iter("testcase"):
            node = f"{case.get('classname')}::{case.get('name')}"
            outcome = next((child.tag for child in case if child.tag in ("failure", "error", "skipped")), "passed")
            counts[family][outcome] += 1
            if outcome != "passed":
                detail = next(child for child in case if child.tag == outcome).get("message", "")
                failures[family].append(f"- `{node}` ({outcome}): {detail.splitlines()[0] if detail else ''}")
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
    if len(sys.argv) not in (2, 3):
        raise SystemExit("usage: acceptance_summary.py RESULTS_DIR [MATRIX_JSON]")
    sys.exit(main(Path(sys.argv[1]), sys.argv[2] if len(sys.argv) == 3 else None))
