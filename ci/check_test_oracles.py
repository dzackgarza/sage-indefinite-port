"""Require pre-existing independent oracles for new/changed mathematical tests.

This checks provenance, not mathematical correctness of the author's
derivation. That independent verification remains required by AGENTS.md.
"""

import argparse
import ast
from collections import defaultdict
from pathlib import Path
import re
import subprocess
import sys

FIELDS = ("ORACLE-FILE", "ORACLE-LOCATOR", "ORACLE-EXCERPT",
          "ORACLE-HYPOTHESES", "ORACLE-EXPECTED", "ORACLE-DERIVATION")
META = re.compile(r"(?m)^\s*(ORACLE-(?:FILE|LOCATOR|EXCERPT|HYPOTHESES|EXPECTED|DERIVATION)):\s*(\S.*?)\s*$")
HUNK = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@")
TEST = re.compile(r"(?m)^\s*def (test_\w+)\s*\(")


def git(*args, optional=False):
    result = subprocess.run(["git", *args], capture_output=True, text=True, check=False)
    if result.returncode:
        if optional:
            return None
        raise RuntimeError(f"git {' '.join(args)}: {result.stderr.strip()}")
    return result.stdout


def changed_lines(diff):
    changed = defaultdict(set)
    path = line = None
    for item in diff.splitlines():
        if item.startswith("+++ b/"):
            path, line = item[6:], None
        elif item.startswith("+++ /dev/null"):
            path = None
        elif item.startswith("@@ "):
            match = HUNK.match(item)
            line = int(match.group(1)) if match else None
        elif path is not None and line is not None:
            if item.startswith("+") and not item.startswith("+++"):
                changed[path].add(line)
                line += 1
            elif item.startswith(" "):
                line += 1
    return {p: lines for p, lines in changed.items()
            if p.startswith("tests/") and p.endswith((".py", ".sage"))}


def tests_in(path, content):
    if path.endswith(".py"):
        tree = ast.parse(content, filename=path)
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_"):
                yield node.name, node.lineno, node.end_lineno or node.lineno, ast.get_docstring(node)
        return
    # Sage files can require preparsing: inspect original function blocks.
    matches = list(TEST.finditer(content))
    lines = content.splitlines()
    for index, item in enumerate(matches):
        start = content.count("\n", 0, item.start()) + 1
        end = (content.count("\n", 0, matches[index + 1].start())
               if index + 1 < len(matches) else len(lines))
        block = "\n".join(lines[start - 1:end])
        doc = re.search(r'(?s)"""(.*?)"""', block) or re.search(r"(?s)'''(.*?)'''", block)
        yield item.group(1), start, end, doc.group(1) if doc else None


def inspect(base, diff, read_file, node=None):
    issues = []
    selected = node.split("::") if node else None
    for path, changed in changed_lines(diff).items():
        if selected and path != selected[0]:
            continue
        contents = read_file(path)
        if contents is None:
            continue
        try:
            items = list(tests_in(path, contents))
        except SyntaxError as error:
            issues.append(f"{path}: invalid test syntax: {error}")
            continue
        for name, start, end, doc in items:
            if selected and name != selected[1].split("[")[0]:
                continue
            if not any(start <= n <= end for n in changed):
                continue
            prefix = f"{path}::{name}"
            metadata = dict(META.findall(doc or ""))
            missing = [key for key in FIELDS if not metadata.get(key)]
            if missing:
                issues.append(f"{prefix}: missing {', '.join(missing)}")
                continue
            source = metadata["ORACLE-FILE"]
            if (not source.startswith(("references/vendor/", "references/oracles/"))
                    or ".." in Path(source).parts):
                issues.append(f"{prefix}: oracle is not an independent pinned source")
                continue
            original = git("show", f"{base}:{source}", optional=True)
            if original is None:
                issues.append(f"{prefix}: {source} must be committed BEFORE this test")
                continue
            excerpt = " ".join(metadata["ORACLE-EXCERPT"].split())
            if len(excerpt) < 15 or excerpt not in " ".join(original.split()):
                issues.append(f"{prefix}: exact oracle excerpt not found in {source}@{base}")
            if (len(metadata["ORACLE-LOCATOR"]) < 6
                    or len(metadata["ORACLE-HYPOTHESES"]) < 12
                    or len(metadata["ORACLE-EXPECTED"]) < 12
                    or len(metadata["ORACLE-DERIVATION"]) < 18):
                issues.append(f"{prefix}: locator, hypotheses, independent expectation or derivation is incomplete")
    return issues


def verify(args):
    if args.revision_range:
        problems = []
        for revision in git("rev-list", "--reverse", args.revision_range).splitlines():
            base = revision + "^"
            patch = git("diff", "--unified=0", base, revision, "--", "tests/")
            problems.extend(inspect(base, patch,
                                    lambda p, r=revision: git("show", f"{r}:{p}", optional=True),
                                    args.node))
        return problems
    if args.staged:
        patch = git("diff", "--cached", "--unified=0", "HEAD", "--", "tests/")
        return inspect("HEAD", patch, lambda p: git("show", f":{p}", optional=True), args.node)
    patch = git("diff", "--unified=0", "HEAD", "--", "tests/")
    for p in git("ls-files", "--others", "--exclude-standard", "--", "tests/").splitlines():
        if p.endswith((".py", ".sage")):
            lines = Path(p).read_text().splitlines()
            patch += f"\n+++ b/{p}\n@@ -0,0 +1,{len(lines)} @@\n"
            patch += "".join(f"+{line}\n" for line in lines)
    return inspect("HEAD", patch, lambda p: Path(p).read_text() if Path(p).is_file() else None,
                   args.node)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--staged", action="store_true")
    mode.add_argument("--worktree", action="store_true")
    mode.add_argument("--range", dest="revision_range")
    parser.add_argument("--node")
    args = parser.parse_args(argv)
    try:
        problems = verify(args)
    except (RuntimeError, ValueError, OSError) as error:
        problems = [f"cannot verify oracle provenance: {error}"]
    for problem in problems:
        print("ORACLE REQUIRED: " + problem, file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
