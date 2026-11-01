#!/usr/bin/env python3
"""Run one acceptance file's nodes, one pytest invocation per node, under a memory guard and a cap.

One node per process is deliberate: it bounds peak memory, and it gives every node its own verdict instead
of letting one heavy case take the others down with it. Results are appended to the output log as they land,
so a long sweep leaves usable evidence even if it is stopped.
"""
import argparse
import os
import re
import signal
import subprocess
import sys
import time

HOME = os.path.expanduser("~")
IMAGE = os.path.join(HOME, "sage-image")
REPO = "/workspace/nucleus/work/projects/p-165a6922b36b40789ddc90db34bb23c7/repo"
CI = os.path.join(HOME, "ai-review-ci")
PROV = os.environ.get("SAGE_EXTENSIONS_TREE", os.path.join(HOME, "sage-extensions-cdafca84"))
PYTHONPATH = "/ci/tool-artifacts/pytest_plugins:/work/src:/work:/probe"


def proot(binds_extra=(), argv=()):
    args = [os.path.join(IMAGE, "run.sh"), "-w", "/work",
            "-b", f"{REPO}:/work", "-b", f"{CI}:/ci", "-b", f"{PROV}:/repos/sage-extensions"]
    for extra in binds_extra:
        args += ["-b", extra]
    args += ["--"]
    args += list(argv)
    return args


class Finished:
    def __init__(self, returncode, stdout, stderr):
        self.returncode, self.stdout, self.stderr = returncode, stdout, stderr


def run_bounded(cmd, env, cap, grace=60):
    """Run cmd, and enforce the cap on its whole process GROUP.

    `timeout` alone does not bound these runs: it kills the ptrace supervisor, the traced python survives,
    is reparented, and keeps the captured pipes open, so a naive subprocess.run blocks well past the cap.
    Own the session and kill the group.
    """
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                               env=env, cwd=IMAGE, start_new_session=True)
    try:
        out, err = process.communicate(timeout=cap + grace)
        return Finished(process.returncode, out, err)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(process.pid), signal.SIGKILL)
        except ProcessLookupError:
            pass
        out, err = process.communicate()
        return Finished(124, out, err + "\n[driver] cap enforced on the process group")


def collect(test_file):
    argv = ["/usr/bin/env", f"PYTHONPATH={PYTHONPATH}", "python", "-u", "-m", "pytest",
            "--collect-only", "-q", "-p", "qc_sage_session", "-o", "addopts=", test_file]
    result = subprocess.run(proot(argv=argv), capture_output=True, text=True, cwd=IMAGE)
    ids = [line.strip() for line in result.stdout.splitlines()
           if line.strip().startswith(test_file) and "::" in line]
    if not ids:
        sys.stderr.write(result.stdout[-3000:] + result.stderr[-3000:])
    return ids


def run_node(node, cap, guard, label, extra_binds, pytest_extra=()):
    inner = ["/usr/bin/env", f"PYTHONPATH={PYTHONPATH}", "python", "-u", "-m", "pytest", node,
             "-p", "qc_sage_session", "-q", "--no-header", "-p", "no:cacheprovider", "-o", "addopts=",
             *pytest_extra]
    cmd = ["python3", os.path.join(IMAGE, "peakrss.py"), "--",
           "timeout", "-s", "KILL", str(cap)] + proot(extra_binds, inner)
    env = dict(os.environ, PEAKRSS_LABEL=label, PEAKRSS_LIMIT_MIB=str(guard))
    started = time.monotonic()
    result = run_bounded(cmd, env, cap)
    elapsed = round(time.monotonic() - started)
    output = result.stdout + result.stderr
    peak = 0
    for match in re.finditer(r"peak tree RSS\s*:\s*(\d+) MiB", output):
        peak = int(match.group(1))
    verdict = "passed" if result.returncode == 0 else (
        "CAPPED" if result.returncode == 124 else f"exit{result.returncode}")
    failed = [line for line in output.splitlines() if line.startswith("FAILED")]
    return verdict, peak, elapsed, failed


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--file", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--cap", type=int, default=1800)
    parser.add_argument("--guard", type=int, default=6144)
    parser.add_argument("--out", required=True)
    parser.add_argument("--ids", default="")
    parser.add_argument("--pytest-extra", default="", help="comma-separated extra pytest arguments")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=0, help="stop after this many nodes (0 = all)")
    parser.add_argument("--stride", type=int, default=1, help="take every Nth node (a documented systematic sample)")
    parser.add_argument("--provider", default="", help="host path of the sage-extensions tree to bind")
    parser.add_argument("--whole-file", action="store_true",
                        help="one pytest invocation over the whole file (for parametrisations too large to run per node)")
    parser.add_argument("--skip", default="", help="comma-separated substrings; a node whose id contains one is skipped")
    parser.add_argument("--bind", action="append", default=[])
    arguments = parser.parse_args()

    global PROV
    if arguments.provider:
        PROV = arguments.provider
    handles = [open(path, "a", buffering=1) for path in (arguments.out, "/dev/stdout")]

    def emit(line):
        for handle in handles:
            handle.write(line + "\n")

    emit(f"=== {arguments.label} | file {arguments.file} | src {os.popen(f'git -C {REPO} rev-parse HEAD:src').read().strip()[:12]} "
         f"| provider {os.path.basename(PROV)} | cap {arguments.cap}s | guard {arguments.guard} MiB | nodes from {arguments.start}")

    if arguments.whole_file:
        inner = ["/usr/bin/env", f"PYTHONPATH={PYTHONPATH}", "python", "-u", "-m", "pytest", arguments.file,
                 "-p", "qc_sage_session", "-q", "--no-header", "-p", "no:cacheprovider", "-o", "addopts="]
        cmd = ["python3", os.path.join(IMAGE, "peakrss.py"), "--",
               "timeout", "-s", "KILL", str(arguments.cap)] + proot(arguments.bind, inner)
        env = dict(os.environ, PEAKRSS_LABEL=arguments.label, PEAKRSS_LIMIT_MIB=str(arguments.guard))
        started = time.monotonic()
        result = run_bounded(cmd, env, arguments.cap)
        elapsed = round(time.monotonic() - started)
        output = result.stdout + result.stderr
        peak = 0
        for match in re.finditer(r"peak tree RSS\s*:\s*(\d+) MiB", output):
            peak = int(match.group(1))
        summary = [line for line in output.splitlines()
                   if re.match(r"^=+ .*(passed|failed|error).* =+$", line.strip()) or line.startswith("FAILED")]
        emit(f"{arguments.label} WHOLE exit={result.returncode} peak={peak} MiB elapsed={elapsed}s")
        for line in summary[:40]:
            emit(f"    {line[:200]}")
        emit(f"{arguments.label}_DONE")
        return 0

    if arguments.ids:
        nodes = [node.strip() for node in arguments.ids.split(",") if node.strip()]
    else:
        nodes = collect(arguments.file)
    emit(f"nodes: {len(nodes)}")
    if not nodes:
        return 1

    skips = [piece for piece in arguments.skip.split(",") if piece]
    selected = list(enumerate(nodes, start=1))[arguments.start:]
    if arguments.stride > 1:
        selected = selected[::arguments.stride]
    for index, node in selected:
        if any(piece in node for piece in skips):
            emit(f"{arguments.label} {index}/{len(nodes)} {node} SKIPPED (excluded by --skip)")
            continue
        if arguments.limit and index - arguments.start > arguments.limit * arguments.stride:
            break
        verdict, peak, elapsed, failed = run_node(node, arguments.cap, arguments.guard,
                                                  f"{arguments.label}{index}", arguments.bind,
                                                  [piece for piece in arguments.pytest_extra.split(",") if piece])
        emit(f"{arguments.label} {index}/{len(nodes)} {node} exit={0 if verdict == 'passed' else verdict} "
             f"peak={peak} MiB verdict={1 if verdict == 'passed' else 'none'} {verdict}  {elapsed}s")
        for line in failed[:3]:
            emit(f"    {line[:200]}")
    emit(f"{arguments.label}_DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
