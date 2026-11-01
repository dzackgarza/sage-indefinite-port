"""Run a command and sample the peak resident memory of its process tree.

Purpose: the square [0] node was SIGKILLed at 373 s against a 7200 s cap with no
allocation error in the log. A bare exit 137 is not OOM evidence, so this records what
would be evidence: the peak summed RSS of the process tree, and the low-water mark of
the host's MemAvailable, sampled every SAMPLE_EVERY seconds.

Usage: python peakrss.py -- <command ...>
"""
import os
import subprocess
import sys
import time

SAMPLE_EVERY = 2.0
LABEL = os.environ.get("PEAKRSS_LABEL", "run")


def read_status(pid):
    try:
        with open(f"/proc/{pid}/status") as handle:
            fields = {}
            for line in handle:
                key, _, rest = line.partition(":")
                fields[key] = rest.strip()
            return int(fields["VmRSS"].split()[0]), int(fields["VmHWM"].split()[0])
    except (FileNotFoundError, ProcessLookupError, KeyError, PermissionError):
        return None


def descendants(root):
    children = {}
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(f"/proc/{entry}/stat") as handle:
                stat = handle.read()
            ppid = int(stat.rsplit(")", 1)[1].split()[1])
        except (FileNotFoundError, ProcessLookupError, IndexError, ValueError):
            continue
        children.setdefault(ppid, []).append(int(entry))
    seen, stack = {root}, [root]
    while stack:
        for child in children.get(stack.pop(), ()):
            if child not in seen:
                seen.add(child)
                stack.append(child)
    return seen


def mem_available_mib():
    try:
        with open("/proc/meminfo") as handle:
            for line in handle:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) // 1024
    except OSError:
        pass
    return -1


command = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else sys.argv[1:]
started_at = time.monotonic()
process = subprocess.Popen(command)
peak_tree, low_available, samples = 0, 10 ** 9, []
while process.poll() is None:
    pids = descendants(process.pid)
    tree = sum(filter(None, ((read_status(p) or (0, 0))[0] for p in pids)))
    hwm = sum(filter(None, ((read_status(p) or (0, 0))[1] for p in pids)))
    available = mem_available_mib()
    peak_tree = max(peak_tree, tree)
    low_available = min(low_available, available if available >= 0 else low_available)
    samples.append((round(time.monotonic() - started_at, 1), len(pids), tree // 1024, hwm // 1024, available))
    time.sleep(SAMPLE_EVERY)
status = process.wait()
START = samples[0][0] if samples else 0.0

shell = status if status >= 0 else 128 - status
print(f"[peakrss] {LABEL}: child waited={status} shell_code={shell} samples={len(samples)}")
print(f"[peakrss] peak tree RSS      : {peak_tree // 1024} MiB")
print(f"[peakrss] low-water MemAvail : {low_available} MiB")
peak_process = max((s[3] for s in samples), default=0)
print(f"[peakrss] peak single VmHWM  : {peak_process} MiB")
for sample in sorted(samples, key=lambda s: -s[2])[:5]:
    print(f"[peakrss]   t={sample[0]:8.1f}s procs={sample[1]:3d} rss={sample[2]:6d} MiB hwm={sample[3]:6d} MiB avail={sample[4]:6d} MiB")

# Propagate the child's status: a wrapper that swallows it turns a SIGKILL into a reported pass.
sys.exit(status if status >= 0 else 128 - status)
