"""Execute a pytest shard in fresh Sage interpreters without dropping test nodes."""
import argparse
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sage', required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--batch-size', type=int, default=80)
    p.add_argument('pytest_args', nargs=argparse.REMAINDER)
    a = p.parse_args()
    if a.batch_size < 1:
        p.error('batch-size must be positive')
    args = a.pytest_args
    if args and args[0] == '--':
        args = args[1:]
    prefix = [a.sage, '-python', '-c', 'import sage.all,sys,pytest;sys.exit(pytest.main(sys.argv[1:]))', '-p', 'qc_sage_session']
    collected = subprocess.run(prefix + args + ['--collect-only', '-q'], capture_output=True, text=True)
    if collected.returncode not in (0, 5):
        sys.stdout.write(collected.stdout + collected.stderr)
        return collected.returncode
    # pytest --collect-only -q emits one nodeid per line; exclude summary lines.
    nodes = [line for line in collected.stdout.splitlines() if line.startswith('tests/') and '::' in line]
    if not nodes:
        print('No selected test nodes', file=sys.stderr)
        return 1
    print(f'Executing {len(nodes)} collected nodes in fresh batches of {a.batch_size}', flush=True)
    suite = ET.Element('testsuites')
    exitcode = 0
    for offset in range(0, len(nodes), a.batch_size):
        batch = nodes[offset:offset+a.batch_size]
        report = a.out / f'batch-{offset//a.batch_size}.xml'
        command = prefix + batch + ['--runxfail', '-rA', '--tb=long', '--instafail', f'--junitxml={report}']
        result = subprocess.run(command)
        if report.exists():
            root = ET.parse(report).getroot()
            if root.tag == 'testsuite':
                suite.append(root)
            else:
                suite.extend(root.findall('testsuite'))
        if result.returncode:
            exitcode = result.returncode
            if not report.exists():
                print(f'Batch {offset//a.batch_size} terminated without JUnit report', file=sys.stderr)
                break
    ET.ElementTree(suite).write(a.out / 'junit.xml', encoding='unicode', xml_declaration=True)
    return exitcode


if __name__ == '__main__':
    sys.exit(main())
