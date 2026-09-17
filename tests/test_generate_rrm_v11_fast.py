#!/usr/bin/env python3
"""Compare sequential output with v11; optionally time three fresh runs each.

Requires Python 3.8+ and GAP. Logs and generated files remain in a new directory
under /tmp. Additional preprocessed GRRM inputs can be passed with --input.
"""

import argparse
import filecmp
import json
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile
import time


ROOT = Path(__file__).resolve().parents[1]
SOURCES = {"v11": "generate_rrm_v11.g", "fast": "generate_rrm_v11_fast.g"}
LABELS = ["", "true", "false", "true,true", "true,false", "false,true", "false,false"]


def run_gap(args, out, tag, body, timeout=300):
    script = out / (tag + ".g")
    script.write_text(body + '\nPrint("RRM_TEST_DONE\\n");\nQUIT;\n')
    started = time.perf_counter()
    with (out / (tag + ".log")).open("w") as log:
        result = subprocess.run(
            [args.gap, "-T", "-b", "-q", "-r", "-m", args.memory, str(script)],
            cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=timeout, check=False,
        )
    elapsed = time.perf_counter() - started
    text = (out / (tag + ".log")).read_text()
    if result.returncode or "RRM_TEST_DONE" not in text.splitlines() or re.search(
        r"^(?:Syntax (?:error|warning)|Error,)", text, re.MULTILINE
    ):
        raise RuntimeError("GAP failed: {}\n{}".format(tag, text))
    return elapsed, text


def generate_body(source, gfile, calls):
    lines = ["Read({});".format(json.dumps(str(ROOT / SOURCES[source]))),
             "Read({});".format(json.dumps(str(gfile)))]
    for vertex, edge, labels in calls:
        lines.append("generate_rrm({},{},symc,ur,urt,ss,org_eq,org_ts{});".format(
            json.dumps(str(vertex)), json.dumps(str(edge)), "," + labels if labels else ""
        ))
    return "\n".join(lines)


def compare_pair(reference, candidate):
    for expected, actual in zip(reference, candidate):
        if not filecmp.cmp(expected, actual, shallow=False):
            raise AssertionError("Output differs: {} != {}".format(expected, actual))
    if not candidate[0].stat().st_size:
        raise AssertionError("Empty vertex output: {}".format(candidate[0]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gap", default=os.environ.get("GAP", "gap"))
    parser.add_argument("--memory", default="512m", help="GAP initial workspace (-m)")
    parser.add_argument("--input", action="append", type=Path, default=[])
    parser.add_argument("--benchmark", action="store_true")
    args = parser.parse_args()
    out = Path(tempfile.mkdtemp(prefix="rrm-sequential-test-", dir="/tmp"))
    print("Output: {}".format(out), flush=True)
    _, index_log = run_gap(args, out, "vertex_index",
                           'Read("tests/test_rrm_vertex_index.g");')
    print(index_log.strip(), flush=True)
    fixtures = [ROOT / "data/Au5Ag_AFIR.g", ROOT / "data/AuCu4_AFIR.g",
                ROOT / "tests/fixtures/small_labels.g", ROOT / "tests/fixtures/zero_ts.g",
                ROOT / "tests/fixtures/synthetic_n6.g"]
    fixtures += [path.resolve() for path in args.input]
    references = {}
    for case, gfile in enumerate(fixtures):
        outputs, logs = {}, {}
        for source in SOURCES:
            calls = [(out / "{}_{}_{}_v.dat".format(case, source, label),
                      out / "{}_{}_{}_e.dat".format(case, source, label), labels)
                     for label, labels in enumerate(LABELS)]
            _, logs[source] = run_gap(args, out, "{}_{}".format(case, source),
                                     generate_body(source, gfile, calls))
            outputs[source] = [call[:2] for call in calls]
        for reference, candidate in zip(outputs["v11"], outputs["fast"]):
            compare_pair(reference, candidate)
        # Compare warning identities as well as files: AuCu4 must still continue.
        warnings = [re.findall(r"^TS\d+\*?:$", logs[source], re.MULTILINE)
                    for source in SOURCES]
        if warnings[0] != warnings[1]:
            raise AssertionError("Pechukas diagnostics differ: {}".format(gfile))
        if case == 1 and warnings[1] != ["TS8:", "TS16:"] * len(LABELS):
            raise AssertionError("Expected both AuCu4 Pechukas diagnostics")
        if case == 3 and any(pair[1].stat().st_size for pair in outputs["fast"]):
            raise AssertionError("Zero-TS input produced edges")
        references[gfile] = outputs["v11"][3]  # Both labels explicitly enabled.
        print("PASS: {} (7 label forms, byte-identical)".format(gfile), flush=True)

    # A corrupt/missing index lookup must not silently emit 'fail' as a vertex ID.
    bad_script = out / "missing_vertex.g"
    bad_script.write_text('Read("generate_rrm_v11_fast.g");\n'
                          'G:=SymmetricGroup(4);; U:=TrivialSubgroup(G);;\n'
                          'b:=RrmBuildTransversals(G,[U]);;\n'
                          'RrmVertexIndex(b.offset,b.idx,[U],1,(1,2,3,4,5));\n'
                          'Print("UNREACHABLE\\n");\nQUIT;\n')
    bad = subprocess.run([args.gap, "-T", "-b", "-q", "-r", str(bad_script)],
                         cwd=ROOT, capture_output=True, text=True, timeout=30)
    (out / "missing_vertex.log").write_text(bad.stdout + bad.stderr)
    if bad.returncode != 1 or "UNREACHABLE" in bad.stdout or "vertex lookup failed" not in bad.stdout:
        raise AssertionError("Missing vertex lookup did not stop with its diagnostic")
    print("PASS: missing vertex lookup stops", flush=True)

    measurements = []
    if args.benchmark:
        for case in (0, 4):
            gfile = fixtures[case]
            times = {source: [] for source in SOURCES}
            for repeat in range(3):
                # Alternate the order to reduce systematic ordering effects.
                order = list(SOURCES) if repeat % 2 == 0 else list(reversed(SOURCES))
                for source in order:
                    tag = "bench_{}_{}_{}".format(case, repeat, source)
                    pair = (out / (tag + "_v.dat"), out / (tag + "_e.dat"))
                    wall, _ = run_gap(args, out, tag, generate_body(
                        source, gfile, [(*pair, "true,true")]))
                    compare_pair(references[gfile], pair)
                    times[source].append(wall)
            speedup = statistics.median(times["v11"]) / statistics.median(times["fast"])
            measurements.append({"input": str(gfile), "wall_seconds": times, "speedup": speedup})
            print("BENCH: {} v11={:.3f}s fast={:.3f}s speedup={:.2f}x".format(
                gfile.name, statistics.median(times["v11"]),
                statistics.median(times["fast"]), speedup), flush=True)
        (out / "timings.json").write_text(json.dumps(measurements, indent=2) + "\n")
    print("PASS: {} file pairs; logs and outputs: {}".format(len(fixtures) * 7, out))


if __name__ == "__main__":
    main()
