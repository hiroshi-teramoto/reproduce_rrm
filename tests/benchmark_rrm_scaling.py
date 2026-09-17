#!/usr/bin/env python3
"""Scale identical-atom count with eight EQs and 30 TSs held fixed.

The synthetic S6 fixture describes six identical atoms plus a distinct atom at
point 1. Extend its ambient group to S7/S8, keeping the same EQ/TS subgroups and
endpoint permutations (new points are fixed by those permutations). This is a
controlled group-theory workload, not a molecular geometry or GRRM calculation.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import statistics
import subprocess
import tempfile

from test_generate_rrm_v11_fast import ROOT, SOURCES, compare_pair, generate_body, run_gap


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gap", default=os.environ.get("GAP", "gap"))
    parser.add_argument("--memory", default="512m")
    parser.add_argument("--timeout", type=int, default=300, help="seconds per GAP process")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    out = Path(tempfile.mkdtemp(prefix="rrm-scaling-", dir="/tmp"))
    print("Output: {}".format(out), flush=True)
    template_path = ROOT / "tests/fixtures/synthetic_n6.g"
    template = template_path.read_text()
    template = template[template.index("sym:="):]
    report = {"memory": args.memory, "timeout_seconds": args.timeout,
              "repetitions": 3, "labels": [True, True],
              "source_sha256": {name: sha256(ROOT / name) for name in SOURCES.values()},
              "template_sha256": sha256(template_path), "cases": []}
    for identical in (6, 7, 8):
        gfile = out / "synthetic_n{}.g".format(identical)
        gfile.write_text("# Synthetic workload; one distinct atom plus {} identical atoms.\n".format(identical)
                        + template.replace("SymmetricGroup([2..7])",
                                           "SymmetricGroup([2..{}])".format(identical + 1)))
        factor = math.factorial(identical) // math.factorial(6)
        case = {"identical_atoms": identical, "total_atoms": identical + 1,
                "eqs": 8, "ts": 30, "vertices": 4440 * factor, "edges": 15600 * factor,
                "input_sha256": sha256(gfile), "runs": [], "reference_compared": False}
        report["cases"].append(case)
        outputs, times, timed_out = {}, {name: [] for name in SOURCES}, set()
        for repeat in range(3):
            order = ("fast", "v11") if repeat % 2 == 0 else ("v11", "fast")
            for source in order:
                if source in timed_out:
                    continue
                tag = "n{}_{}_{}".format(identical, source, repeat)
                pair = (out / (tag + "_v.dat"), out / (tag + "_e.dat"))
                run = {"source": source, "repeat": repeat + 1}
                try:
                    wall, log = run_gap(args, out, tag, generate_body(
                        source, gfile, [(*pair, "true,true")]), timeout=args.timeout)
                except subprocess.TimeoutExpired:
                    timed_out.add(source)
                    run.update(status="timeout", timeout_seconds=args.timeout)
                    print("TIMEOUT: {} after {} s; remaining repetitions skipped".format(
                        tag, args.timeout), flush=True)
                else:
                    if "Violation of Pechukus theorem" in log:
                        raise AssertionError("Unexpected Pechukas violation: " + tag)
                    with pair[0].open() as stream:
                        vertices = sum(bool(re.match(r"^\d+\[", line)) for line in stream)
                    with pair[1].open() as stream:
                        edges = sum(1 for line in stream)
                    if (vertices, edges) != (case["vertices"], case["edges"]):
                        raise AssertionError("Incorrect output counts: " + tag)
                    for previous in outputs.values():
                        compare_pair(previous, pair)
                    outputs[source] = pair
                    case["reference_compared"] = len(outputs) == 2
                    times[source].append(wall)
                    run.update(status="ok", wall_seconds=wall,
                               vertices_sha256=sha256(pair[0]), edges_sha256=sha256(pair[1]))
                    print("PASS: {} {:.3f} s, {} vertices / {} edges".format(
                        tag, wall, vertices, edges), flush=True)
                case["runs"].append(run)
                (out / "timings.json").write_text(json.dumps(report, indent=2) + "\n")
        case["median_wall_seconds"] = {
            source: statistics.median(values) for source, values in times.items() if values}
        if all(len(times[source]) == 3 for source in SOURCES):
            case["speedup"] = statistics.median(times["v11"]) / statistics.median(times["fast"])
        (out / "timings.json").write_text(json.dumps(report, indent=2) + "\n")
        print("CASE: {}".format(json.dumps({key: value for key, value in case.items()
                                          if key != "runs"})), flush=True)
    print("Finished; inspect completion and reference_compared per case: {}".format(out))


if __name__ == "__main__":
    main()
