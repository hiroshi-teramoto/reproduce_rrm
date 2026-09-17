# Sequential v11 speedup: validation and measurements

Measured on qc4 on 2026-09-17 against upstream commit
`d64568d003788815e6d0ce501c643c1050e3f83c`.

## Implementation

`generate_rrm_v11_fast.g` is an optional implementation of the existing
`generate_rrm` function. For each EQ, a dictionary maps its canonical right-coset
representatives to their positions in the original `RightTransversal` order.
Adding the preceding EQs' vertex counts gives the same ID as v11's
`Position(vertices, [eq, representative])`. This avoids scanning the vertex list
for every edge endpoint and every vertex written.

The same-EQ Pechukas check acts on the cached transversal and reuses one action
homomorphism per EQ. Edges are written in the original TS/transversal order
through an open stream, without keeping the complete edge list in memory.
Vertex output also uses an open stream. The API, label defaults and Pechukas
diagnostic/continuation policy follow v11. Missing vertex lookups or failures to
open output files terminate with status 1 rather than producing invalid IDs.

## Reproduction

From the repository root, with GAP 4.13.0 and Python 3.8+:

```bash
python3 tests/test_generate_rrm_v11_fast.py --gap /path/to/gap --benchmark
```

The default initial GAP workspace is `512m`; override it with `--memory`.
The runner has no third-party Python dependencies. It writes GAP scripts, logs,
dat files and `timings.json` to a new directory under `/tmp`, printed at startup.
The regression runner gives each GAP invocation a 300-second timeout; a failure
or timeout fails that run. The separate scaling benchmark below records timeouts.

## Correctness

The supplied suite passed 72 vertex-index comparisons and 35 byte comparisons
of vertex/edge file pairs: each of these five fixtures with all seven forms of
the optional label arguments (omitted, one Boolean, or two Booleans).

| Fixture | Vertices | Edges | Coverage |
|---|---:|---:|---|
| Au5Ag | 1,704 | 10,020 | Bundled GRRM sample |
| AuCu4 | 102 | 474 | TS8 and TS16 diagnostics followed by continued output |
| Small labeled graph | 15 | 30 | Inversion labels, self-loops and multiple edges |
| Zero TS | 6 | 0 | Empty edge output |
| Synthetic S6 | 4,440 | 15,600 | Larger vertex index; synthetic, not a chemical model |

The suite also verifies that an invalid vertex lookup exits with status 1.
Every timed output below passed another byte comparison against v11.

An additional local run used `--input` for five preprocessed organic inputs:
CH4, CH5+, HCOOH, C3NH and a saved n-butane dataset. All seven label forms matched
for each input, bringing the correctness run to **70 file pairs / 140 files**.
CH5+ emitted the same Pechukas diagnostic and continued in both implementations,
without a continuation flag. Agreement with v11 on a flagged input does not
establish the physical validity of its map. These additional organic inputs are
not distributed in this PR; their source and preprocessing are recorded in the
[development fork's report](https://github.com/mochimochi-brian/reproduce_rrm/blob/55d8cc6/docs/results/organic-v11-byte-comparison.md).

## Performance

Environment: AMD EPYC 9654 (96 available CPUs; one GAP process at a time),
Linux 6.8.0-139-generic, GAP 4.13.0, Python 3.11.3. Source on NFS; outputs on
local XFS under `/tmp`. Both labels enabled, GAP `-T -b -q -r -m 512m`.
Three fresh processes per implementation and input, with alternating order.
Wall time includes process startup, input loading, generation and exit, measured
by Python's monotonic `perf_counter`. Preprocessing and post-run comparisons are
outside the timed interval. Measurements follow the correctness suite, so this
is a warm-cache comparison.

| Input | v11 wall times (s) | Fast wall times (s) | Median v11 / fast (s) | Speedup |
|---|---|---|---|---:|
| Au5Ag | 4.424, 4.475, 4.474 | 1.618, 1.618, 1.618 | 4.474 / 1.618 | 2.77x |
| Synthetic S6 | 2.770, 2.720, 2.720 | 1.367, 1.367, 1.367 | 2.720 / 1.367 | 1.99x |

These are measurements for these inputs and this environment, not a scaling
guarantee. Peak memory was not measured in this extraction.

## Controlled atom-count scaling

The scaling benchmark extends the existing synthetic fixture from S6 to S7 and
S8. The model has one distinct atom at point 1 plus 6, 7 or 8 identical atoms,
so the total atom counts are 7, 8 and 9. The eight EQ stabilizers, 30 TS
stabilizers and endpoint permutations are held fixed; their permutations fix the
additional points. Only the ambient group `sym = symc = S_n` grows. Inversion
symmetry is not added. This is a controlled group-theory workload, not a GRRM
calculation on a molecular geometry. Real expansion sizes also depend on element
composition, stabilizers and the connected component being reconstructed.

For fixed stabilizers, increasing n from 6 to 7 multiplies both output counts by
7, and increasing n from 7 to 8 multiplies them again by 8. The runner checks
these expected counts for every completed run. Both labels are enabled, with
the same GAP version, machine, initial workspace and output filesystem as above.
Three runs are requested for each implementation, alternating order and starting
with the fast implementation. No cache flush or separate warmup is performed.
Completed outputs are compared byte for byte against all earlier completed
outputs for that input, including v11 whenever a reference finishes.

```bash
python3 tests/benchmark_rrm_scaling.py --gap /path/to/gap --memory 512m --timeout 300
```

The timeout is per process. After an implementation times out, its remaining
repetitions are skipped; an unfinished run is not included in a median or
reported as a completed runtime. The output directory contains the exact
extended inputs, generated GAP commands, logs, dat files and `timings.json`.

| Total atoms (identical + distinct) | Vertices | Edges | v11 median wall (s) | Fast median wall (s) | Speedup | Byte comparison with v11 |
|---|---:|---:|---:|---:|---:|---|
| 7 (6 + 1) | 4,440 | 15,600 | 2.720 | 1.367 | 1.99x | Passed, all runs |
| 8 (7 + 1) | 31,080 | 109,200 | 57.497 | 2.520 | 22.82x | Passed, all runs |
| 9 (8 + 1) | 248,640 | 873,600 | Timed out at 300 s | 11.392 | Not calculated | Unavailable (timeout) |

The eight-atom case demonstrates a **22.82x** measured speedup with complete
reference agreement. In the nine-atom case, only the fast implementation
completed: 11.392, 12.091 and 11.391 seconds (median 11.392 seconds). Its three
outputs are byte-identical to each other and have the expected 248,640 vertices
and 873,600 edges. The original v11 attempt timed out at 300 seconds and was
not repeated. There is no complete reference comparison or measured speedup
for that case.

[Per-run timings, statuses and SHA-256 digests](data/sequential-atom-scaling.json)
include every attempt; timed-out output has no completed-file checksum.
Raw scaling output on qc4: `/tmp/rrm-scaling-7_aij_y6/`.

SHA-256 of the measured scripts:

```text
a61bf68fd8253306561c283559dc2c484ec27ee76c6e6770301eba0eda7f10d5  generate_rrm_v11.g
1abdfabe2c21faca76d55e022cd1d63db041c22e6280402556a4a87a2b294265  generate_rrm_v11_fast.g
```

Raw logs, outputs and timings on qc4: `/tmp/rrm-sequential-test-fo72id8s/`.
