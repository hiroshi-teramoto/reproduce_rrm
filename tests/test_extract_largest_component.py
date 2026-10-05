#!/usr/bin/env python3
"""Run with python3 tests/test_extract_largest_component.py (stdlib only)."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

SCRIPT = Path(__file__).resolve().parents[1] / "extract_largest_component.py"


def main():
    with tempfile.TemporaryDirectory(prefix="rrm-component-test-") as directory:
        root = Path(directory)
        prefix = root / "molecule"
        eq_path = root / "molecule_EQ_list.log"
        ts_path = root / "molecule_TS_list.log"
        payload = b"H   0.000000000001  1.0  2.0\r\nEnergy = -123.456789012345\r\n\r\n"

        def eq(number):
            return f"# Geometry of EQ {number}, SYMMETRY = C1\r\n".encode() + payload

        def ts(number, a, b):
            return (f"# Geometry of TS {number}, SYMMETRY = C1\r\n".encode() + payload
                    + f"CONNECTION : {a} - {b}\r\n\r\n".encode())

        # Unsorted, gapped IDs; parallel edges, retained and excluded self loops.
        original_eq = b"List of Equilibrium Structures\r\n\r\n" + b"".join(eq(i) for i in [90, 20, 7, 100])
        original_ts = b"List of Transition Structures\r\n\r\n" + b"".join([
            ts(80, 90, 90), ts(12, 7, 7), ts(10, 20, 7), ts(4, 7, 20)])
        eq_path.write_bytes(original_eq)
        ts_path.write_bytes(original_ts)
        # Missing TS80 is fine: that component is excluded.
        for i in [4, 10, 12]:
            (root / f"molecule_TS{i}.log").write_bytes(f"Original TS{i}\r\n".encode() + payload)

        def run(name, expected_error=None):
            output = root / name
            result = subprocess.run([sys.executable, str(SCRIPT), str(prefix), "--out", str(output)],
                                    text=True, capture_output=True)
            if expected_error is None:
                assert result.returncode == 0, result.stderr
            else:
                assert result.returncode != 0 and expected_error in result.stderr, result.stderr
                assert not output.exists(), name
            return output

        output = run("selected")
        data = json.loads((output / "mapping.json").read_text())
        assert data["eq_new_to_old"] == {"0": 7, "1": 20}
        assert data["ts_new_to_old"] == {"0": 4, "1": 10, "2": 12}
        assert data["input_counts"] == {"eq": 4, "ts": 4}
        assert data["output_counts"] == {"eq": 2, "ts": 3}
        assert data["excluded_components"] == [
            {"eq_ids": [90], "ts_ids": [80], "eq_count": 1, "ts_count": 1},
            {"eq_ids": [100], "ts_ids": [], "eq_count": 1, "ts_count": 0}]
        assert (output / "EQ_list.log").read_bytes() == b"List of Equilibrium Structures\r\n\r\n" + eq(0) + eq(1)
        assert (output / "TS_list.log").read_bytes() == b"List of Transition Structures\r\n\r\n" + ts(0, 0, 1) + ts(1, 1, 0) + ts(2, 0, 0)
        for new, old in enumerate([4, 10, 12]):
            assert (output / f"TS{new}.log").read_bytes() == (root / f"molecule_TS{old}.log").read_bytes()
        assert eq_path.read_bytes() == original_eq and ts_path.read_bytes() == original_ts

        retained_log = root / "molecule_TS4.log"
        retained_bytes = retained_log.read_bytes()
        retained_log.unlink()
        run("missing", "missing or empty")
        retained_log.write_bytes(b"")
        run("empty", "missing or empty")
        retained_log.write_bytes(retained_bytes)
        for name, data, error in [
            ("dc", original_ts.replace(b"20 - 7", b"DC0 - 7"), "numeric CONNECTION"),
            ("unknown", original_ts.replace(b"20 - 7", b"999 - 7"), "missing EQ"),
            ("duplicate_ts", original_ts + ts(4, 7, 20), "duplicate TS4"),
            ("double_connection", original_ts + b"CONNECTION : 7 - 20\n", "numeric CONNECTION"),
            ("bad_header", original_ts.replace(b"Geometry of TS 4,", b"Geometry of TS x,"), "geometry header"),
        ]:
            ts_path.write_bytes(data)
            run(name, error)
        ts_path.write_bytes(original_ts)
        eq_path.write_bytes(original_eq + eq(7))
        run("duplicate_eq", "duplicate EQ7")
        eq_path.write_bytes(original_eq)
        ts_path.write_bytes(original_ts + ts(81, 90, 100))
        run("tie", "tied")
        ts_path.write_bytes(original_ts)
        existing = root / "existing"
        existing.mkdir()
        (existing / "sentinel").write_text("keep")
        result = subprocess.run([sys.executable, str(SCRIPT), str(prefix), "--out", str(existing)], capture_output=True)
        assert result.returncode != 0 and (existing / "sentinel").read_text() == "keep"

        # Sorting must not concatenate headers onto an unterminated last block.
        eq_path.write_bytes(eq(20) + eq(7).rstrip(b"\r\n"))
        ts_path.write_bytes(ts(10, 20, 7) + ts(4, 7, 20).rstrip(b"\r\n"))
        eof_output = run("no_final_newline")
        assert (eof_output / "EQ_list.log").read_bytes() == eq(0).rstrip(b"\r\n") + b"\r\n" + eq(1)
        assert (eof_output / "TS_list.log").read_bytes() == ts(0, 0, 1).rstrip(b"\r\n") + b"\r\n" + ts(1, 1, 0)

        # A graph that is already connected, including the one-EQ / zero-TS case.
        eq_path.write_bytes(eq(7) + eq(20))
        ts_path.write_bytes(ts(4, 7, 20))
        assert json.loads((run("connected") / "mapping.json").read_text())["excluded_components"] == []
        eq_path.write_bytes(eq(7))
        ts_path.write_bytes(b"List of Transition Structures\n")
        assert json.loads((run("single") / "mapping.json").read_text())["output_counts"] == {"eq": 1, "ts": 0}
    print("PASS: extraction, numbering, bytes, loops, parallel edges, validation, and overwrite protection")


if __name__ == "__main__":
    main()
