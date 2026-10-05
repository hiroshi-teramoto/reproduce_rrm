#!/usr/bin/env python3
"""Extract the unique largest GRRM EQ/TS component for rrm_reconstruction_v18.py.

Uses only Python 3.8+ standard libraries; molecular data stay intact.
"""
import argparse
import json
from pathlib import Path
import re
import shutil


GEOMETRY = re.compile(
    rb"^([ \t]*#[ \t]*Geometry of (EQ|TS)[ \t]+)([0-9]+)(?=[ \t,])", re.M
)
CONNECTION = re.compile(
    rb"^([ \t]*CONNECTION[ \t]*:[ \t]*)([0-9]+)([ \t]*-[ \t]*)"
    rb"([0-9]+)([ \t]*\r?$)", re.M
)


def read_blocks(path, kind):
    data = path.read_bytes()
    headers = list(GEOMETRY.finditer(data))
    if len(headers) != len(re.findall(rb"^[ \t]*#[ \t]*Geometry\b", data, re.M)):
        raise ValueError(f"{path}: malformed geometry header")
    blocks = {}
    for index, header in enumerate(headers):
        number = int(header[3])
        if header[2] != kind.encode():
            raise ValueError(f"{path}: expected {kind} geometry headers")
        if number in blocks:
            raise ValueError(f"{path}: duplicate {kind}{number}")
        end = headers[index + 1].start() if index + 1 < len(headers) else len(data)
        blocks[number] = data[header.start():end]
    if not blocks and kind == "EQ":
        raise ValueError(f"{path}: no EQ structures")
    return data[:headers[0].start()] if headers else data, blocks


def read_connections(eq_blocks, ts_blocks):
    connections = {}
    for number, block in ts_blocks.items():
        matches = list(CONNECTION.finditer(block))
        if len(matches) != 1 or len(re.findall(rb"^[ \t]*CONNECTION\b", block, re.M)) != 1:
            raise ValueError(f"TS{number}: expected one numeric CONNECTION; DC/unknown connections are unsupported")
        a, b = int(matches[0][2]), int(matches[0][4])
        if a not in eq_blocks or b not in eq_blocks:
            raise ValueError(f"TS{number}: CONNECTION references missing EQ ({a}, {b})")
        connections[number] = (a, b)
    return connections


def connected_components(eq_blocks, connections):
    adjacency = {number: set() for number in eq_blocks}
    for a, b in connections.values():
        adjacency[a].add(b)
        adjacency[b].add(a)
    unseen, components = set(adjacency), []
    while unseen:
        node = unseen.pop()
        component, todo = {node}, [node]
        while todo:
            neighbors = adjacency[todo.pop()] & unseen
            unseen.difference_update(neighbors)
            component.update(neighbors)
            todo.extend(neighbors)
        components.append(sorted(component))
    return sorted(components, key=lambda nodes: (-len(nodes), nodes[0]))


def renumber(block, number, eq_map=None):
    header = GEOMETRY.match(block)
    block = block[:header.start(3)] + str(number).encode() + block[header.end(3):]
    if eq_map is not None:
        block = CONNECTION.sub(
            lambda m: m[1] + str(eq_map[int(m[2])]).encode() + m[3]
            + str(eq_map[int(m[4])]).encode() + m[5], block, count=1
        )
    return block


def join_blocks(header, blocks):
    # A final source block may lack a newline and move earlier after sorting.
    return header + b"".join(
        block if block.endswith(b"\n") else block + (b"\r\n" if b"\r\n" in block else b"\n")
        for block in blocks
    )


def extract(prefix, output):
    prefix, output = Path(prefix).resolve(), Path(output).absolute()
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Refusing to overwrite existing output: {output}")
    eq_header, eq_blocks = read_blocks(Path(str(prefix) + "_EQ_list.log"), "EQ")
    ts_header, ts_blocks = read_blocks(Path(str(prefix) + "_TS_list.log"), "TS")
    connections = read_connections(eq_blocks, ts_blocks)
    components = connected_components(eq_blocks, connections)
    if len(components) > 1 and len(components[0]) == len(components[1]):
        raise ValueError("Largest component is tied by EQ count; selection is ambiguous")
    eq_component = {eq: index for index, nodes in enumerate(components) for eq in nodes}
    component_ts = [[] for _ in components]
    for ts in sorted(connections):
        component_ts[eq_component[connections[ts][0]]].append(ts)
    kept_eq, kept_ts = components[0], component_ts[0]
    eq_map = {old: new for new, old in enumerate(kept_eq)}
    ts_map = {old: new for new, old in enumerate(kept_ts)}
    ts_logs = {old: Path(str(prefix) + f"_TS{old}.log") for old in kept_ts}
    for old, path in ts_logs.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"Retained TS{old}: missing or empty path log: {path}")
    mapping = {
        "source_prefix": str(prefix),
        "input_counts": {"eq": len(eq_blocks), "ts": len(ts_blocks)},
        "output_counts": {"eq": len(kept_eq), "ts": len(kept_ts)},
        "eq_new_to_old": {str(new): old for old, new in eq_map.items()},
        "ts_new_to_old": {str(new): old for old, new in ts_map.items()},
        "excluded_components": [
            {"eq_ids": nodes, "ts_ids": component_ts[index],
             "eq_count": len(nodes), "ts_count": len(component_ts[index])}
            for index, nodes in enumerate(components) if index != 0
        ],
    }

    # mkdir reserves a fresh destination; cleanup owns only this directory,
    # not its potentially shared parents (which may remain after a failure).
    output.mkdir(parents=True)
    try:
        (output / "EQ_list.log").write_bytes(
            join_blocks(eq_header, (renumber(eq_blocks[old], new) for old, new in eq_map.items()))
        )
        (output / "TS_list.log").write_bytes(
            join_blocks(ts_header, (renumber(ts_blocks[old], new, eq_map) for old, new in ts_map.items()))
        )
        for old, new in ts_map.items():
            shutil.copyfile(ts_logs[old], output / f"TS{new}.log")

        _, out_eq = read_blocks(output / "EQ_list.log", "EQ")
        _, out_ts = read_blocks(output / "TS_list.log", "TS")
        out_connections = read_connections(out_eq, out_ts)
        expected = {new: tuple(eq_map[eq] for eq in connections[old]) for old, new in ts_map.items()}
        if (list(out_eq) != list(range(len(kept_eq)))
                or list(out_ts) != list(range(len(kept_ts)))
                or out_connections != expected
                or len(connected_components(out_eq, out_connections)) != 1):
            raise ValueError("Extracted graph failed connectivity/renumbering verification")
        (output / "mapping.json").write_text(json.dumps(mapping, indent=2) + "\n", encoding="utf-8")
    except Exception:
        shutil.rmtree(output)
        raise
    return mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", type=Path, help="GRRM prefix, e.g. /data/Pt5_AFIR")
    parser.add_argument("--out", required=True, type=Path, help="New output directory (must not exist)")
    args = parser.parse_args()
    try:
        result = extract(args.prefix, args.out)
    except (OSError, ValueError) as error:
        parser.exit(1, f"error: {error}\n")
    before, after = result["input_counts"], result["output_counts"]
    print(f"{before['eq']} EQ / {before['ts']} TS -> {after['eq']} EQ / {after['ts']} TS")
    print(f"Excluded {len(result['excluded_components'])} components; output: {args.out}")


if __name__ == "__main__":
    main()
