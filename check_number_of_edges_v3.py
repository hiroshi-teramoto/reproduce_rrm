#!/usr/bin/env python3
import re
import sys
from collections import defaultdict

def parse_dot_file(filename):
    with open(filename, 'r') as f:
        content = f.read()

    # Split the file into statements using ";" as a separator.
    statements = re.split(r'\n', content)

    # Only process statements that do not contain an edge definition ("--")
    vertex_pattern = re.compile(r'^\s*(\d+)\s*\[label="([^"]+)"\]')
    vertices = {}
    for stmt in statements:
        if '--' in stmt:
            continue
        m = vertex_pattern.search(stmt)
        if m:
            vid = m.group(1)
            label = m.group(2)
            # The number is assumed to be the first token of the label.
            number = label.split()[0]
            vertices[vid] = {'label': label, 'num': number, 'degree': 0}

    # Parse edges in the full content.
    edge_pattern = re.compile(r'(\d+)\s*--\s*(\d+)')
    for match in edge_pattern.finditer(content):
        v1, v2 = match.group(1), match.group(2)
        if v1 == v2:
            if v1 in vertices:
                vertices[v1]["degree"] += 1
        else:
            if v1 in vertices:
                vertices[v1]['degree'] += 1
            if v2 in vertices:
                vertices[v2]['degree'] += 1
    return vertices

def check_consistency(vertices):
    # Group vertices by the number part of their label.
    groups = defaultdict(list)
    for vid, info in vertices.items():
        groups[info['num']].append((vid, info['degree']))

    # Check if each group has consistent degree counts.
    errors = {}
    for num, verts in groups.items():
        degs = set(deg for vid, deg in verts)
        if len(degs) > 1:
            errors[num] = verts
    return errors

def main():
    if len(sys.argv) != 2:
        print("Usage: {} <dot_file>".format(sys.argv[0]))
        sys.exit(1)

    filename = sys.argv[1]
    vertices = parse_dot_file(filename)
    errors = check_consistency(vertices)
    if errors:
        print("Error: Inconsistent degrees for vertices with the same number found:")
        for num, verts in errors.items():
            vert_info = ", ".join("vertex {} (degree {})".format(vid, deg) for vid, deg in verts)
            print("Number {}: {}".format(num, vert_info))
        sys.exit(1)
    else:
        print("All vertices with the same number have consistent degrees.")

if __name__ == '__main__':
    main()
