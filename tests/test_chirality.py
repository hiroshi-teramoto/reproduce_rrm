"""Run with the rrm-recon Python environment on a compute node."""
import ast
import math
from pathlib import Path

import numpy as np
from pymatgen.core import Molecule
from pymatgen.symmetry.analyzer import PointGroupAnalyzer


def main():
    source = Path(__file__).resolve().parents[1] / "rrm_reconstruction_v18.py"
    # Load the helpers without executing the script's command-line file reads.
    tree = ast.parse(source.read_text())
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef)
                 and n.name in {"dist", "perm", "wchiral", "rotation_permutations"}]
    namespace = {"np": np, "math": math, "tol": 0.1}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), namespace)
    cases = [
        ("HCN", ["H", "C", "N"], [[0, 0, -1], [0, 0, 0], [0, 0, 1.2]], "C*v", False),
        ("CCO", ["C", "C", "O"], [[-1.3, 0, 0], [0, 0, 0], [1.2, 0, 0]], "C*v", False),
        ("OCO", ["O", "C", "O"], [[0, -1.2, 0], [0, 0, 0], [0, 1.2, 0]], "D*h", False),
        ("water", ["O", "H", "H"], [[0, 0, 0], [1, 0, 0], [0, 1, 0]], "C2v", False),
        ("CHFClBr", ["C", "H", "F", "Cl", "Br"],
         [[0, 0, 0], [1, 1, 1], [-1, -1, 1], [-1, 1, -1], [1, -1, -1]], "C1", True),
    ]
    for name, species, coordinates, point_group, expected in cases:
        # Check the original and its inversion, including a truly chiral control.
        for sign in (1, -1):
            molecule = Molecule(species, sign * np.array(coordinates)).get_centered_molecule()
            analyzer = PointGroupAnalyzer(molecule, tolerance=0.1)
            assert analyzer.sch_symbol == point_group, (name, analyzer.sch_symbol)
            actual = namespace["wchiral"](analyzer.get_symmetry_operations(), analyzer.sch_symbol)
            assert actual is expected, (name, sign, actual)
            permutations = {tuple(p) for p in namespace["rotation_permutations"](
                molecule, analyzer.get_symmetry_operations(), analyzer.sch_symbol)}
            expected_permutations = {tuple(range(1, len(species) + 1))}
            if name == "OCO":
                expected_permutations.add((3, 2, 1))  # (1,3) must be in ur/urt.
            elif name == "water":
                expected_permutations.add((1, 3, 2))
            assert permutations == expected_permutations, (name, sign, permutations)

    # Nonlinear Td must still exclude improper operations: 12 even permutations,
    # rather than all 24 permutations of the four identical ligands.
    methane = Molecule(["C", "H", "H", "H", "H"], cases[-1][2])
    analyzer = PointGroupAnalyzer(methane, tolerance=0.1)
    assert analyzer.sch_symbol == "Td"
    permutations = {tuple(p) for p in namespace["rotation_permutations"](
        methane, analyzer.get_symmetry_operations(), analyzer.sch_symbol)}
    assert len(permutations) == 12, permutations
    for p in permutations:
        assert p[0] == 1
        assert sum(p[i] > p[j] for i in range(5) for j in range(i + 1, 5)) % 2 == 0, p
    print("PASS: chirality and rotation permutations, including D*h endpoint exchange and nonlinear Td")


if __name__ == "__main__":
    main()
