#!/usr/bin/env python3
"""Checks the ParaView collections written by the polycrystal.

The command given after `--` (the polycrystal and its options, one time
step) is run three times, each in a directory of its own:

- displacement alone (--enable-paraview-export): collection
  OutputFile-Uniaxial-polycrystal, field u;
- von Mises stress alone (--disable-paraview-export
  --enable-export-von-Mises-stress): collection PolycrystalParaviewOutput,
  field vonMisesStress;
- both, in single precision (--enable-paraview-export
  --enable-export-von-Mises-stress --enable-paraview-single-precision):
  collection PolycrystalParaviewResults, fields u and vonMisesStress in the
  same files, with one mesh, all the reals in Float32.

The cycle 0 of each collection is the end of the time step.

In each case, no other collection is written.
usage: check-paraview-exports.py <directory> -- <command...>
"""
import os
import re
import subprocess
import sys

COLLECTIONS = ("OutputFile-Uniaxial-polycrystal", "PolycrystalParaviewOutput",
               "PolycrystalParaviewResults")

CASES = (
    ("displacement", ["--enable-paraview-export"],
     "OutputFile-Uniaxial-polycrystal", (0,), ("u",), "Float64"),
    ("von-Mises", ["--disable-paraview-export",
                   "--enable-export-von-Mises-stress"],
     "PolycrystalParaviewOutput", (0,), ("vonMisesStress",), "Float64"),
    ("both-single-precision", ["--enable-paraview-export",
                               "--enable-export-von-Mises-stress",
                               "--enable-paraview-single-precision"],
     "PolycrystalParaviewResults", (0,), ("u", "vonMisesStress"),
     "Float32"),
)


def check_case(root, command, case):
    """runs a case and returns the list of the errors found"""
    name, options, collection, cycles, fields, real = case
    directory = os.path.join(root, name)
    os.makedirs(directory, exist_ok=True)
    run = subprocess.run(command + options, cwd=directory,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True, check=False)
    if run.returncode != 0:
        return [f"{name}: the command failed ({run.returncode}):\n"
                + run.stdout[-2000:]]
    errors = []
    for other in COLLECTIONS:
        if (other != collection) and os.path.exists(
                os.path.join(directory, other)):
            errors.append(f"{name}: collection {other} written")
    for c in cycles:
        f = os.path.join(directory, collection, f"Cycle{c:06d}",
                         "proc000000.vtu")
        if not os.path.exists(f):
            errors.append(f"{name}: {f} not written")
            continue
        with open(f, "rb") as vtu:
            text = vtu.read().decode("latin-1")
        names = re.findall(r'<DataArray type="Float(?:32|64)" Name="([^"]+)"',
                           text)
        types = set(re.findall(r'type="(Float(?:32|64))"', text))
        if text.count("<Points>") != 1:
            errors.append(f"{name}, cycle {c}: {text.count('<Points>')} meshes")
        if sorted(names) != sorted(fields):
            errors.append(f"{name}, cycle {c}: fields {names}, "
                          f"{list(fields)} expected")
        if types != {real}:
            errors.append(f"{name}, cycle {c}: reals {sorted(types)}, "
                          f"{real} expected")
    return errors


def main():
    if (len(sys.argv) < 4) or (sys.argv[2] != "--"):
        sys.exit(__doc__)
    errors = []
    for case in CASES:
        errors += check_case(sys.argv[1], sys.argv[3:], case)
    for e in errors:
        print(e)
    print("success" if not errors else "failure")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
