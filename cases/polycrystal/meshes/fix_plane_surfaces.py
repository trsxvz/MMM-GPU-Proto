#!/usr/bin/env python3
"""Plane surfaces of a geometry of Merope that are not planar, defined as
surfaces filling their curve loop: usage:
fix_plane_surfaces.py in.geo out.geo [tolerance]
fix_plane_surfaces.py --log gmsh.log in.geo out.geo

When the Laguerre tessellation has a tiny face seen by one cell only, Merope
merges its points into one ("An entire face should be merged", without
checks): the faces around this point are no longer planar (1e-4 for edges
of 1 to 8 in a geometry of 12,800 grains), OpenCASCADE refuses them ("Could
not add plane surface") and the volumes they bound are left without
elements. These plane surfaces are rewritten with `Surface` instead of
`Plane Surface`: OpenCASCADE then builds a surface through the loop, as flat
as the loop. The other lines are copied unchanged.

By default, the plane surfaces whose vertices are farther than `tolerance`
(default 1e-8, absolute, below the tolerance of OpenCASCADE, 1e-7) from
their least squares plane are rewritten, which needs no run of gmsh. With
--log, the lines given by the errors of a parsing of in.geo by gmsh
(option -parse_and_exit) are rewritten instead."""
import re
import sys

import numpy as np


def from_log(log):
    bad = set()
    with open(log, errors="replace") as f:
        for line in f:
            m = re.search(r"line (\d+): Could not add plane surface", line)
            if m:
                bad.add(int(m.group(1)))
    return bad


def non_planar(src, tolerance):
    """line numbers of the plane surfaces whose loops are not planar, and
    the largest distance to the plane"""
    points, lines, loops, surfaces = {}, {}, {}, []
    number = re.compile(r"[-+0-9.eE]+")
    with open(src) as f:
        for n, line in enumerate(f, start=1):
            if "= newp;" in line:
                name = line.split("=", 1)[0].strip()
                values = line.split("{", 1)[1].split("}", 1)[0].split(",")
                points[name] = [float(v) for v in values[:3]]
            elif "= newc;" in line and "Line(" in line:
                name = line.split("=", 1)[0].strip()
                a, b = line.split("{", 1)[1].split("}", 1)[0].split(",")
                lines[name] = (a.strip(), b.strip())
            elif "= newcl;" in line:
                name = line.split("=", 1)[0].strip()
                loops[name] = [c.strip() for c in
                               line.split("{", 1)[1].split("}", 1)[0].split(",")]
            elif "Plane Surface(" in line:
                loop = line.split("{", 1)[1].split("}", 1)[0].strip()
                surfaces.append((n, loop))
    bad, largest = set(), 0.0
    for n, loop in surfaces:
        ids = []
        for c in loops[loop]:
            a, b = lines[c.lstrip("-")]
            ids.append(b if c.startswith("-") else a)
        xyz = np.array([points[i] for i in ids])
        centered = xyz - xyz.mean(axis=0)
        # distance to the least squares plane: smallest singular value
        normal = np.linalg.svd(centered)[2][-1]
        distance = np.abs(centered @ normal).max()
        largest = max(largest, distance)
        if distance > tolerance:
            bad.add(n)
    return bad, largest, len(surfaces)


def main():
    args = sys.argv[1:]
    if args and args[0] == "--log":
        bad = from_log(args[1])
        src, dst = args[2], args[3]
        how = f"errors of {args[1]}"
    else:
        src, dst = args[0], args[1]
        tolerance = float(args[2]) if len(args) > 2 else 1e-8
        bad, largest, count = non_planar(src, tolerance)
        how = (f"{count} plane surfaces, distance to their plane above "
               f"{tolerance:g} (largest {largest:.3g})")
    fixed = 0
    with open(src) as f, open(dst, "w") as out:
        for n, line in enumerate(f, start=1):
            if n in bad:
                if "Plane Surface(" not in line:
                    sys.exit(f"{src}:{n}: not a plane surface: {line.strip()}")
                line = line.replace("Plane Surface(", "Surface(", 1)
                fixed += 1
            out.write(line)
    print(f"{fixed} plane surfaces defined as surfaces ({how})")
    for n in sorted(bad):
        print(f"  line {n}")


if __name__ == "__main__":
    main()
