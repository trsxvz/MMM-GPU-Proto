#!/usr/bin/env python3
"""Compares the ParaView files written by two computations of the polycrystal.

The files are the ones written by MFEM (`.vtu` of each process, binary,
compressed by zlib). Each data array is decoded and compared value by
value:

- `split <before> <after>`: the displacement and the von Mises stress
  written in two collections (`OutputFile-Uniaxial-polycrystal` and
  `PolycrystalParaviewOutput`) against the same fields written in one
  collection (`PolycrystalParaviewResults`): at each cycle (the end of a
  time step), the mesh (coordinates of the nodes, connectivity, offsets,
  types, attributes), the displacement and the von Mises stress must be the
  same bit for bit;
- `same <collection a> <collection b>`: two collections, cycle by cycle,
  bit for bit;
- `float <double collection> <single collection>`: each real of the second
  collection must be the rounding to the nearest single precision value of
  the one of the first collection, the other arrays being the same.

usage: compare-vtu.py split|same|float <directory a> <directory b> [tolerance]
(directories of the computations for `split`, of the collections
otherwise). With a tolerance (`split` and `same`), the reals are compared
with this tolerance relative to the largest absolute value of their array,
the largest relative differences being printed: for computations which
are not reproducible bit for bit. The exit status is 1 if a difference is
found.
"""
import base64
import glob
import math
import os
import re
import sys
import zlib

import numpy

TYPES = {"Float64": "<f8", "Float32": "<f4", "Int32": "<i4", "Int64": "<i8",
         "UInt8": "u1", "UInt32": "<u4", "UInt64": "<u8"}


def decode(text, dtype):
    """values of a data array written in base64, compressed by zlib, with
    a header of 32 bits integers (number of blocks, size of the blocks,
    size of the last block, compressed sizes)"""
    text = b"".join(text.split())
    n = int(numpy.frombuffer(base64.b64decode(text[:8])[:4], "<u4")[0])
    length = 4 * math.ceil(4 * (3 + n) / 3)
    header = numpy.frombuffer(base64.b64decode(text[:length]), "<u4")
    sizes = header[3:3 + n]
    payload = base64.b64decode(text[length:])
    chunks, offset = [], 0
    for s in sizes:
        chunks.append(zlib.decompress(payload[offset:offset + int(s)]))
        offset += int(s)
    return numpy.frombuffer(b"".join(chunks), dtype)


def read(path):
    """data arrays of a .vtu file, by name (`Points` for the coordinates)"""
    with open(path, "rb") as f:
        data = f.read()
    if b'compressor="vtkZLibDataCompressor"' not in data[:400]:
        sys.exit(f"{path}: not compressed by zlib")
    arrays = {}
    for m in re.finditer(rb"<DataArray([^>]*)>(.*?)</DataArray>", data, re.S):
        attributes = dict(re.findall(rb'(\w+)="([^"]*)"', m.group(1)))
        name = attributes.get(b"Name", b"Points").decode()
        kind = attributes[b"type"].decode()
        arrays[name] = (kind, decode(m.group(2), TYPES[kind]))
    return arrays


def files(collection, cycle):
    """files of the processes for a cycle of a collection"""
    return sorted(glob.glob(os.path.join(collection, f"Cycle{cycle:06d}",
                                         "proc*.vtu")))


def cycles(collection):
    """cycles of a collection"""
    return sorted(int(os.path.basename(d)[5:])
                  for d in glob.glob(os.path.join(collection, "Cycle*")))


class Comparison:
    """counts of the compared values and of the differences"""

    def __init__(self, tolerance=None):
        self.arrays = 0
        self.values = 0
        self.differences = []
        self.tolerance = tolerance
        # largest relative difference of each array name
        self.largest = {}

    def same(self, label, a, b):
        """compares two arrays bit for bit (same type, same values), or
        within the tolerance for the reals"""
        self.arrays += 1
        self.values += a[1].size
        if (self.tolerance is not None) and a[0].startswith("Float") and \
                (a[0] == b[0]) and (a[1].shape == b[1].shape):
            x, y = a[1].astype(numpy.float64), b[1].astype(numpy.float64)
            scale = max(float(numpy.abs(x).max(initial=0)), 1e-300)
            r = float(numpy.abs(x - y).max(initial=0)) / scale
            name = label.split(", ")[-1]
            self.largest[name] = max(self.largest.get(name, 0), r)
            if r > self.tolerance:
                self.differences.append(f"{label}: relative difference "
                                        f"{r:.3e}")
            return
        if (a[0] != b[0]) or (a[1].shape != b[1].shape) or \
                (a[1].tobytes() != b[1].tobytes()):
            self.report(label, a, b)

    def rounded(self, label, d, s):
        """checks that an array of reals in single precision holds the
        roundings to the nearest of the ones of an array in double
        precision"""
        self.arrays += 1
        self.values += d[1].size
        if (d[0] != "Float64") or (s[0] != "Float32") or \
                (d[1].astype(numpy.float32).tobytes() != s[1].tobytes()):
            self.report(label, d, s)

    def report(self, label, a, b):
        """records a difference"""
        if a[1].shape != b[1].shape:
            self.differences.append(
                f"{label}: {a[0]} {a[1].shape} against {b[0]} {b[1].shape}")
            return
        d = numpy.abs(a[1].astype(numpy.float64) - b[1].astype(numpy.float64))
        self.differences.append(f"{label}: {a[0]} against {b[0]}, "
                                f"{int(numpy.count_nonzero(d))} values differ, "
                                f"largest difference {d.max():.6e}")

    def end(self):
        """prints the outcome and returns the exit status"""
        for name, r in sorted(self.largest.items()):
            print(f"{name}: largest relative difference {r:.3e}")
        for d in self.differences[:40]:
            print(d)
        print(f"{self.arrays} arrays and {self.values} values compared, "
              f"{len(self.differences)} differences")
        return 1 if self.differences else 0


def compare_collections(a, b, check, tolerance=None):
    """compares the arrays of two collections, cycle by cycle"""
    c = Comparison(tolerance)
    ca, cb = cycles(a), cycles(b)
    if (not ca) or (ca != cb):
        sys.exit(f"cycles {ca} of {a} and {cb} of {b}")
    for k in ca:
        fa, fb = files(a, k), files(b, k)
        if [os.path.basename(f) for f in fa] != \
                [os.path.basename(f) for f in fb]:
            sys.exit(f"cycle {k}: other files")
        for x, y in zip(fa, fb):
            ra, rb = read(x), read(y)
            if sorted(ra) != sorted(rb):
                c.differences.append(f"{y}: arrays {sorted(rb)}, "
                                     f"{sorted(ra)} expected")
                continue
            for name in ra:
                label = f"cycle {k}, {os.path.basename(x)}, {name}"
                if (check == "float") and ra[name][0].startswith("Float"):
                    c.rounded(label, ra[name], rb[name])
                else:
                    c.same(label, ra[name], rb[name])
    print(f"cycles {ca[0]} to {ca[-1]}")
    return c.end()


def compare_split(before, after, tolerance=None):
    """compares two collections of one field each to one collection with
    both fields"""
    u = os.path.join(before, "OutputFile-Uniaxial-polycrystal")
    vm = os.path.join(before, "PolycrystalParaviewOutput")
    both = os.path.join(after, "PolycrystalParaviewResults")
    c = Comparison(tolerance)
    cu, cvm, cboth = cycles(u), cycles(vm), cycles(both)
    if (not cboth) or (cu != cboth) or (cvm != cboth):
        sys.exit(f"cycles {cu} (displacement), {cvm} (von Mises stress), "
                 f"{cboth} (both)")
    mesh = ("Points", "connectivity", "offsets", "types", "attribute")
    for k in cboth:
        for fb in files(both, k):
            p = os.path.basename(fb)
            rb = read(fb)
            ru = read(os.path.join(u, f"Cycle{k:06d}", p))
            for name in mesh + ("u",):
                c.same(f"cycle {k}, {p}, {name}", ru[name], rb[name])
            rvm = read(os.path.join(vm, f"Cycle{k:06d}", p))
            for name in mesh + ("vonMisesStress",):
                c.same(f"cycle {k}, {p}, {name}", rvm[name], rb[name])
    print(f"cycles {cboth[0]} to {cboth[-1]}")
    return c.end()


def main():
    if (len(sys.argv) not in (4, 5)) or \
            (sys.argv[1] not in ("split", "same", "float")) or \
            ((len(sys.argv) == 5) and (sys.argv[1] == "float")):
        sys.exit(__doc__)
    tolerance = float(sys.argv[4]) if len(sys.argv) == 5 else None
    if sys.argv[1] == "split":
        sys.exit(compare_split(sys.argv[2], sys.argv[3], tolerance))
    sys.exit(compare_collections(sys.argv[2], sys.argv[3], sys.argv[1],
                                 tolerance))


if __name__ == "__main__":
    main()
