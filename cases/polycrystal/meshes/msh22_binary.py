#!/usr/bin/env python3
"""Gmsh 2.2 mesh in ASCII to the same mesh in binary: usage:
msh22_binary.py in.msh out.msh

The nodes and the elements keep their numbers, their order and their
values (the coordinates read in the ASCII file are written as the same
doubles), the other sections ($PhysicalNames, $Periodic, which stay in
ASCII in the binary format 2.2) are copied. MFEM reads both formats
(mesh/gmsh.cpp of MFEM 4.10): the binary file gives the same mesh,
read much faster (the reading of the ASCII mesh of 6,400 grains, 570 MB,
took about 50 s on each process)."""
import struct
import sys


def main():
    src, dst = sys.argv[1], sys.argv[2]
    with open(src) as f, open(dst, "wb") as out:
        line = f.readline()
        while line:
            tag = line.strip()
            if tag == "$MeshFormat":
                version, file_type, size = f.readline().split()
                if version != "2.2" or file_type != "0" or size != "8":
                    sys.exit(f"{src}: not an ASCII Gmsh 2.2 mesh")
                f.readline()  # $EndMeshFormat
                out.write(b"$MeshFormat\n2.2 1 8\n")
                out.write(struct.pack("<i", 1))
                out.write(b"\n$EndMeshFormat\n")
            elif tag == "$Nodes":
                n = int(f.readline())
                out.write(b"$Nodes\n%d\n" % n)
                pack = struct.Struct("<iddd").pack
                for _ in range(n):
                    w = f.readline().split()
                    out.write(pack(int(w[0]), float(w[1]), float(w[2]),
                                   float(w[3])))
                f.readline()  # $EndNodes
                out.write(b"\n$EndNodes\n")
            elif tag == "$Elements":
                n = int(f.readline())
                out.write(b"$Elements\n%d\n" % n)
                # blocks of consecutive elements of the same type and number
                # of tags, in the order of the file
                block, key = [], None

                def flush():
                    if block:
                        t, ntags = key
                        out.write(struct.pack("<iii", t, len(block), ntags))
                        out.write(struct.pack(
                            "<%di" % sum(len(e) for e in block),
                            *[v for e in block for v in e]))
                        block.clear()

                for _ in range(n):
                    w = [int(v) for v in f.readline().split()]
                    number, t, ntags = w[0], w[1], w[2]
                    if key != (t, ntags):
                        flush()
                        key = (t, ntags)
                    block.append([number] + w[3:])
                    if len(block) == 65536:
                        flush()
                flush()
                f.readline()  # $EndElements
                out.write(b"\n$EndElements\n")
            elif tag.startswith("$"):
                # sections kept in ASCII ($PhysicalNames, $Periodic, ...)
                out.write(line.encode())
                end = "$End" + tag[1:]
                for line in f:
                    out.write(line.encode())
                    if line.strip() == end:
                        break
            line = f.readline()


if __name__ == "__main__":
    main()
