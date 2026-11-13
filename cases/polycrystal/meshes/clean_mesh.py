#!/usr/bin/env python3
"""Removes the short edges of a periodic polycrystal mesh (gmsh msh 2.2,
tetrahedra and triangles, $Periodic section), which give flat tetrahedra:
the short edges of the Laguerre tessellation (their length decreases with
the number of grains) are meshed by edges of the same length, between
tetrahedra of the size of the mesh.

usage: clean_mesh.py in.msh eps out.msh

The edges shorter than eps are collapsed (union of their nodes, propagated to
the periodic images). Each merged class takes the position of its most
constrained node: on the box (corner, then edge, then face), then on the
grains (quadruple point, triple line, grain boundary), so that the box and
the boundaries of the grains move as little as possible (the collapsed
short edges of the tessellation are removed, as in the regularization of
Neper). The periodic nodes are then put exactly at the position of their
master plus the translation. The degenerate tetrahedra and triangles are
removed. The script stops (status 1) if a tetrahedron is reversed or loses
more than 99 % of its volume. Prints the quality before and after.
"""
import os
import subprocess
import sys
import tempfile

import numpy as np


def quality(xyz, tets, label):
    p = [xyz[tets[:, i]] for i in range(4)]
    sv = np.einsum("ij,ij->i", p[1] - p[0], np.cross(p[2] - p[0], p[3] - p[0])) / 6
    vol = np.abs(sv)

    def area(a, b, d):
        return 0.5 * np.linalg.norm(np.cross(b - a, d - a), axis=1)

    A = np.maximum.reduce([area(p[0], p[1], p[2]), area(p[0], p[1], p[3]),
                           area(p[0], p[2], p[3]), area(p[1], p[2], p[3])])
    h = 3 * vol / np.where(A > 0, A, 1)
    e = np.array([np.linalg.norm(p[i] - p[j], axis=1)
                  for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))])
    asp = e.max(axis=0) / np.where(h > 0, h, 1e-300)
    print(f"{label}: {len(tets)} tets, volume {vol.sum():.6g}, min edge "
          f"{e.min():.3g}, min height {h.min():.3g}, aspect max {asp.max():.4g}, "
          f">30: {(asp > 30).sum()}, >100: {(asp > 100).sum()}, "
          f">1000: {(asp > 1000).sum()}")
    return sv, asp


def split(path, work):
    """Nodes, tetrahedra, triangles and periodic section in text files."""
    script = r'''
/^\$Nodes/ {getline n; for (i = 0; i < n; i++) {getline; print > (w "/nodes.txt")}; next}
/^\$Elements/ {getline n; for (i = 0; i < n; i++) {getline;
  if ($2 == 4 && $3 == 2) print $1, $4, $5, $6, $7, $8, $9 > (w "/tets.txt");
  else if ($2 == 2 && $3 == 2) print $1, $4, $5, $6, $7, $8 > (w "/tris.txt");
  else if ($2 == 4 || $2 == 2) {print "unexpected tags: " $0 > "/dev/stderr"; exit 1}};
  next}
/^\$Periodic/ {p = 1}
p {print > (w "/periodic.txt")}
/^\$EndPeriodic/ {p = 0}
/^\$MeshFormat/ {getline; print > (w "/format.txt")}
'''
    subprocess.run(["awk", "-v", f"w={work}", script, path], check=True)


def main():
    path, eps, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    work = tempfile.mkdtemp(prefix="clean-mesh-")
    split(path, work)
    fmt = open(os.path.join(work, "format.txt")).read().split()
    if fmt[0] != "2.2" or fmt[1] != "0":
        sys.exit(f"msh 2.2 ASCII expected, got {fmt}")
    nodes = np.fromfile(os.path.join(work, "nodes.txt"), sep=" ").reshape(-1, 4)
    ids = nodes[:, 0].astype(np.int64)
    xyz = nodes[:, 1:].copy()
    n = len(ids)
    idx = np.full(ids.max() + 1, -1, dtype=np.int64)
    idx[ids] = np.arange(n)
    T = np.fromfile(os.path.join(work, "tets.txt"), sep=" ",
                    dtype=np.int64).reshape(-1, 7)
    F = np.fromfile(os.path.join(work, "tris.txt"), sep=" ",
                    dtype=np.int64).reshape(-1, 6)
    tets = idx[T[:, 3:]]
    tris = idx[F[:, 3:]]
    grain = T[:, 1]
    # periodic section: entities (header, affine line, pairs)
    plines = open(os.path.join(work, "periodic.txt")).read().split("\n")
    assert plines[0] == "$Periodic"
    entities = []
    q = 2
    while q < len(plines) and plines[q] != "$EndPeriodic":
        header = plines[q]
        q += 1
        affine = None
        if plines[q].startswith("Affine"):
            affine = plines[q]
            q += 1
        m = int(plines[q])
        pairs = np.array([l.split() for l in plines[q + 1:q + 1 + m]],
                         dtype=np.int64).reshape(-1, 2)
        q += 1 + m
        entities.append((header, affine, pairs))
    print(f"{n} nodes, {len(tets)} tets, {len(tris)} triangles, "
          f"{len(entities)} periodic entities, "
          f"{sum(len(p) for _, _, p in entities)} periodic pairs")
    sv0, asp0 = quality(xyz, tets, "before")
    # period: the RVE of Merope is a union of whole grains (its periodic
    # boundary follows the grains), of the volume of the cube of side L
    period = np.abs(sv0).sum() ** (1 / 3)
    pair_list = np.vstack([idx[p] for _, _, p in entities])  # (slave, master)
    offsets = np.round((xyz[pair_list[:, 0]] - xyz[pair_list[:, 1]]) /
                       period).astype(np.int64)
    residual = np.abs(xyz[pair_list[:, 0]] - xyz[pair_list[:, 1]] -
                      offsets * period).max()
    print(f"period {period:.10g}, translation residual of the pairs {residual:.3g}")
    # constraint of each node: periodic boundary (number of distinct
    # translations to its images: 1 on a face, more on edges and corners),
    # then grains (number of grains around the node)
    images_count = {}
    for (s_, m_), off in zip(pair_list.tolist(), offsets.tolist()):
        images_count.setdefault(s_, set()).add(tuple(off))
        images_count.setdefault(m_, set()).add(tuple(-o for o in off))
    on_box = np.zeros(n, dtype=np.int64)
    for a, offs in images_count.items():
        on_box[a] = len(offs)
    key = np.unique(np.stack([tets.ravel(), np.repeat(grain, 4)], axis=1), axis=0)
    ngrains = np.bincount(key[:, 0], minlength=n)
    priority = on_box * 16 + np.minimum(ngrains, 15)
    # short edges
    edges = np.vstack([np.sort(tets[:, [a, b]], axis=1)
                       for a, b in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))])
    edges = np.unique(edges, axis=0)
    length = np.linalg.norm(xyz[edges[:, 0]] - xyz[edges[:, 1]], axis=1)
    short = edges[length < eps]
    print(f"edges shorter than {eps}: {len(short)} of {len(edges)} "
          f"(shortest {length.min():.3g})")
    images = {}
    for (s_, m_), off in zip(pair_list.tolist(), offsets.tolist()):
        images.setdefault(m_, []).append((s_, tuple(off)))
        images.setdefault(s_, []).append((m_, tuple(-o for o in off)))

    def collapse(short, verbose):
        """Classes of the nodes merged by the given short edges (with their
        periodic images), new positions, kept elements."""
        parent = np.arange(n)

        def find(a):
            root = a
            while parent[root] != root:
                root = parent[root]
            while parent[a] != root:
                parent[a], a = root, parent[a]
            return root

        def union(a, b):
            a, b = find(a), find(b)
            if a == b:
                return False
            parent[max(a, b)] = min(a, b)
            return True

        for a, b in short:
            union(a, b)
        # periodic propagation: the images by the same translation of merged
        # nodes are merged
        passes = 0
        while True:
            passes += 1
            changed = False
            merged = set(short.ravel().tolist()) | {
                a for a in np.where(parent != np.arange(n))[0].tolist()}
            classes = {}
            for a in merged:
                classes.setdefault(find(a), set()).add(a)
            for r, members in classes.items():
                by_offset = {}
                for a in members:
                    for b, off in images.get(a, ()):
                        by_offset.setdefault(off, []).append(b)
                for bs in by_offset.values():
                    for b in bs[1:]:
                        changed |= union(bs[0], b)
            if not changed:
                break
        rep = np.array([find(a) for a in range(n)])
        moved = np.where(rep != np.arange(n))[0]
        if verbose:
            print(f"nodes merged: {len(moved)} (periodic passes: {passes})")
        # position of each class: its most constrained node
        new_xyz = xyz.copy()
        for r in np.unique(rep[moved]):
            members = np.where(rep == r)[0]
            best = members[np.argmax(priority[members])]
            new_xyz[r] = xyz[best]
        # exact periodicity: slave = master + translation (chains repeated)
        reps = rep[pair_list]
        for _ in range(4):
            target = new_xyz[reps[:, 1]] + offsets * period
            if np.abs(new_xyz[reps[:, 0]] - target).max() == 0:
                break
            new_xyz[reps[:, 0]] = target
        err = np.abs(new_xyz[reps[:, 0]] - new_xyz[reps[:, 1]] -
                     offsets * period).max()
        new_tets = rep[tets]
        s_ = np.sort(new_tets, axis=1)
        keep_t = np.all(s_[:, 1:] != s_[:, :-1], axis=1)
        new_tris = rep[tris]
        s_ = np.sort(new_tris, axis=1)
        keep_f = np.all(s_[:, 1:] != s_[:, :-1], axis=1)
        return rep, new_xyz, new_tets, new_tris, keep_t, keep_f, err

    # collapses; the short edges around elements reversed, crushed or made
    # much worse are given up (with their periodic images), until none
    forbidden = np.zeros(n, dtype=bool)
    for attempt in range(1, 21):
        use = short[~(forbidden[short[:, 0]] | forbidden[short[:, 1]])]
        rep, new_xyz, new_tets, new_tris, keep_t, keep_f, err = collapse(
            use, attempt == 1)
        p = [new_xyz[new_tets[:, i]] for i in range(4)]
        sv1 = np.einsum("ij,ij->i", p[1] - p[0],
                        np.cross(p[2] - p[0], p[3] - p[0])) / 6
        bad = keep_t & ((np.sign(sv1) != np.sign(sv0)) |
                        (np.abs(sv1) < 1e-2 * np.abs(sv0)))
        # aspect of the kept tetrahedra
        q = [new_xyz[new_tets[keep_t, i]] for i in range(4)]
        vol = np.abs(sv1[keep_t])

        def area(a, b, d):
            return 0.5 * np.linalg.norm(np.cross(b - a, d - a), axis=1)

        A = np.maximum.reduce([area(q[0], q[1], q[2]), area(q[0], q[1], q[3]),
                               area(q[0], q[2], q[3]), area(q[1], q[2], q[3])])
        h = 3 * vol / np.where(A > 0, A, 1)
        e = np.array([np.linalg.norm(q[i] - q[j], axis=1)
                      for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3),
                                   (2, 3))]).max(axis=0)
        asp = e / np.where(h > 0, h, 1e-300)
        worse = np.zeros(len(tets), dtype=bool)
        worse[np.where(keep_t)[0]] = (asp > 200) & (asp > 2 * asp0[keep_t])
        bad |= worse
        print(f"attempt {attempt}: {len(use)} short edges collapsed, "
              f"{bad.sum()} elements reversed, crushed or made much worse")
        if not bad.any():
            break
        # nodes of the classes of the bad elements, and their images
        bad_reps = np.unique(new_tets[bad])
        nodes_bad = np.isin(rep, bad_reps)
        for a in np.where(nodes_bad)[0].tolist():
            for b, _ in images.get(a, ()):
                nodes_bad[b] = True
        forbidden |= nodes_bad
    else:
        sys.exit("invalid mesh after 20 attempts: no output written")
    print(f"short edges kept: {len(short) - len(use)}")
    print(f"periodicity error after the update: {err:.3g}")
    print(f"tets removed: {(~keep_t).sum()}, triangles removed: {(~keep_f).sum()}")
    quality(new_xyz, new_tets[keep_t], "after ")
    # output
    # nodes renumbered contiguously from 1
    kept = np.where(rep == np.arange(n))[0]
    newid = np.zeros(n, dtype=np.int64)
    newid[kept] = np.arange(1, len(kept) + 1)
    with open(out, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n$Nodes\n")
        f.write(f"{len(kept)}\n")
        np.savetxt(f, np.column_stack([newid[kept], new_xyz[kept]]),
                   fmt=["%d", "%.17g", "%.17g", "%.17g"])
        f.write("$EndNodes\n$Elements\n")
        tri_out = np.column_stack([F[keep_f, 1:3], newid[new_tris[keep_f]]])
        tet_out = np.column_stack([T[keep_t, 1:3], newid[new_tets[keep_t]]])
        nf, nt = len(tri_out), len(tet_out)
        f.write(f"{nf + nt}\n")
        np.savetxt(f, np.column_stack([np.arange(1, nf + 1), np.full(nf, 2),
                                       np.full(nf, 2), tri_out]), fmt="%d")
        np.savetxt(f, np.column_stack([np.arange(nf + 1, nf + nt + 1),
                                       np.full(nt, 4), np.full(nt, 2), tet_out]),
                   fmt="%d")
        f.write("$EndElements\n$Periodic\n")
        f.write(f"{len(entities)}\n")
        for header, affine, pairs in entities:
            p = rep[idx[pairs]]
            p = np.unique(p[p[:, 0] != p[:, 1]], axis=0)
            f.write(header + "\n")
            if affine:
                f.write(affine + "\n")
            f.write(f"{len(p)}\n")
            np.savetxt(f, newid[p], fmt="%d")
        f.write("$EndPeriodic\n")
    print(f"written: {out}")


if __name__ == "__main__":
    main()
