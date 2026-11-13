import sys, numpy as np, subprocess, tempfile, os
def vols(path):
    w = tempfile.mkdtemp()
    subprocess.run(["awk", "-v", f"w={w}", r'''/^\$Nodes/ {getline n; for (i=0;i<n;i++) {getline; print > (w "/n.txt")}; next}
/^\$Elements/ {getline n; for (i=0;i<n;i++) {getline; if ($2==4) print $4, $6, $7, $8, $9 > (w "/t.txt")}; next}''', path], check=True)
    nd = np.fromfile(w + "/n.txt", sep=" ").reshape(-1, 4); t = np.fromfile(w + "/t.txt", sep=" ", dtype=np.int64).reshape(-1, 5)
    idx = np.full(int(nd[:, 0].max()) + 1, -1); idx[nd[:, 0].astype(int)] = np.arange(len(nd)); x = nd[:, 1:]
    c = idx[t[:, 1:]]; p = [x[c[:, i]] for i in range(4)]
    v = np.abs(np.einsum("ij,ij->i", p[1]-p[0], np.cross(p[2]-p[0], p[3]-p[0]))) / 6
    return np.bincount(t[:, 0], weights=v)
a, b = vols(sys.argv[1]), vols(sys.argv[2])
nz = a > 0
r = np.abs(b[nz] - a[nz]) / a[nz]
print(f"grains {nz.sum()} -> {(b > 0).sum()}, relative change of the grain volumes: max {r.max():.2e}, mean {r.mean():.2e}")
