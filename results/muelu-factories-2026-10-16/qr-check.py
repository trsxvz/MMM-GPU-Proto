# replicates LocalQRDecompFunctor of MueLu 17.2.1 (TentativePFactory_kokkos)
# and checks the tentative prolongator block Q[:, :n] against B = Q R
import numpy as np

def muelu_kokkos_qr(B):
    m, n = B.shape
    r = B.copy()
    qt = np.eye(m)
    for k in range(n):
        s = np.sum(r[k + 1:, k] ** 2)
        norm = np.sqrt(r[k, k] ** 2 + s)
        r[k, k] -= norm
        norm_x = np.sqrt(r[k, k] ** 2 + s)
        if norm_x == 0:
            r[k, k] = norm
            continue
        r[k:, k] /= norm_x
        for j in range(k + 1, n):
            si = np.dot(r[k:, k], r[k:, j])
            r[k:, j] -= 2 * si * r[k:, k]
        for j in range(k, m):
            si = np.dot(r[k:, k], qt[k:, j])
            qt[k:, j] -= 2 * si * r[k:, k]
        r[k, k] = norm
        r[k + 1:, k] = 0
    R = np.triu(r[:n, :n])
    return qt, R  # the functor uses q(j, k) = qt(j, k) as Q(j, k)

def rigid_body_modes(x):
    # 3D, unknowns ordered by node: (u_x, u_y, u_z) per node
    nn = x.shape[0]
    B = np.zeros((3 * nn, 6))
    c = x - x.mean(axis=0)
    for i in range(nn):
        for d in range(3):
            B[3 * i + d, d] = 1
        X, Y, Z = c[i]
        B[3 * i + 0, 3], B[3 * i + 1, 3] = -Y, X
        B[3 * i + 1, 4], B[3 * i + 2, 4] = -Z, Y
        B[3 * i + 0, 5], B[3 * i + 2, 5] = Z, -X
    return B

rng = np.random.default_rng(0)
for nodes in (2, 3, 5, 9, 27):
    x = rng.random((nodes, 3))
    B = rigid_body_modes(x)
    q, R = muelu_kokkos_qr(B)
    n = B.shape[1]
    P = q[:, :n]
    err = np.linalg.norm(P @ R - B) / np.linalg.norm(B)
    Pt = q.T[:, :n]
    errt = np.linalg.norm(Pt @ R - B) / np.linalg.norm(B)
    orth = np.linalg.norm(P.T @ P - np.eye(n))
    print(f"aggregate of {nodes:2d} nodes: |Q R - B|/|B| = {err:.2e} "
          f"(with the transpose of q: {errt:.2e}), |Q^T Q - I| = {orth:.2e}")

# controls: the same reflections, Q^T accumulated on all its columns then
# transposed (correct), and the reflections stored, Q accumulated backward
# on the trailing blocks (correct, as LAPACK's dorgqr)
def householder(B, mode):
    m, n = B.shape
    r = B.copy()
    qt = np.eye(m)
    vs = []
    for k in range(n):
        s = np.sum(r[k + 1:, k] ** 2)
        norm = np.sqrt(r[k, k] ** 2 + s)
        r[k, k] -= norm
        norm_x = np.sqrt(r[k, k] ** 2 + s)
        if norm_x == 0:
            r[k, k] = norm
            vs.append(None)
            continue
        v = r[k:, k] / norm_x
        vs.append(v)
        for j in range(k + 1, n):
            r[k:, j] -= 2 * np.dot(v, r[k:, j]) * v
        if mode == "full":
            for j in range(m):
                qt[k:, j] -= 2 * np.dot(v, qt[k:, j]) * v
        r[k, k] = norm
        r[k + 1:, k] = 0
    R = np.triu(r[:n, :n])
    if mode == "full":
        return qt.T, R
    q = np.eye(m)
    for k in reversed(range(n)):
        v = vs[k]
        if v is None:
            continue
        for j in range(k, m):
            q[k:, j] -= 2 * np.dot(v, q[k:, j]) * v
    return q, R

for mode in ("full", "backward"):
    for nodes in (2, 9, 27):
        x = rng.random((nodes, 3))
        B = rigid_body_modes(x)
        q, R = householder(B, mode)
        P = q[:, :6]
        print(f"{mode}: {nodes:2d} nodes: |Q R - B|/|B| = "
              f"{np.linalg.norm(P @ R - B) / np.linalg.norm(B):.2e}, "
              f"|Q^T Q - I| = {np.linalg.norm(P.T @ P - np.eye(6)):.2e}")
