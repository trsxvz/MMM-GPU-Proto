# Writes the matrix of the 7-point Laplacian on a n x n x n grid in the
# Matrix Market format, the input of the amgx_capi example of AmgX.
# usage: python3 poisson.py n file
import sys

n = int(sys.argv[1])


def index(i, j, k):
    return (i * n + j) * n + k


entries = []
for i in range(n):
    for j in range(n):
        for k in range(n):
            r = index(i, j, k)
            entries.append((r, r, 6.0))
            for di, dj, dk in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                               (0, 0, 1), (0, 0, -1)):
                a, b, c = i + di, j + dj, k + dk
                if 0 <= a < n and 0 <= b < n and 0 <= c < n:
                    entries.append((r, index(a, b, c), -1.0))
with open(sys.argv[2], 'w') as f:
    f.write('%%MatrixMarket matrix coordinate real general\n')
    f.write(f'{n**3} {n**3} {len(entries)}\n')
    for r, c, v in entries:
        f.write(f'{r + 1} {c + 1} {v}\n')
