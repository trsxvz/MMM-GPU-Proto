# Scaling of the beam, first pass (11/10/2026)

First pass of `cases/measure-beam-scaling.sh`, with the default memory
pool of hypre (4 GiB) and without the matrix free variant. The
computations on 8 cores (and their rerun, files `*-rerun.txt`) are 13 to
43 % slower than those of the next pass, for a cause not established
(compilations of another session were seen during these computations, but
the rerun, at a load of 1.1, gave the same times). Reference results and
analysis: `../scaling-beam-2026-10-11-0247`.
