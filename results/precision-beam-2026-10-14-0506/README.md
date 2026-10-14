# Float against double precision, first run (14/10/2026)

Same setup as `results/precision-beam-2026-10-14-0514`, with the Newton
criterion of the case (1e-8) and 1.2e-5 for both precisions. The six
computations in float fail: the relative tolerance, bounded to 100 ε =
1.2e-5 by mfem-mgis in float, is below the floor of the relative residual
at r4 (about 1.3e-5 from the sixth step on, see
`float-case-criterion-1.txt`), and Newton stops after 50 iterations. The
computations in double are the reference of the second run.
