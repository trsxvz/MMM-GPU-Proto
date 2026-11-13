# Variants of the behaviour of the polycrystal

Variants of `mm-opera-hpc/polycrystal/law/Mono_UO2_CosH_Jaco3.mfront` (same
behaviour name, same equations and tolerances) tested on 14/11/2026 for the
robustness of the polycrystal series (see
`~/Bureau/MMM-GPU-Proto-reflexion-robustesse-divergence-2026-11-14.md` and
`JOURNAL.md`):

- **A** (`Mono_UO2_CosH_Jaco3-A.mfront`): predictor of the volumetric part
  of the elastic strain. The plastic flow being isochoric, the volumetric
  part of the elastic strain at the end of the step is the one of the
  elastic prediction; the bound of 1e-5 per iteration and per component of
  the elastic strain made the local Newton method reach its 1000 iterations
  as soon as the volumetric increment reached about 1 %.
- **B** (`Mono_UO2_CosH_Jaco3-B.mfront`): A, and the bound of 1e-5 per
  component replaced by a bound of the increment of the resolved shear
  stresses (tau0 of each system per iteration), the direction of the
  Newton correction being kept.

At one point (MTest, standard TFEL 5.1), the stresses of A and B differ
from the ones of the law by 2e-7 at most (local tolerance 1e-11), and the
failures of the integration disappear up to increments of 2 % (A) or 10 %
(B). They are built for the GPU (generic-parallel, hip) and the CPU
(generic) by `adastra-fm/build-laws.slurm`.
