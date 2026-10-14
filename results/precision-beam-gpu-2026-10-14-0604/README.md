# Précision float contre double sur le GPU, poutre (14/10/2026)

Poutre r4 du prototype (111 843 inconnues), 10 pas, Newton (critère 1e-4,
et celui du cas, 1e-8, en double), PCG de hypre préconditionné par
BoomerAMG (tolérance 1e-4), tout sur le GPU (`--execution gpu`, loi du
backend `cuda`), un processus, GPU NVIDIA RTX 5080 libre (la session
voisine en pause). Script : `PERIPHERIQUE=gpu cases/measure-precision-beam.sh`.
Double : build/cas-9b (mfem-mgis gpu-proto installé à part, environnement
mmm-gpu). Float : build/cas-float-gpu (environnement mmm-float-gpu, MFEM et
hypre en simple avec CUDA). Même version de mfem-mgis.

| Calcul | Critère | Temps min. (s) | Médiane (s) | Itérations | Norme des déplacements |
|---|---|---|---|---|---|
| double | 1e-8 | 2,91 | 3,02 | 32 | 48,83704586 |
| double | 1e-4 | 1,86 | 1,86 | 16 | 48,83715008 |
| float | 1e-4 | 1,42 | 1,42 | 16 | 48,85159302 et 48,84963608 (2 calculs sur 3) |

Temps par phase (premier passage à 1e-4) : résolution linéaire 1,13 s en
double, 0,81 s en float ; prédiction 0,55 et 0,49 s ; le reste est petit.

- Le float est 1,31 fois plus rapide quand il converge, le gain venant
  surtout de la résolution linéaire.
- Mais le chemin float sur le GPU n'est pas fiable : un calcul sur trois a
  divergé (Newton au dixième pas), les calculs convergés diffèrent entre
  eux (4e-5) et du double (3e-4, contre 8e-5 pour le float sur le CPU).
  Analyse (`scratch`, 14/10) : 4 échecs sur 8 avec l'assemblage atomique,
  8 sur 8 par noyaux avec l'assemblage ordonné et un solveur itératif (un
  calcul avec MUMPS converge) ; le calcul converge si les éléments
  finis sont sur le CPU, les lois et le solveur restant sur le GPU, et
  diverge si les éléments finis sont sur le GPU, quels que soient les
  placements des lois et du solveur (PCG et BoomerAMG). Les mêmes noyaux
  en float sur le CPU convergent. La convergence de Newton se dégrade de
  pas en pas sur le GPU.
  En double, au critère du cas (1e-8), les trois calculs sont identiques
  au bit près ; au critère 1e-4, ils diffèrent de 4e-7 (assemblage
  atomique).
- Diagnostic (results/diagnostic-float-gpu-2026-10-14) : les contractions
  en multiplications-additions fusionnées de nvcc dans les noyaux de raideur
  faisaient diverger le Newton en float. Corrigé dans mfem-mgis (en float,
  noyaux de raideur compilés pour le GPU sans contraction). Mesure refaite :
  results/precision-beam-gpu-2026-10-14-0843 (float 1,37 fois plus
  rapide que le double, 3 calculs sur 3 convergés, à 3,4e-5 du double).
