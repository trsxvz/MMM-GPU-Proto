# Précision float contre double sur le GPU, poutre, après correction (14/10/2026)

Même mesure que results/precision-beam-gpu-2026-10-14-0604, après la
correction de mfem-mgis : en float, les noyaux des matrices de raideur sont
compilés pour le GPU sans contraction en multiplications-additions
fusionnées ; en double, le code est inchangé (voir
results/diagnostic-float-gpu-2026-10-14).

Poutre r4 du prototype (111 843 inconnues), 10 pas, Newton (critère 1e-4, et
celui du cas, 1e-8, en double), PCG de hypre préconditionné par BoomerAMG
(tolérance 1e-4), tout sur le GPU (`--execution gpu`, assemblage atomique,
loi du backend `cuda`), un processus sur le cœur 0, GPU NVIDIA RTX 5080
libre. Aucun autre calcul pendant la mesure (la charge de 3,5 sur la
minute précédente vient des suites de tests finies juste avant). Script : `PERIPHERIQUE=gpu cases/measure-precision-beam.sh`.
Double : build/cas-9b (environnement mmm-gpu). Float : build/cas-float-gpu
(environnement mmm-float-gpu, MFEM et hypre en simple avec CUDA).

| Calcul | Critère | Temps min. (s) | Médiane (s) | Itérations | Norme des déplacements |
|---|---|---|---|---|---|
| double | 1e-8 | 2,89 | 2,94 | 31 à 32 | 48,83704586 |
| double | 1e-4 | 1,85 | 1,85 | 16 | 48,83715811 à 48,83716422 |
| float | 1e-4 | 1,35 | 1,35 | 16 | 48,83538818 à 48,83585739 |

Temps par phase (premier passage à 1e-4) : résolution linéaire 1,12 s en
double, 0,80 s en float ; prédiction 0,55 et 0,42 s ; assemblage de la
tangente 97 et 83 ms pour 16 assemblages (100 et 84 à 85 ms avant la
correction).

- Les trois calculs float convergent, avec les mêmes itérations que le
  double. Avant la correction, un sur trois divergeait.
- Le float est 1,37 fois plus rapide que le double (1,31 avant la
  correction, sur les calculs qui convergeaient).
- La norme des déplacements en float s'écarte de 2,7e-5 à 3,6e-5 du double
  au même critère, et de 3,4e-5 au plus du double convergé (3e-4 avant la
  correction ; le double au critère 1e-4 en est à 2,4e-6).
- Les calculs en double diffèrent entre eux de 1,3e-7, les calculs en
  float de 9,6e-6 : l'assemblage atomique ne garantit pas l'ordre des
  additions (l'assemblage ordonné,
  `MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered`, donne des calculs
  identiques au bit).
