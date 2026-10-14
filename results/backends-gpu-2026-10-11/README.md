# Backends GPU de l'interface generic-parallel dans mfem-mgis (11/10/2026)

Script : `cases/measure-backends-gpu.sh`. Poutre r4 (111 843 inconnues, 10
pas) sur un processus et le GPU, tout par noyaux, tangente assemblée,
Newton et PCG + BoomerAMG ; loi IsotropicLinearHardeningPlasticity générée
par l'interface generic-parallel avec chaque configuration de TFEL pour le
GPU. mfem-mgis 517e33cde4 (`build/cas-autres` ; la ligne « mfem-mgis » de
`summary.txt` donne la version des sources, pas celle de la bibliothèque
liée). Bibliothèques de la mesure générées à la main avec les mêmes
configurations (variable `libs` du script), puis construites par CMake
dans `build/cas-autres/poutre/generic-parallel-*`. Machine calme, deux
passages.

## Quels backends fonctionnent dans mfem-mgis

| Configuration | Résultat |
|---|---|
| cuda-nvcc | fonctionne |
| cuda-clang | fonctionne |
| hip-clang-nvidia-platform | fonctionne |
| hip-hipcc | fonctionne |
| kokkos | fonctionne si le programme initialise Kokkos (option `--kokkos` de la poutre, `cases/common/KokkosRuntime`) ; sinon « Kokkos is not initialized » |
| sycl-default-acpp | refusé : « no SYCL queue given », le programme doit fournir une file SYCL dans les options d'exécution |
| stlpar-parunseq-nvhpc-gpu | échoue (« libgomp: TODO ») : le support d'exécution d'OpenMP de nvhpc entre en conflit avec libgomp, chargé par MFEM et mfem-mgis compilés par GCC |
| openmp-offload-nvhpc-gpu | idem |

Les cinq backends qui fonctionnent donnent la norme attendue au test court
(r1, 3,264041277, celle de la loi générée par l'interface `generic`) et,
en r4, 48,83704586, comme les mesures précédentes de la poutre avec Newton
et PCG + BoomerAMG (48,83704551 à 48,83704606 selon le nombre de processus
et les réglages, tolérance du solveur itératif). Mêmes itérations de Newton (32 ; une permutation de deux
pas pour kokkos au premier passage, même total).

## Temps (secondes)

| Configuration | Boucle sur les pas (deux passages) | Intégration (42 appels) |
|---|---|---|
| cuda-nvcc | 3,04 et 3,28 | 0,05 |
| cuda-clang | 3,15 et 3,08 | 0,05 |
| hip-clang-nvidia-platform | 3,10 et 3,42 | 0,05 à 0,06 |
| hip-hipcc | 3,10 et 3,12 | 0,05 |
| kokkos | 3,17 et 3,09 | 0,04 |

- Sur ce cas, la loi coûte environ 1 ms par intégration (un peu moins avec
  kokkos : 0,039 et 0,041 s pour les 42 intégrations, contre 0,050 à
  0,055 s, aux deux passages) : les écarts de la boucle sur les pas
  (jusqu'à 13 %, de 3,04 à 3,42 s) sont du bruit, la résolution linéaire
  dominant (2,2 à 2,5 s).
- Pour départager les backends sur l'intégration, il faut une loi coûteuse
  (polycristal, ex8).
