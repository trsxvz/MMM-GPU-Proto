# Backends GPU sur le polycristal (loi coûteuse, 11/10/2026)

Script : `cases/measure-backends-polycrystal.sh`. Polycristal de mm-opera-hpc
(5 grains, loi de monocristal de l'UO2, orthotrope, grandes
transformations, implicite), 10 pas, un processus et le GPU, tout par
noyaux. Loi intégrée sur l'hôte (interface `generic`, référence) ou sur le
GPU avec la bibliothèque générée par l'interface `generic-parallel` pour
chaque backend qui fonctionne sans modifier le programme
(`build/cas-autres/polycristal`). Exécutable de mm-opera-hpc de
`build/mm-opera-hpc`, construit le 10/10 au soir et lié statiquement à
mfem-mgis (de2e76a913 ou antérieur ; la ligne « mfem-mgis » de
`summary.txt` donne la version des sources). Machine calme, deux passages.

Contraintes macroscopiques : écart relatif au plus 2,4e-13 par rapport au
premier calcul avec la loi `generic` (70 valeurs). Mêmes itérations de
Newton (98) dans les dix calculs.

| Bibliothèque | Total | Intégration (158 appels) | Résolution |
|---|---|---|---|
| `generic` (hôte) | 109,73 et 108,67 | 102,18 et 101,19 | 5,53 et 5,52 |
| cuda (nvcc) | 12,81 et 12,83 | 5,57 et 5,58 | 5,53 et 5,54 |
| hip (hipcc, plateforme NVIDIA) | 12,86 et 12,84 | 5,58 et 5,58 | 5,57 et 5,56 |
| cuda (clang) | 14,92 et 14,87 | 7,67 et 7,65 | 5,56 et 5,51 |
| hip (clang, plateforme NVIDIA) | 14,90 et 14,98 | 7,67 et 7,68 | 5,54 et 5,58 |

- Sur le GPU, l'intégration est 18 fois plus rapide (13 fois avec clang)
  que sur un cœur de l'hôte (un processus, loi `generic`).
- Les deux chaînes qui passent par nvcc (cuda par nvcc, hip par hipcc, qui
  appelle nvcc sur la plateforme NVIDIA) donnent le même temps ; les deux
  qui passent par clang (cuda et hip) prennent 38 % de temps en plus sur
  l'intégration (nvcc 27 % de moins que clang), et 16 % de bout en bout. Le langage de la source (cuda ou
  hip) ne compte pas, le compilateur oui.
- kokkos n'est pas mesuré ici : le programme de mm-opera-hpc
  n'initialise pas Kokkos (voir `results/backends-gpu-2026-10-11`). Il
  l'est dans `results/kokkos-polycrystal-2026-10-11`, avec mfem-mgis
  construit avec Kokkos (option `enable-kokkos`), qui l'initialise.
