# Backends GPU sur le polycristal (loi coûteuse, 11/10/2026)

Script : `cas/mesures-backends-polycristal.sh`. Polycristal de mm-opera-hpc
(5 grains, loi de monocristal de l'UO2, orthotrope, grandes
transformations, implicite), 10 pas, un processus et le GPU, tout par
noyaux. Loi intégrée sur l'hôte (interface `generic`, référence) ou sur le
GPU avec la bibliothèque générée par l'interface `generic-parallel` pour
chaque backend qui fonctionne sans modifier le programme
(`build/cas-autres/polycristal`). Exécutable de mm-opera-hpc de
`build/mm-opera-hpc` (mfem-mgis installé). Machine calme, deux passages.

Contraintes macroscopiques : écart relatif au plus 2,4e-13 entre tous les
calculs (70 valeurs). Mêmes itérations de Newton (98).

| Bibliothèque | Total | Intégration (158 appels) | Résolution |
|---|---|---|---|
| `generic` (hôte) | 109,73 et 108,67 | 102,18 et 101,19 | 5,53 |
| cuda (nvcc) | 12,81 et 12,83 | 5,57 et 5,58 | 5,53 |
| hip (hipcc, plateforme NVIDIA) | 12,86 et 12,84 | 5,58 et 5,58 | 5,56 |
| cuda (clang) | 14,92 et 14,87 | 7,67 et 7,65 | 5,56 |
| hip (clang, plateforme NVIDIA) | 14,90 et 14,98 | 7,67 et 7,68 | 5,56 |

- Sur le GPU, l'intégration est 18 fois plus rapide que sur l'hôte.
- Les deux chaînes qui passent par nvcc (cuda par nvcc, hip par hipcc, qui
  appelle nvcc sur la plateforme NVIDIA) donnent le même temps ; les deux
  qui passent par clang (cuda et hip) sont 38 % plus lentes sur
  l'intégration, et 16 % de bout en bout. Le langage de la source (cuda ou
  hip) ne compte pas, le compilateur oui.
- kokkos n'est pas mesuré ici : le programme de mm-opera-hpc
  n'initialise pas Kokkos (voir `results/backends-gpu-2026-10-11`).
