# Backend kokkos sur le polycristal (11/10/2026)

Script : `cases/measure-backends-polycrystal.sh`, avec les variables
`polycristal=build/mm-opera-hpc-kokkos/polycrystal`,
`cas=build/cas-kokkos` et `backends="cuda:cuda kokkos:kokkos"`.

Polycristal de mm-opera-hpc (5 grains, loi de monocristal de l'UO2,
orthotrope, grandes transformations, implicite), 10 pas, un processus et
le GPU, tout par noyaux. Même exécutable pour tous les calculs : celui de
mm-opera-hpc construit contre l'installation `install/mfem-mgis-kokkos`
(branche kokkos de mfem-mgis, aed33abb8a, option `enable-kokkos`). La
ligne « mfem-mgis » de `summary.txt` donne la version du sous-module, pas
celle de cette installation.

Lois comparées :

- `generic` : interface generic, intégrée sur l'hôte (référence des
  contraintes macroscopiques) ;
- cuda : interface generic-parallel, configuration `config-cuda-nvcc.json` ;
- kokkos : interface generic-parallel, configuration `config-kokkos.json`
  (Kokkos 5.2.2 construit pour CUDA, compilé par `nvcc_wrapper`). Kokkos
  est initialisé par mfem-mgis (variable `MFEM_MGIS_KOKKOS=1`), sur le GPU
  du périphérique de MFEM.

Machine calme (navigateur seul), deux passages.

Contraintes macroscopiques : écart relatif au plus 2,8e-13 par rapport au
premier calcul avec la loi `generic` (70 valeurs). Mêmes itérations de
Newton (98) dans les six calculs.

| Bibliothèque | Total (s) | Intégration, 158 appels (s) | Résolution (s) |
|---|---|---|---|
| `generic` (hôte) | 107,65 et 107,48 | 99,97 et 99,84 | 5,73 et 5,66 |
| cuda (nvcc) | 12,88 et 12,79 | 5,61 et 5,61 | 5,64 et 5,56 |
| kokkos (nvcc) | 12,97 et 13,22 | 5,62 et 5,69 | 5,70 et 5,88 |

- kokkos et cuda donnent le même temps d'intégration, à 1,4 % près. Les
  deux chaînes passent par nvcc. Comme pour cuda et hip
  (`results/backends-polycrystal-2026-10-11`), c'est le compilateur qui
  compte, pas le modèle de programmation.
- L'écart sur le total (0,1 et 0,4 s) vient surtout de la résolution, qui
  ne dépend pas de la loi. L'initialisation de Kokkos, faite par
  `mfem_mgis::initialize`, précède la résolution du problème et n'est pas
  comptée dans le total.
- Il fallait que le programme initialise Kokkos. C'est maintenant le rôle
  de mfem-mgis : aucune modification de mm-opera-hpc n'a été nécessaire.
