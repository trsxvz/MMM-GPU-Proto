# MMM-GPU-Proto

Prototype de la chaîne MFront → MGIS → mfem-mgis sur GPU et CPU
parallèles. Première étape (en cours de finalisation) : TFEL/MFront.

- **Interface `generic-parallel` des lois et des modèles** : MFront génère,
  pour chaque loi et chaque hypothèse de modélisation, une fonction qui
  intègre tout un lot de points d'intégration.
  - Couverture : petites et grandes déformations, fonctions
    d'initialisation, post-traitements, rotations des lois orthotropes,
    monocristaux, variables de comportement, brique
    StandardElastoViscoPlasticity, solveur TDLS.
  - Backends : CUDA, HIP, algorithmes parallèles de la bibliothèque
    standard (g++/TBB, nvc++ sur GPU), Kokkos, SYCL (AdaptiveCpp) et
    OpenMP.
- **Interface `generic-parallel` des propriétés matériau** (existante dans
  TFEL) : backends Kokkos et OpenMP ajoutés.
- **Tests** : comparaison systématique avec l'interface `generic`, qui sert
  de référence, pour chaque backend.

## Contenu

| Chemin | Contenu |
|---|---|
| `JOURNAL.md` | journal daté, plan cochable, état et objectifs : **point d'entrée pour reprendre** |
| `DESIGN.md` | conception de l'interface `generic-parallel` (données, exécution par point, backends) |
| `tfel/` | sources complètes de TFEL modifié, prêtes à compiler (copie du dernier commit de la branche de développement, sans l'historique de TFEL) |
| `patches/tfel/` | les mêmes modifications sous forme de série de patchs depuis le master `fec0ade99` (`git am`), et son README |
| `build-tfel.sh` | configuration, construction et installation de TFEL (chemins de la machine de développement) |
| `scripts/` | environnement (`env.sh`), mise à jour de `tfel/` et des patchs (`sync-tfel.sh`), matrice des backends, balayage de compilation de tous les tests MFront, contrôle des annotations hôte/device par clang CUDA |
| `results/` | résultats des matrices et des balayages, défauts de compilateurs rencontrés avec leurs reproducteurs (`results/nvhpc-bugs/`) |
| `scratch/` | scripts de validation, configurations de développement, proposition de conception pour MGIS (`d1-mgis-design.md`) |
| `deps/spack-sycl/` | spécification spack du LLVM utilisé par AdaptiveCpp |

La documentation utilisateur des interfaces est dans
`tfel/docs/web/generic-parallel-behaviour-interface.md` et
`tfel/docs/web/generic-parallel-material-property-interface.md`.

## Compiler

`build-tfel.sh` compile `tfel/` dans `build/tfel` et l'installe dans
`install/tfel`. Pour repartir du dépôt officiel de TFEL avec son historique,
appliquer les patchs :

```bash
git clone <dépôt officiel de TFEL> tfel
cd tfel
git checkout -b generic-parallel-behaviours fec0ade99
git am ../patches/tfel/*.patch
```

`build-tfel.sh` et `scripts/env.sh` décrivent les compilateurs et les
dépendances utilisés : nvcc 13.0, hipcc 5.7 sur la plateforme NVIDIA,
g++ 14, nvc++ 25.11, Kokkos 5.0.2, AdaptiveCpp 25.10 et TDLS.

## Bibliothèques suivantes

MGIS (`818f0f9`) et mfem-mgis (`4c8a91d815`) ne sont pas encore modifiés.
Ils seront adaptés une fois TFEL finalisé (voir `JOURNAL.md`, plan D1/D2).
