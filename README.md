# MMM-GPU-Proto

Prototype de la chaîne MFront (TFEL), MGIS et mfem-mgis, sur CPU et sur
GPU.

## En bref

- **TFEL/MFront** génère, pour chaque loi de comportement, modèle et
  propriété matériau, des fonctions qui traitent tout un lot de points
  d'intégration. Le même fichier MFront compile pour CUDA, HIP, Kokkos,
  SYCL, OpenMP et les algorithmes parallèles du C++. Les résultats sont
  comparés systématiquement à ceux de l'interface `generic`.
- **MGIS** charge ces bibliothèques et calcule des lots sans copier les
  données. Ses fonctions existantes basculent seules sur le calcul par
  lot. L'appelant choisit la mémoire, la file d'exécution, et le calcul
  synchrone (par défaut) ou asynchrone. Les données peuvent être rangées
  composante par composante, ou appartenir à l'appelant et être décrites
  par des vues, sur GPU compris, fin de pas comprise.
- **mfem-mgis** intègre les lois par lots, sur CPU et sur GPU. Les données
  aux points vivent dans la mémoire de MFEM. Gradients, forces
  intérieures, tangente (assemblée ou sans matrice), prédiction,
  thermique, micromorphe, FBar et post-traitements sont calculés par
  noyaux. Des solveurs adaptés au GPU s'ajoutent à Newton : point fixe
  accéléré, NewtonK0, Newton inexact, solveur direct cuDSS. Toute la suite
  de tests passe sur le GPU, sauf les deux tests qui exigent l'absence de
  GPU. Le plan suivi est découpé en lots
  ([docs/plan-mfem-mgis-detaille.md](docs/plan-mfem-mgis-detaille.md)).
- **mfem-mgis-examples** et **mm-opera-hpc** : tous les exemples et les cas
  (polycristal, cermet, bulle) tournent par lots et sur le GPU. Sur une
  RTX 5080, la poutre prend 1,96 s avec NewtonK0, contre 9,5 s sur 8 cœurs
  avec la même méthode (13,6 s avec Newton), et ex8 2,6 s contre 19,2 s
  sur 2 cœurs
  ([docs/vue-ensemble.md](docs/vue-ensemble.md#état-davancement)).
- Les changements ne sont pas encore proposés aux dépôts officiels. Leur
  découpage en contributions est décrit dans
  [docs/integration-amont.md](docs/integration-amont.md).

## Différences avec les versions officielles

| Bibliothèque | Base officielle | Ce qui change | Détails |
|---|---|---|---|
| TFEL/MFront | master `ff0959ece` | interface `generic-parallel` des lois et des modèles (6 backends, 14 configurations), interface des propriétés matériau étendue, code GPU dans TFEL et dans le code généré, corrections du code existant ; 116 patchs | [docs/tfel.md](docs/tfel.md) |
| MGIS | master `e8ef4b4` | lots de lois, de modèles et de propriétés matériau, politique d'exécution, emplacement des tableaux, ressources mémoire, vues et rangements, fin de pas sur GPU, corrections ; 68 patchs | [docs/mgis.md](docs/mgis.md) |
| mfem-mgis | master `8fac67336e` | intégration par lots sur CPU et GPU (données aux points dans la mémoire de MFEM), noyaux des gradients, des forces intérieures et de la tangente (assemblée ou sans matrice ; mécanique et thermique stationnaire), prédiction par noyaux, solveur `FixedPoint` (point fixe, accélérations de MTest, NewtonK0), Newton inexact, options de BoomerAMG, lisseurs sans matrice, réutilisation de l'analyse de MUMPS, solveur direct cuDSS sur le GPU, assemblage ordonné de la tangente (mêmes résultats à chaque passage sur le GPU), réglages d'exécution (`--execution cpu|gpu|mixed`, placement de chaque phase en mode `mixed` par `--finite-elements`, `--behaviours`, `--linear-solver-device` et `--post-processings`, placement de la loi de chaque matériau par le paramètre `BehaviourExecutionLocation`, copies hôte/GPU comptées, bibliothèques de lois pour le CPU et le GPU), noyaux compilés pour CUDA ou HIP (GPU AMD : compilation vérifiée, sans essai), corrections ; 154 patchs | [docs/plan-mfem-mgis-detaille.md](docs/plan-mfem-mgis-detaille.md), [JOURNAL.md](JOURNAL.md) |
| mfem-mgis-examples | master `600da70` (dépôt `latug0/mfem-mgis-examples`) | options des exemples (bibliothèque de lois, périphérique, solveurs, profilage), synchronisations sur l'hôte pour le GPU, choix des solveurs directs, lois pour le CPU et le GPU choisies par `--execution`, tests sur GPU ; 31 patchs | `patches/mfem-mgis-examples/` |
| mm-opera-hpc | main `49147f2` (dépôt `rprat-pro/mm-opera-hpc`) | polycristal : choix du solveur non linéaire (Newton, NewtonK0), de la tolérance adaptative et d'AmgX, loi intégrée sur le GPU ; bulle et cermet : options de mfem-mgis et de la loi, synchronisations sur l'hôte, profilage des post-traitements, solveurs directs, lois pour le CPU et le GPU choisies par `--execution`, tests sur GPU ; 26 patchs | `patches/mm-opera-hpc/` ; mesure du polycristal : [results/L5-polycristal-2026-10-10-0502/README.md](results/L5-polycristal-2026-10-10-0502/README.md) |

## Documentation

À lire dans cet ordre :

1. [docs/vue-ensemble.md](docs/vue-ensemble.md) : but, rôle de chaque
   bibliothèque, principes, vocabulaire, état d'avancement.
2. [docs/tfel.md](docs/tfel.md) : changements de TFEL et de MFront, et leurs
   justifications.
3. [docs/mgis.md](docs/mgis.md) : changements de MGIS.
4. [docs/validation.md](docs/validation.md) : tests, balayages, outils et
   résultats.
5. [docs/integration-amont.md](docs/integration-amont.md) : découpage en
   contributions pour les dépôts officiels.
6. [docs/plan-mfem-mgis.md](docs/plan-mfem-mgis.md) : plan pour mfem-mgis,
   avec les questions à trancher, et son découpage en lots
   ([docs/plan-mfem-mgis-detaille.md](docs/plan-mfem-mgis-detaille.md)).
7. [docs/plan-execution.md](docs/plan-execution.md) : choix de
   l'exécution sur CPU, sur GPU ou mixte (`--execution`), réalisé.
8. [docs/plan-precision.md](docs/plan-precision.md) : précision float ou
   double sur toute la chaîne, plan à valider.
9. [docs/plan-ci-site.md](docs/plan-ci-site.md) : intégration continue et
   site du prototype, plan à valider.
10. [docs/solveurs-gpu.md](docs/solveurs-gpu.md) : solveurs,
    préconditionneurs et techniques CPU et GPU retenus pour l'étude de
    mfem-mgis.
11. [JOURNAL.md](JOURNAL.md) : cheminement jour par jour.

La documentation des interfaces, en anglais, est aussi dans l'arbre de
chaque bibliothèque : `tfel/docs/web/generic-parallel-*.md` pour TFEL,
`mgis/docs/web/` pour MGIS.

## Contenu du dépôt

| Chemin | Contenu |
|---|---|
| `docs/` | documentation du prototype ; `docs/historique/` garde les documents de travail dépassés (relecture du 06/10, conception initiale de MGIS) |
| `JOURNAL.md` | journal daté du travail |
| `tfel/` | sources complètes de TFEL modifié, prêtes à compiler (copie du dernier commit, sans l'historique de TFEL) |
| `patches/tfel/` | les mêmes changements en série de patchs depuis le master de TFEL (`git am --keep-cr`) |
| `patches/mgis/` | changements de MGIS en série de patchs depuis son master |
| `mfem-mgis/`, `mfem-mgis-examples/`, `mm-opera-hpc/` | sources de ces dépôts (copie du dernier commit, sans leur historique) |
| `patches/mfem-mgis/`, `patches/mfem-mgis-examples/`, `patches/mm-opera-hpc/` | leurs changements en séries de patchs depuis leur base |
| `build-tfel.sh` | configuration, construction et installation de TFEL |
| `build-mgis.sh`, `build-mfem-mgis.sh` | construction et tests de MGIS, puis de mfem-mgis, des exemples, de mm-opera-hpc et des cas de `cas/`, contre le TFEL installé par `build-tfel.sh` |
| `cas/` | cas de mesure (poutre, thermique, polycristal, bulle, cermet, ssna303, exemples) et scripts de mesure `mesures-*.sh`, voir `cas/README.md` |
| `scripts/` | environnement (`env.sh`), synchronisation des sources et des patchs (`sync-tfel.sh`, `sync-mgis.sh`, `sync-repository.sh`) et leur vérification (`verify-patches.sh`), matrice des tests, suite de mfem-mgis sur le GPU (`mfem-mgis-gpu-suite.sh`), construction pour un GPU AMD (`build-hip-amd.sh`), en-têtes GPU de MGIS (`check-mgis-gpu-headers.sh`), annotations hôte et GPU des en-têtes de TFEL (`clang-cuda-header-scan.sh`, `annotate-function-templates.py`), balayages de compilation et d'exécution, contrôle des arrondis, garde-fou mémoire des longues compilations |
| `scratch/` | outils de développement suivis : validation complète (`validate3.sh`), tests avec l'arbre de développement (`dev-test.sh`, `targeted-test.sh`), configurations de développement (`b2/`) |
| `reproducers/` | petits programmes autonomes qui vérifient un point précis |
| `results/` | résultats des matrices et des balayages ; défauts de compilateurs et leurs reproducteurs (`results/nvhpc-bugs/`) |
| `deps/spack/` | environnement spack du prototype (`mmm-gpu`) et recettes complémentaires |
| `deps/spack-sycl/` | spécification du LLVM de l'ancien environnement |

Les dépôts de travail sont des clones locaux. `tfel/`, `mfem-mgis/`,
`mfem-mgis-examples/` et `mm-opera-hpc/` sont aussi copiés ici, sans leur
historique, par `scripts/sync-tfel.sh` et `scripts/sync-repository.sh`.
`tfel-dev/` (branche de développement de TFEL) et `mgis/` ne sont pas
copiés : MGIS est suivi par sa série de patchs.

## Construire

Les chemins des scripts sont ceux de la machine de développement.
`scripts/env.sh` décrit l'environnement : compilateurs, Kokkos,
AdaptiveCpp, HIP, CUDA.

TFEL, depuis les sources de ce dépôt :

```bash
./build-tfel.sh install   # construit tfel/ et l'installe dans install/tfel-mmm-gpu
```

TFEL, depuis le dépôt officiel :

```bash
git clone <dépôt officiel de TFEL> tfel && cd tfel
git checkout -b generic-parallel ff0959ece
git am --keep-cr ../patches/tfel/*.patch   # certains fichiers de TFEL sont en CRLF
```

MGIS, depuis le dépôt officiel, contre ce TFEL :

```bash
git clone <dépôt officiel de MGIS> mgis && cd mgis
git checkout -b gpu-proto e8ef4b4
git am ../patches/mgis/*.patch
cd .. && source scripts/env.sh
cmake -S mgis -B build/mgis -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PROTO_TFEL_INSTALL" -DCMAKE_CXX_COMPILER=g++-14 \
  -Denable-generic-parallel-cuda-tests=ON \
  -DCMAKE_CUDA_COMPILER="$(which nvcc)" -DCMAKE_CUDA_HOST_COMPILER=g++-14
cmake --build build/mgis --target check -j16   # construit et lance les tests
```

Les scripts de la racine font ces étapes avec les dépôts locaux :
`./build-mgis.sh` construit MGIS (dossier `mgis/`), le teste et
l'installe dans `install/mgis` ; `./build-mfem-mgis.sh` construit, teste
et installe mfem-mgis, puis construit et teste les exemples, mm-opera-hpc
et les cas de `cas/` (voir l'en-tête de chaque script). Les dépendances
(MFEM, hypre, solveurs, compilateurs) viennent de l'environnement spack
`mmm-gpu`, décrit dans [deps/spack/README.md](deps/spack/README.md) et
installé par `deps/spack/install-mmm-gpu.sh`.

Les tests et leur lecture sont décrits dans
[docs/validation.md](docs/validation.md).

Compilation pour un GPU AMD (compilation et édition de liens seulement :
le GPU AMD intégré de la machine de développement n'est accessible
qu'aux membres du groupe render) : `scripts/build-hip-amd.sh`
construit hypre (pour l'hôte), MFEM compilé avec HIP pour gfx90a,
mfem-mgis et les exemples avec ce MFEM, et lance les tests qui ne
demandent pas de GPU. ROCm est installé sans droits d'administrateur
(paquets Python de ROCm, voir l'en-tête du script).
