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
  synchrone (par défaut) ou asynchrone.
- **mfem-mgis** n'est pas encore modifié. Son plan est à discuter.

## Différences avec les versions officielles

| Bibliothèque | Base officielle | Ce qui change | Détails |
|---|---|---|---|
| TFEL/MFront | master `fec0ade99` | interface `generic-parallel` des lois et des modèles (6 backends, 14 configurations), interface des propriétés matériau étendue, code GPU dans TFEL et dans le code généré, corrections du code existant | [docs/tfel.md](docs/tfel.md) |
| MGIS | master `818f0f9` | lots de lois, de modèles et de propriétés matériau, politique d'exécution, emplacement des tableaux, ressources mémoire, corrections | [docs/mgis.md](docs/mgis.md) |
| mfem-mgis | master `4c8a91d815` | rien pour l'instant | [docs/plan-mfem-mgis.md](docs/plan-mfem-mgis.md) |

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
   avec les questions à trancher.
7. [JOURNAL.md](JOURNAL.md) : cheminement jour par jour.

La documentation des interfaces, en anglais, est aussi dans l'arbre de
chaque bibliothèque : `tfel/docs/web/generic-parallel-*.md` pour TFEL,
`mgis/docs/web/` pour MGIS.

## Contenu du dépôt

| Chemin | Contenu |
|---|---|
| `docs/` | documentation du prototype ; `docs/historique/` garde les documents de travail dépassés (relecture du 06/10, conception initiale de MGIS) |
| `JOURNAL.md` | journal daté du travail |
| `tfel/` | sources complètes de TFEL modifié, prêtes à compiler (copie du dernier commit, sans l'historique de TFEL) |
| `patches/tfel/` | les mêmes changements en série de patchs depuis le master de TFEL (`git am`) |
| `patches/mgis/` | changements de MGIS en série de patchs depuis son master |
| `build-tfel.sh` | configuration, construction et installation de TFEL |
| `scripts/` | environnement (`env.sh`), synchronisation des patchs (`sync-tfel.sh`, `sync-mgis.sh`), matrice des tests, balayages de compilation et d'exécution, contrôle des arrondis, garde-fou mémoire des longues compilations |
| `scratch/` | outils de développement suivis : validation complète (`validate3.sh`), tests avec l'arbre de développement (`dev-test.sh`, `targeted-test.sh`), configurations de développement (`b2/`) |
| `reproducers/` | petits programmes autonomes qui vérifient un point précis |
| `results/` | résultats des matrices et des balayages ; défauts de compilateurs et leurs reproducteurs (`results/nvhpc-bugs/`) |
| `deps/spack/` | environnement spack du prototype (`mmm-gpu`) et recettes complémentaires |
| `deps/spack-sycl/` | spécification du LLVM de l'ancien environnement |

Les dépôts de travail `tfel-dev/` (branche de développement de TFEL),
`mgis/` et `mfem-mgis/` sont des clones locaux, non suivis par ce dépôt.

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
git checkout -b generic-parallel fec0ade99
git am ../patches/tfel/*.patch
```

MGIS, depuis le dépôt officiel, contre ce TFEL :

```bash
git clone <dépôt officiel de MGIS> mgis && cd mgis
git checkout -b gpu-proto 818f0f9
git am ../patches/mgis/*.patch
cd .. && source scripts/env.sh
cmake -S mgis -B build/mgis -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$PROTO_TFEL_INSTALL" -DCMAKE_CXX_COMPILER=g++-14 \
  -Denable-generic-parallel-cuda-tests=ON \
  -DCMAKE_CUDA_COMPILER="$(which nvcc)" -DCMAKE_CUDA_HOST_COMPILER=g++-14
cmake --build build/mgis --target check -j16   # construit et lance les tests
```

Les tests et leur lecture sont décrits dans
[docs/validation.md](docs/validation.md).
