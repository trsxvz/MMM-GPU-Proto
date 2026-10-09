# Environnement spack du prototype

Contenu :

- `mmm-gpu/spack.yaml` : environnement. NVIDIA HPC SDK 25.11 et son CUDA
  13.0.2 en externes ; Kokkos 5.2.2, HIP 7.2.3, LLVM 20.1.8, AdaptiveCpp
  25.10.0, TBB 2023.1.0 (chaînes de TFEL) ; MFEM 4.9.0 avec CUDA et les
  bibliothèques de solveurs de l'étude (hypre 3.2.0, Ginkgo 1.11.0 avec
  MPI, PETSc 3.26.0 avec Kokkos, SUNDIALS, AmgX, SuperLU_DIST, STRUMPACK,
  MUMPS, UMFPACK), cuDSS 0.7.1, OpenMPI 5, OpenBLAS. Le choix des
  solveurs est expliqué dans [docs/solveurs-gpu.md](../../docs/solveurs-gpu.md) ;
- `repos/spack_repo/mmm_gpu` : dépôt de paquets complémentaire, prioritaire
  sur celui de spack :
  - `adaptivecpp` : AdaptiveCpp 25.10 (le paquet `hipsycl` de spack
    s'arrête à 24.10, qui ne compile pas avec CUDA 13) ;
  - `amgx` : version `2.5.0.20260701`, commit de la branche principale
    avec deux corrections complémentaires pour CUDA 13, postérieures à la
    version 2.5.0 (C++17 avec CUDA 13.1 et plus, tuples de CUDA 13) ;
  - `petsc` : copie de la recette de spack, avec deux changements. La
    variante `cuda_lambda` de Kokkos n'existe plus en Kokkos 5 (les
    lambdas CUDA sont toujours actives) ; un OpenMPI sans CUDA est
    accepté ;
- `make-nvhpc-cuda-prefix.sh` : préfixe CUDA à la disposition du toolkit
  autonome, construit à partir du SDK (dont les bibliothèques
  mathématiques sont dans un dossier séparé) ;
- `install-mmm-gpu.sh` : installation en deux temps, adaptée à une machine
  de 30 Go de mémoire (voir plus bas).

## Recréer l'environnement

```bash
git clone --depth 1 --branch v1.2.2 https://github.com/spack/spack.git spack-mmm
git clone https://github.com/spack/spack-packages.git spack-packages-mmm
git -C spack-packages-mmm checkout 3b7fdffc24
deps/spack/make-nvhpc-cuda-prefix.sh <sdk>/Linux_x86_64/25.11 13.0 <prefix>
```

Dans `mmm-gpu/spack.yaml`, adapter :

- le chemin du dépôt `builtin` (clone de spack-packages) ;
- les préfixes de `nvhpc` (racine du SDK) et de `cuda` (préfixe construit
  ci-dessus) ;
- le compilateur `gcc` externe ;
- `cuda_arch` : 80 pour A100, 120 pour Blackwell (RTX 5080). Kokkos
  n'accepte qu'une architecture par construction.

Puis :

```bash
export SPACK_DISABLE_LOCAL_CONFIG=1
spack-mmm/bin/spack -e deps/spack/mmm-gpu concretize -f
deps/spack/install-mmm-gpu.sh spack-mmm/bin/spack
```

`SPACK_DISABLE_LOCAL_CONFIG` isole l'environnement de la configuration de
l'utilisateur (`~/.spack`). L'environnement a aussi son propre cache
(`misc_cache`) : le cache de l'utilisateur est partagé avec les autres
installations de spack de la machine, dont le dépôt `builtin` peut être à
un autre commit, ce qui casse l'index des patchs.

`install-mmm-gpu.sh` installe d'abord tout sauf MFEM (deux paquets à la
fois, 4 tâches chacun), puis MFEM seul avec 2 tâches : chaque processus
`ptxas` des noyaux CUDA de MFEM peut prendre 10 Go. Chaque phase tourne
sous `scripts/memory-governor.sh`, qui met la construction en pause quand
la mémoire disponible est trop basse.

## Préfixe CUDA

Dans le préfixe, `bin` est un dossier de liens, pas un lien : nvcc cherche
ses en-têtes dans `bin/../targets`. Avec un lien, `bin/..` désigne le CUDA
du SDK, sans les en-têtes des bibliothèques mathématiques (`cublas_v2.h`
introuvable pour les paquets qui ne passent pas le dossier des en-têtes
explicitement, comme Kokkos Kernels).

## MPI

OpenMPI est compilé sans CUDA. Avec CUDA, il lui faut soit son hwloc
interne, qui ne se compile pas avec CUDA (chemin d'OpenCL invalide), soit
le hwloc de spack avec CUDA, dont le changement reconstruit LLVM et TBB,
utilisés par TFEL. Avec un seul GPU, un MPI compatible GPU n'apporte rien.
PETSc doit être lancé avec l'option `-use_gpu_aware_mpi 0`.

## Pilotes

CUDA 13.0 demande un pilote 580 ou plus récent. Avec un pilote plus
ancien, utiliser le CUDA 12.9 du même SDK (`make-nvhpc-cuda-prefix.sh
... 12.9 ...` et externe `cuda@12.9.1`).
