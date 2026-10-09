# Environnement spack du prototype

Contenu :

- `mmm-gpu/spack.yaml` : environnement. NVIDIA HPC SDK 25.11 et son CUDA
  13.0.2 en externes ; Kokkos 5.2.2, HIP 7.2.3, LLVM 20.1.8, AdaptiveCpp
  25.10.0, TBB 2023.1.0 (chaînes de TFEL) ; MFEM 4.10.0 avec CUDA,
  OpenMP et libCEED, et les
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
  - `openmpi` : hwloc interne sans OpenCL (voir [MPI](#mpi)) ;
  - `mfem` : une construction partielle reprend sans tout recompiler
    (`make config` réécrit `config/config.mk`, dont dépendent tous les
    objets : sa date est gardée quand son contenu ne change pas) ; patch
    `mfem-4.10-convection-ea-cuda13.patch` : avec CUDA 13 (sm_120), cicc ne
    finit pas `fem/integ/bilininteg_convection_ea.cpp` (plus de 2 h 30 et
    17 Go, puis tué faute de mémoire) ; les noyaux 3D de l'assemblage par
    élément de la convection ne sont plus spécialisés au-delà de l'ordre 3,
    le noyau générique du même fichier prend les ordres plus élevés ;
  - `petsc` : copie de la recette de spack. La variante `cuda_lambda` de
    Kokkos n'existe plus en Kokkos 5 (les lambdas CUDA sont toujours
    actives) ;
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
fois, 4 tâches chacun), puis MFEM seul avec 2 tâches. Chaque phase tourne
sous `scripts/memory-governor.sh` : quand la mémoire disponible est trop
basse, il suspend tous les processus de la construction sauf le plus gros
(et ses ancêtres), qui peut ainsi finir et libérer sa mémoire ; tout
reprend ensuite. `ptxas` prend environ 22 Go pour les noyaux de
`fem/lor/lor_batched.cpp` (sm_120), les autres fichiers jusqu'à 8 Go.
Chaque phase est relancée jusqu'à trois fois en cas d'échec (manque de
mémoire), MFEM avec une seule tâche ; une construction interrompue reprend
où elle s'est arrêtée (`--dont-restage` et recette `mfem` de `mmm_gpu`).

## Préfixe CUDA

Dans le préfixe, `bin` est un dossier de liens, pas un lien : nvcc cherche
ses en-têtes dans `bin/../targets`. Avec un lien, `bin/..` désigne le CUDA
du SDK, sans les en-têtes des bibliothèques mathématiques (`cublas_v2.h`
introuvable pour les paquets qui ne passent pas le dossier des en-têtes
explicitement, comme Kokkos Kernels).

## MPI

OpenMPI est compilé avec CUDA (MPI compatible GPU), pour que le code
reste générique, même si la machine n'a qu'un GPU. Il utilise son hwloc
interne, construit sans son greffon OpenCL : avec CUDA, ce hwloc trouve
la bibliothèque OpenCL du toolkit et échoue à l'édition de liens (`cannot
determine absolute directory name of 'no/opencl/lib/'`). La recette
`openmpi` de `mmm_gpu` ajoute `--disable-opencl`. Le hwloc de spack avec
CUDA n'est pas utilisé : partagé avec LLVM et TBB, il les reconstruirait.

## État de l'installation du 09/10/2026

L'installation actuelle de MFEM (hachage `zxpdstv`) a été finie à la main
pour tenir les délais ; son contenu est celui que donnent les recettes,
mais pas encore par une construction directe :

- `fem/integ/bilininteg_convection_ea.cpp` a été compilé depuis la version
  patchée (même contenu que `mfem-4.10-convection-ea-cuda13.patch`), dans le
  dossier de construction ;
- la bibliothèque statique de gslib (`libgs.a`, dans le préfixe de gslib) a
  été recompilée avec `-fPIC` (mêmes sources, mêmes options) : non
  relogeable, elle ne pouvait pas entrer dans `libmfem.so`. L'originale est
  gardée hors du dépôt.

À faire : imposer gslib `+shared` dans `mmm-gpu/spack.yaml`, reconcrétiser
et reconstruire MFEM (environ 3 h), pour une installation entièrement
produite par les recettes.

## Contraintes entre paquets

- PETSc avec `+openmp` : les en-têtes d'un Kokkos compilé avec OpenMP
  exigent OpenMP chez le compilateur (`#error "You enabled Kokkos OpenMP
  support without enabling OpenMP in the compiler!"`).
- MFEM avec `cxxstd=20` : Umpire 2026 (par camp) utilise les concepts du
  C++20 ; MFEM se compile en C++17 par défaut.
- kokkos-kernels partagé exige Kokkos `+pic` : Kokkos a été reconstruit
  avec cette variante (même version, mêmes options sinon).

## Pilotes

CUDA 13.0 demande un pilote 580 ou plus récent. Avec un pilote plus
ancien, utiliser le CUDA 12.9 du même SDK (`make-nvhpc-cuda-prefix.sh
... 12.9 ...` et externe `cuda@12.9.1`).
