# Environnement spack du prototype

Contenu :

- `mmm-gpu/spack.yaml` : environnement (NVIDIA HPC SDK 25.11 et son CUDA
  13.0.2 en externes, Kokkos 5.2.2, Ginkgo 1.11.0, MFEM 4.10.0, HIP 7.2.3,
  LLVM 20.1.8, AdaptiveCpp 25.10.0, TBB 2023.1.0) ;
- `repos/spack_repo/mmm_gpu` : dépôt de paquets complémentaire, avec la
  recette d'AdaptiveCpp 25.10 (le paquet `hipsycl` de spack s'arrête à
  24.10, qui ne compile pas avec CUDA 13) ;
- `make-nvhpc-cuda-prefix.sh` : préfixe CUDA à la disposition du toolkit
  autonome, construit à partir du SDK (dont les bibliothèques
  mathématiques sont dans un dossier séparé).

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
spack-mmm/bin/spack -e deps/spack/mmm-gpu install -j8
```

`SPACK_DISABLE_LOCAL_CONFIG` isole l'environnement de la configuration de
l'utilisateur (`~/.spack`). Sur une machine de 30 Go de mémoire, les
compilations de LLVM et des paquets CUDA peuvent manquer de mémoire :
`scripts/memory-governor.sh <pid>` met l'installation en pause quand la
mémoire disponible est trop basse.

## Pilotes

CUDA 13.0 demande un pilote 580 ou plus récent. Avec un pilote plus
ancien, utiliser le CUDA 12.9 du même SDK (`make-nvhpc-cuda-prefix.sh
... 12.9 ...` et externe `cuda@12.9.1`).
