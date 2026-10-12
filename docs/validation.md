# Validation

Ce document décrit comment les changements de TFEL et de MGIS sont
vérifiés, avec quels outils, et où trouver les résultats.

## Principe

L'interface `generic` sert de référence. Pour chaque loi testée, la même
source MFront est compilée deux fois : avec l'interface `generic` par le
compilateur de l'hôte, et avec l'interface `generic-parallel` pour la
configuration testée. Les deux bibliothèques calculent les mêmes points, et
les résultats sont comparés. Les écarts admis sont des écarts d'arrondi.

Un test qui ne peut pas échouer ne vaut rien. Les tests importants sont
donc vérifiés par une version volontairement fausse du code (indice relatif
au lieu d'absolu, pas d'un champ faux, flux par défaut au lieu du flux de
l'appelant, etc.) : le test doit la détecter.

## Environnement

- Machine : GPU NVIDIA RTX 5080 (architecture `sm_120`), sans GPU AMD.
- Environnement spack `mmm-gpu` (`deps/spack/`) : NVIDIA HPC SDK 25.11
  (nvc++ et CUDA 13.0.2), LLVM 20.1.8, AdaptiveCpp 25.10, Kokkos 5.2.2,
  HIP 7.2.3 (plateforme NVIDIA), hypre et Ginkgo. Voir
  `deps/spack/README.md`.
- `scripts/env.sh` donne l'environnement complet de tous les scripts.
- La gestion hétérogène de la mémoire (HMM) des pilotes NVIDIA permet au
  GPU de lire la mémoire ordinaire de l'hôte. Elle masque donc les erreurs
  d'emplacement. Les validations finales se font sans HMM
  (`uvm_disable_hmm=1`).

## TFEL

### Tests de l'interface des lois et des modèles

`tfel/test-suite/behaviours/generic-parallel/<configuration>/` : un dossier
par configuration (voir [tfel.md](tfel.md#configurations-et-compilateurs)).
Chaque test compile les lois de test avec les deux interfaces, puis
l'outil `compare.cxx` les compare sur des lots variés :

- trois types de calcul (intégration seule, opérateur tangent, prédiction)
  et deux rangements mémoire ;
- valeurs uniformes et valeurs par point, données groupe par groupe ou
  variable par variable ;
- fonctions d'initialisation, post-traitements, rotations, modèles,
  paramètres modifiés à l'exécution, échecs volontaires (indice du premier
  point en échec, statuts par point), calcul asynchrone ;
- solveur TDLS dans les lois implicites (sous-test `tdls`).

### Tests de l'interface des propriétés matériau

`tfel/test-suite/material-properties/generic-parallel/<configuration>/` :
valeurs comparées aux valeurs attendues, arguments uniformes et avec pas,
bornes et politiques, propriétés appelant d'autres propriétés ou définies
par des données, options d'exécution (flux CUDA et HIP non bloquant,
instance Kokkos, file SYCL absente), évaluation asynchrone.

Le test du flux est discriminant : un noyau lent écrit les arguments sur le
flux de l'appelant juste avant l'évaluation. Les valeurs ne sont justes que
si la propriété attend ce noyau. Avec le flux par défaut, le test échoue.

### Balayages

Les lois de test de MFront (`tfel/mfront/tests`) sont bien plus nombreuses
que les lois des tests de l'interface. Des balayages les couvrent toutes :

| Outil | Ce qui est vérifié | Résultat |
|---|---|---|
| `scripts/compile-sweep.sh` | compilation de toutes les lois, modèles et propriétés de test, pour nvcc, clang CUDA (qui refuse tout appel d'une fonction hôte depuis le code GPU), AdaptiveCpp et nvc++ | 541 compilations réussies, aucun code utilisateur bloquant |
| `scripts/runtime-sweep.py` | exécution de toutes les lois et modèles de test, comparée à l'interface `generic`, avec trois amplitudes de chargement | tous les écarts restants sont des écarts d'arrondi, vérifiés |
| `scripts/property-sweep.py` | exécution des 56 propriétés matériau de test | conformes |
| `scripts/roundoff-check.sh`, `scripts/roundoff-triage.py` | rejoue une loi sur l'hôte avec et sans contraction des opérations (FMA), comme nvcc sur GPU | les écarts signalés viennent de lois qui amplifient l'arrondi |
| `scripts/clang-cuda-header-scan.sh` | compilation de chaque en-tête de TFEL par clang CUDA | annotations hôte et GPU cohérentes |

### Validation complète

`scratch/validate3.sh <étiquette>` enchaîne :

1. la reconfiguration et l'installation de TFEL ;
2. la suite de tests complète de TFEL (`ctest`) : 6965 tests, tous verts ;
3. la matrice des tests `generic-parallel` à `-O3` pour toutes les
   configurations (`scripts/run-generic-parallel-tests.sh`) ;
4. les balayages de compilation avec nvcc et clang.

Les résultats des matrices sont dans `results/generic-parallel-tests-*.txt`,
ceux des balayages dans `results/compile-sweep-*.tsv`.

### Tests ciblés pendant le développement

- `scratch/dev-test.sh <configuration>` lance un test avec l'arbre de
  développement (MFront de `scratch/build-dev2`, en-têtes de `tfel-dev`),
  sans installer TFEL. `DEVTEST_KIND=material-properties` choisit les
  propriétés matériau.
- `scratch/targeted-test.sh <configuration> <motif> <loi>...` ne construit
  que les lois données.
- `reproducers/` : petits programmes autonomes qui vérifient un point
  précis (rotations des blocs de l'opérateur tangent, chaînes de caractères
  dans le code GPU, coût d'un appel synchrone). Voir
  `reproducers/README.md`.

## MGIS

MGIS est construit dans `build/mgis` contre le TFEL du prototype. Ses
tests (`ctest`, 50 tests) comprennent ceux d'origine et :

| Test | Ce qui est vérifié |
|---|---|
| `GenericParallelLoadTest` | chargement, refus croisés des fonctions d'un point et de lot, disposition mémoire des copies des structures de TFEL |
| `GenericParallelIntegrateTest` | lots sur CPU (`stlpar`) comparés au calcul point par point : intégration, plages, échecs partiels, fonctions d'initialisation, post-traitements, rotations 2D et 3D, grandes déformations, vitesse du son, toutes les sortes de variables externes, modèle, politiques d'exécution, emplacements et ressources mémoire |
| `GenericParallelEvaluateTest` | propriétés matériau sur CPU : chargement, valeurs, bornes et politiques, emplacements, arguments invalides, politique d'exécution |
| `GenericParallelCUDATest` | mêmes vérifications sur GPU (`cuda`), données en mémoire managée, ressource mémoire, flux non par défaut, modèle et propriétés matériau |
| `MaterialDataViewTest` | vues : `FieldView` et ses contrats, `getPoints` ; intégration, fonctions d'initialisation et post-traitements sur vues dans quatre rangements (point par point, composante par composante, avec ou sans pas plus grand que nécessaire), identiques bit à bit au gestionnaire, avec un espace de travail, avec les politiques synchrone et asynchrone et par `integrate_debug`, chaque calcul partant d'une copie de l'état avant calcul ; les huit types d'intégration, prédictions comprises ; vitesse du son ; modèle ; énergies stockée et dissipée ; ensembles vides ; vue sans blocs de la tangente ; `make_view` ; opérations de l'hôte ; `update`, `revert`, `updateValues` (avec et sans contexte) ; variables d'une vue d'état, blocs de la tangente, extraction ; refus attendus, pour les deux interfaces |
| `MaterialDataManagerLayoutTest` | gestionnaire rangé composante par composante identique bit à bit au rangement point par point : intégration sur des plages, ensembles vides, pool de threads (erreurs comprises), politiques, `revert` ; initialisation et post-traitement sous toutes leurs formes ; fin de pas des propriétés matériau et variables externes (tableaux externes, locaux, valeurs uniformes, politiques de mise à jour) ; conversions des grandes transformations sur vues et sur gestionnaire, en 3D et en déformation plane, contrainte de Piola-Kirchhoff comparée à sa définition ; vues des variables, vues *strided coalesced* ; `MaterialFunctionManager` ; `save` et `restore` |
| `GenericParallelIntegrateTest`, `GenericParallelEvaluateTest` (ajouts) | rotations et `evaluate` sur vues, identiques bit à bit aux surcharges sur tableaux, pour tous les rangements de la source et de la destination, en place ou non, matrices données des cinq façons, sans politique, synchrone et asynchrone |
| `FieldOperationsCUDATest` | opérations CUDA (copies, transpositions, remplissages), cycle complet sur GPU (tableaux de l'appelant composante par composante en mémoire du GPU seulement, intégration et fin de pas sur le GPU, comparés à l'hôte), intégration, initialisation, post-traitement, `update`, `revert` et extraction directement sur une `MaterialDataView` en mémoire du GPU, rotations et `evaluate` sur vues en mémoire du GPU, calculs asynchrones avec un pool de mémoire managée, masse volumique gardée sur l'hôte refusée |
| `AllocationTest` | aucune allocation sur le tas dans le thread appelant pour l'intégration, `update` et `revert` des gestionnaires (deux rangements) et des vues, les rotations et `evaluate`, pour les deux interfaces |
| `RotationMatrixCopyTest` | copies et déplacements des matrices de rotation, stockage local et externe, rotation après destruction de l'original |

Les tests des liaisons C et Fortran sont construits par la cible `check`
de CMake, pas par `build-tests` : sans elle, `ctest` relance d'anciens
exécutables.

La construction `build/mgis-full` ajoute HDF5 et les liaisons C, Fortran et
Python (96 tests), dont les tests du rangement dans chaque liaison
(`LayoutTest` : rangement du gestionnaire, rangement invalide refusé en C,
blocs de la tangente en Fortran dans les deux rangements, pas des tableaux
numpy et tableaux à pas refusés en Python). Le script `scripts/check-mgis-gpu-headers.sh` compile les
en-têtes GPU de MGIS (opérations de copie et ressources mémoire) avec les
chaînes CUDA, HIP, SYCL, Kokkos et OpenMP déporté du prototype, et les
exécute sur le GPU (`reproducers/mgis-gpu-headers/check.cxx`).

Le test CUDA demande l'option CMake `enable-generic-parallel-cuda-tests`.
Les bibliothèques de lots des tests sont construites par les fonctions de
`mgis/cmake/modules/generic-parallel.cmake`.

## mfem-mgis

mfem-mgis est construit dans `build/mfem-mgis` (ou dans une construction
séparée par lot, `build/mfem-mgis-<lot>`, avec les cibles `BehaviourTest` et
`GenericParallelBehaviour-*`, hors de `all` et de `build-tests`). Sa suite
(234 tests au 12/10/2026, dont 52 variantes qui utilisent le GPU) se lance
de trois façons, qui doivent toutes passer :

- telle quelle : `ctest -j 3 -E
  "cuda|gpu|GPU|AmgX|HypreDevicePoolSize|okkos"` (182 tests, 182/182 le
  12/10 sur `gpu-proto` 3ad387b002), puis les variantes qui utilisent le
  GPU, `ctest -j 2 -R "cuda|gpu|GPU|AmgX|HypreDevicePoolSize|okkos"` (52
  tests) ;
- avec tous les calculs forcés aux noyaux du chemin par lots, par les
  variables `MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch`,
  `MFEM_MGIS_GRADIENTS_COMPUTATION=Kernels`,
  `MFEM_MGIS_INNER_FORCES_COMPUTATION=Kernels` et
  `MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION=Kernels` (182/182 le 12/10, hors variantes
  GPU) ;
- sur le GPU, avec la variable `MFEM_MGIS_DEVICE=cuda` (équivalente à
  l'option `--device cuda` de chaque test), avec ou sans les variables
  précédentes, un test à la fois (`ctest -j 1 -E "cuda|gpu|GPU|debug"`,
  182 tests, 180/182 le 12/10 en 624 s) : des tests simultanés se
  partagent le GPU et échouent. Les
  tests à 4 processus demandent aussi
  `MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE=1` : hypre réserve sinon 4 Gio par
  processus sur le GPU (16 Go). Seuls `DeviceForallTest-cpu` et
  `DeviceOptionTest-none`, qui exigent l'absence de périphérique,
  échouent alors. `FixedPointSolverTest-kernels-cpu`, qui garde
  son option `--device cpu`, et `UnitTestParameters`, qui n'initialise pas
  mfem-mgis, ne passent alors pas par le GPU. `FixedPointSolverTest` a
  échoué une fois sur 26 passages sur le GPU avec les noyaux forcés, sans
  reproduction en 25 passages isolés (certaines accélérations ont un
  nombre d'itérations sensible aux arrondis).

L'option `enable-kokkos` (initialisation de Kokkos par mfem-mgis) se
valide dans deux constructions séparées : `build/mfem-mgis-kokkos-off`
(sans Kokkos : `ctest -R KokkosOptionTest`, 6 tests, option refusée et
variable ignorée avec un avertissement) et `build/mfem-mgis-kokkos` (avec
le Kokkos de l'environnement, construit pour CUDA, installé dans
`install/mfem-mgis-kokkos` : `ctest -R "KokkosOptionTest|-kokkos"`, 7
tests dont 4 sur le GPU). Côté clients, mm-opera-hpc
(`build/mm-opera-hpc-kokkos`) et `cas/` (`build/cas-kokkos`, avec
`MM_OPERA_HPC_POLYCRYSTAL_WITH_KOKKOS=ON`) se construisent contre cette
installation sans chercher Kokkos ; leurs tests `poutre-gpu-kokkos` et
`polycristal-batch-kokkos` (avec sa comparaison) utilisent la loi générée
par le backend kokkos.

La gestion hétérogène de la mémoire (HMM) masque ici les données de l'hôte
lues par erreur sur le GPU : la mutation qui donne au noyau d'un modèle
aux points les pointeurs hôte de ses entrées n'est pas détectée avec HMM.

La factorisation incomplète ILU(k) de hypre 3.2.0 de niveau 1 ou plus
plante sur le GPU (reproducteur `reproducers/hypre-ilu-gpu`) : mfem-mgis y
prend le niveau 0 par défaut et refuse les niveaux supérieurs. La
factorisation ILUT (`HypreILUType` 1) y fonctionne. Le préconditionneur
ParaSails, qui ne calcule que sur l'hôte, est refusé quand hypre calcule
sur le GPU.

Les exemples de mfem-mgis-examples et les cas de mm-opera-hpc se valident
dans `cas/` (construction `build/cas-post`, liée à une installation de
mfem-mgis et aux constructions des exemples et de mm-opera-hpc) : chaque
cas garde la vérification de l'exemple (référence ou solution analytique)
et passe point par point, par noyaux, par lots avec chaque backend de
l'hôte et sur le GPU. Sur le CPU, `ctest -E
"cuda|gpu|GPU|AmgX|HypreDevicePoolSize|okkos"` (les tests Kokkos
utilisent le GPU) ; sur le GPU, `ctest -R cuda`, les tests qui partagent
un fichier de sortie ou le GPU étant sérialisés par des verrous.

## Défauts de compilateurs

Les défauts rencontrés sont contournés de façon ciblée, avec un
reproducteur pour chacun : `results/nvhpc-bugs/README.md`. Le script
`results/nvhpc-bugs/run-nvcxx-reproducers.sh` les rejoue, pour vérifier
s'ils sont corrigés dans une nouvelle version d'un compilateur.
