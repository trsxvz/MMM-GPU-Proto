# Journal

Ce journal raconte le travail jour par jour : ce qui a été fait, les
défauts trouvés, les décisions et leurs raisons. Chaque entrée décrit
l'état au moment où elle a été écrite. Des choix ont pu changer ensuite :
l'état actuel est décrit dans [docs/](docs/vue-ensemble.md).

## Phases

| Repère | Sujet |
|---|---|
| Jalons 1 à 3 | code généré compilable pour GPU, interface `generic-parallel` des lois, premiers tests |
| B1 à B7 | couverture de MFront : grandes déformations, fonctions d'initialisation, post-traitements, rotations, monocristaux, modèles, TDLS |
| C1 à C5 | backends : HIP, `stlpar` avec nvc++, Kokkos, SYCL, OpenMP |
| R1 | corrections de la relecture du 06/10 |
| R2 | lois de test portables entre CPU et GPU |
| R3 | essais d'exécution de toutes les lois de test |
| R4 | validations à `-O3`, environnement spack |
| R5 | clang-tidy et relecture indépendante |
| R6, R7 | données partagées par GPU, paramètres des variables de comportement, rotation de tous les blocs |
| R8 | erreurs sans coût en cas de succès |
| R9 | espace mémoire des bibliothèques |
| R10 | options d'exécution et évaluation asynchrone des propriétés matériau |
| R11 | documentation de TFEL et de MGIS à jour |
| D1 | MGIS : lots de lois, de modèles et de propriétés matériau |
| D2 | mfem-mgis : plan à discuter ([docs/plan-mfem-mgis.md](docs/plan-mfem-mgis.md)) |

## Conventions

- Les chemins cités sont ceux du moment. Plusieurs documents ont été
  déplacés depuis :
  - `DESIGN.md` est devenu [docs/tfel.md](docs/tfel.md) ;
  - `UPSTREAM.md` est devenu [docs/integration-amont.md](docs/integration-amont.md) ;
  - `scratch/review-tfel-2026-10-06.md` est dans
    [docs/historique/](docs/historique/relecture-tfel-2026-10-06.md) ;
  - `scratch/d1-mgis-design.md` est dans
    [docs/historique/](docs/historique/conception-mgis-2026-10-08.md) ;
  - `scratch/d2-mfem-mgis-design.md` est devenu
    [docs/plan-mfem-mgis.md](docs/plan-mfem-mgis.md).
- `tfel/` désigne le dépôt TFEL du prototype, `tfel-dev/` sa copie de
  développement (même dépôt, autre arbre de travail).
- Les mesures de temps sont celles de la machine de développement (GPU RTX
  5080).

## Index des entrées

| Entrée |
|---|
| [2026-10-05, état des lieux de TFEL master](#2026-10-05-état-des-lieux-de-tfel-master) |
| [2026-10-05, sonde J2 (DSL par défaut) dans un kernel CUDA](#2026-10-05-sonde-j2-dsl-par-défaut-dans-un-kernel-cuda) |
| [2026-10-05, jalon 1 terminé : code généré compilable pour le GPU](#2026-10-05-jalon-1-terminé--code-généré-compilable-pour-le-gpu) |
| [2026-10-05, jalon 2 : interface `generic-parallel` des lois](#2026-10-05-jalon-2--interface-generic-parallel-des-lois) |
| [2026-10-05, jalon 3 et balayage des lois de test](#2026-10-05-jalon-3-et-balayage-des-lois-de-test) |
| [2026-10-05, lois de matériau, valeurs propres et pile sur GPU](#2026-10-05-lois-de-matériau-valeurs-propres-et-pile-sur-gpu) |
| [2026-10-05, réorganisation, commits locaux, piège d'environnement](#2026-10-05-réorganisation-commits-locaux-piège-denvironnement) |
| [2026-10-05, HIP sur GPU NVIDIA (C1, commit 0007)](#2026-10-05-hip-sur-gpu-nvidia-c1-commit-0007) |
| [2026-10-05, `stlpar` sur GPU avec nvc++ 25.11 (C2, en cours)](#2026-10-05-stlpar-sur-gpu-avec-nvc-2511-c2-en-cours) |
| [2026-10-05, backends Kokkos et SYCL, factorisation, corrections CMake](#2026-10-05-backends-kokkos-et-sycl-factorisation-corrections-cmake) |
| [2026-10-05, matrice complète des backends, corrections](#2026-10-05-matrice-complète-des-backends-corrections) |
| [2026-10-05 (soir), B5 grandes déformations, B2 initialisations et post-traitements](#2026-10-05-soir-b5-grandes-déformations-b2-initialisations-et-post-traitements) |
| [2026-10-05 (soir, suite), B3 rotations des lois orthotropes](#2026-10-05-soir-suite-b3-rotations-des-lois-orthotropes) |
| [2026-10-05 (nuit), B4 monocristaux, défauts trouvés par nvcc](#2026-10-05-nuit-b4-monocristaux-défauts-trouvés-par-nvcc) |
| [2026-10-05 (nuit, suite), B1 couverture numérique, nouveaux défauts nvcc et clang](#2026-10-05-nuit-suite-b1-couverture-numérique-nouveaux-défauts-nvcc-et-clang) |
| [2026-10-05 (nuit, suite), B6 modèles autonomes](#2026-10-05-nuit-suite-b6-modèles-autonomes) |
| [2026-10-05 (nuit, suite), B7 solveur TDLS sur GPU](#2026-10-05-nuit-suite-b7-solveur-tdls-sur-gpu) |
| [2026-10-05 (nuit, suite), variables de comportement et balayage complet](#2026-10-05-nuit-suite-variables-de-comportement-et-balayage-complet) |
| [2026-10-05 (nuit, suite), brique StandardElastoViscoPlasticity et en-têtes de TFEL](#2026-10-05-nuit-suite-brique-standardelastoviscoplasticity-et-en-têtes-de-tfel) |
| [2026-10-05 (nuit, suite), commits 0016 à 0018](#2026-10-05-nuit-suite-commits-0016-à-0018) |
| [2026-10-05 (soir), C5 backend OpenMP, points d'entrée factorisés](#2026-10-05-soir-c5-backend-openmp-points-dentrée-factorisés) |
| [2026-10-05 (soir, suite), propriétés matériau : backends Kokkos et OpenMP](#2026-10-05-soir-suite-propriétés-matériau--backends-kokkos-et-openmp) |
| [2026-10-05 (soir, suite), matrice c, OpenMP, modèles compilés par nvc++](#2026-10-05-soir-suite-matrice-c-openmp-modèles-compilés-par-nvc) |
| [2026-10-06 (nuit), versions de NVHPC et de CUDA, défauts d'OpenMP déporté](#2026-10-06-nuit-versions-de-nvhpc-et-de-cuda-défauts-dopenmp-déporté) |
| [2026-10-06 (nuit, suite), champs par variable, symboles, rôles](#2026-10-06-nuit-suite-champs-par-variable-symboles-rôles) |
| [2026-10-06 (nuit, suite), pilotes de l'interface generic, OpenMP en mode direct](#2026-10-06-nuit-suite-pilotes-de-linterface-generic-openmp-en-mode-direct) |
| [2026-10-06 (matin), configurations clang CUDA et HIP](#2026-10-06-matin-configurations-clang-cuda-et-hip) |
| [2026-10-06 (matin, suite), balayage clang CUDA de tous les tests de MFront](#2026-10-06-matin-suite-balayage-clang-cuda-de-tous-les-tests-de-mfront) |
| [2026-10-06 (matin, suite), DSL non exécutés sur GPU et deux défauts de nvc++ OpenMP](#2026-10-06-matin-suite-dsl-non-exécutés-sur-gpu-et-deux-défauts-de-nvc-openmp) |
| [2026-10-06 (matin, fin), validation de la branche fusionnée](#2026-10-06-matin-fin-validation-de-la-branche-fusionnée) |
| [2026-10-06 (après-midi), corrections de la relecture (R1)](#2026-10-06-après-midi-corrections-de-la-relecture-r1) |
| [2026-10-06 (soir), lois de test portables (R2)](#2026-10-06-soir-lois-de-test-portables-r2) |
| [2026-10-07 (nuit), essais d'exécution de toutes les lois de test (R3)](#2026-10-07-nuit-essais-dexécution-de-toutes-les-lois-de-test-r3) |
| [2026-10-07 (matin), grandes validations à -O3 (R4, première partie)](#2026-10-07-matin-grandes-validations-à--o3-r4-première-partie) |
| [2026-10-07 (après-midi), passage à l'environnement spack mmm-gpu (R4, en cours)](#2026-10-07-après-midi-passage-à-lenvironnement-spack-mmm-gpu-r4-en-cours) |
| [2026-10-07 (soir), environnement mmm-gpu validé (fin de R4)](#2026-10-07-soir-environnement-mmm-gpu-validé-fin-de-r4) |
| [2026-10-08 (nuit), clang-tidy et relecture indépendante (R5)](#2026-10-08-nuit-clang-tidy-et-relecture-indépendante-r5) |
| [2026-10-08 (nuit, suite), données partagées, symboles, duplication (R6)](#2026-10-08-nuit-suite-données-partagées-symboles-duplication-r6) |
| [2026-10-08 (nuit, fin), paramètres des variables de comportement, rotations (R7)](#2026-10-08-nuit-fin-paramètres-des-variables-de-comportement-rotations-r7) |
| [2026-10-08 (matin), gestion des erreurs homogène (R8)](#2026-10-08-matin-gestion-des-erreurs-homogène-r8) |
| [2026-10-08 (matin, suite), MGIS : chargement et calcul de lots (D1, étapes 1 à 3)](#2026-10-08-matin-suite-mgis--chargement-et-calcul-de-lots-d1-étapes-1-à-3) |
| [2026-10-08 (après-midi), espace mémoire (R9) et relecture de MGIS (D1, étape 4b)](#2026-10-08-après-midi-espace-mémoire-r9-et-relecture-de-mgis-d1-étape-4b) |
| [2026-10-08 (après-midi, suite), MGIS : propriétés matériau et modèles de lots (D1, étape 6)](#2026-10-08-après-midi-suite-mgis--propriétés-matériau-et-modèles-de-lots-d1-étape-6) |
| [2026-10-08 (soir), options d'exécution des propriétés matériau (R10)](#2026-10-08-soir-options-dexécution-des-propriétés-matériau-r10) |
| [2026-10-08 (soir, suite), évaluation asynchrone des propriétés matériau (R10, suite)](#2026-10-08-soir-suite-évaluation-asynchrone-des-propriétés-matériau-r10-suite) |
| [2026-10-08 (soir, fin), files SYCL et blocs de fils des propriétés matériau (R10, fin)](#2026-10-08-soir-fin-files-sycl-et-blocs-de-fils-des-propriétés-matériau-r10-fin) |
| [2026-10-08 (nuit), documentation de TFEL et de MGIS (R11)](#2026-10-08-nuit-documentation-de-tfel-et-de-mgis-r11) |

## Entrées

### 2026-10-05, état des lieux de TFEL master

- Interface `generic-parallel` des propriétés matériau (#966) :
  backends CUDA, HIP, C++ standard, SYCL ; deux prototypes (avec et sans pas
  de mémoire) ; paramètres dans une structure POD passée par valeur au kernel ;
  erreurs dans un tableau ; métadonnées exportées ; tests dans le répertoire
  d'installation (`test-suite/material-properties/generic-parallel`).
- Lois de comportement : le code généré est déjà presque compatible GPU
  (méthodes `TFEL_HOST_DEVICE`, aucune exception). Branche
  `411-introduce-a-derived-generic-interface-for-gpus` : ébauche de 2023,
  abandonnée.

### 2026-10-05, sonde J2 (DSL par défaut) dans un kernel CUDA

Erreurs ou avertissements trouvés, dans l'ordre :
1. `integrate` généré par le DSL par défaut non annoté `TFEL_HOST_DEVICE`.
2. Constructeurs depuis `mfront::gb::BehaviourData` (loi, données, données
   d'intégration) non annotés.
3. Destructeur virtuel de `MechanicalBehaviour` (TFEL) non annoté, donc le
   destructeur généré (`= default`) est hôte seulement.
4. Paramètres lus via le singleton `...ParametersInitializer::get()` (hôte).
5. Vérifications `raise_if(smflag != STANDARDTANGENTOPERATOR, ...)` générées
   dans `integrate` et `computePredictionOperator`.

`tfel::reportContractViolation` ne convient pas pour 5 : il appelle `abort`
sur l'hôte, ce qui changerait le comportement CPU. Décision : macro
`TFEL_DEVICE_CODE` (compilation du code GPU seulement) et branche dans le
code généré, le CPU restant identique.

### 2026-10-05, jalon 1 terminé : code généré compilable pour le GPU

- `integrate`, `computePredictionOperator`, constructeurs et destructeur
  annotés `TFEL_HOST_DEVICE` ; erreurs « impossibles » (`raise`) remplacées
  par `writeDeviceCompatibleRaise` : exception sur l'hôte (inchangé), échec
  retourné sur le GPU.
- Paramètres : le constructeur d'une interface peut lire les paramètres dans
  un argument (`getParametersSourceInBehaviourConstructor`) au lieu des
  singletons hôtes.
- Bornes (`BoundsCheck.hxx` réécrit, même API) : fonctions
  `TFEL_HOST_DEVICE` qui renvoient `bool`. Hôte : exception ou avertissement
  comme avant. GPU : `false` si la politique est stricte. `checkBounds()`
  renvoie `bool`, les bornes physiques vérifiées après intégration font
  `return FAILURE;` sur GPU. Les trois spécialisations `BoundsCheck<1,2,3>`
  remplacées par un seul modèle.
- Le générateur n'écrit qu'une fois les blocs identiques produits par les
  interfaces (constructeurs des données, `exportStateData`...) : `generic` et
  `generic-parallel` peuvent être utilisées ensemble.

### 2026-10-05, jalon 2 : interface `generic-parallel` des lois

- Conception (voir DESIGN.md) : même squelette que l'interface des
  propriétés matériau (fabrique de backends, options dans les
  fichiers de configuration, bibliothèque `GenericParallel<lib>-<backend>`,
  symboles exportés avec l'interface `GenericParallel`). Différence : la
  logique par point est écrite une fois dans des en-têtes `TFEL_HOST_DEVICE`
  (`MFront/GenericParallelBehaviour/Integrate.hxx`, `GPU.hxx`,
  `ParallelSTL.hxx`) ; le code généré par hypothèse fait une quinzaine de
  lignes (structure des paramètres, tailles, appel).
- `Integrate.hxx` : équivalent sans exception de `mfront::gb::integrate`.
  Deux chemins : accès direct si tous les champs sont rangés point par point,
  sinon copie dans des tableaux locaux puis recopie (rangement par
  composante, accès coalescés sur GPU).
- `GPU.hxx` : un thread par point ; bilan (échecs, résultats peu fiables,
  `rdt` minimal) réduit par warp, atomiques seulement si un point du warp a
  quelque chose à signaler ; `cudaMallocAsync` sur le flux de l'appelant.
- Valeurs de retour : 1 succès, 0 `rdt` < 0,99 quelque part, -1 échec en
  au moins un point, -2 arguments invalides, -3 erreur d'exécution.

Pièges trouvés :
1. nvcc ne sait pas générer l'amorce hôte d'un kernel dont un paramètre de
   modèle non-type est une classe (C++20) : tailles passées par un type.
2. **`std::string_view` construit depuis une chaîne C dans du code GPU
   (libstdc++ 14 + nvcc 13) vide le kernel entier, sans erreur ni
   avertissement** (`__builtin_strlen` traité comme inatteignable). Les
   kernels de Norton faisaient 4 registres et ne calculaient rien.
   `BoundsCheck` prend maintenant des `const char*`. Seule la comparaison
   numérique l'a détecté : contrôle à garder (nombre de registres, tests).
3. libstdc++ : `std::for_each(par_unseq, ...)` sur `std::views::iota`
   s'exécute sur **un seul thread** (catégorie C++17 de l'itérateur pas
   « random access »). Le backend `stlpar` des propriétés matériau de master
   est donc séquentiel avec gcc. Backend des lois : itérateur de comptage.
4. MFront plaçait `$(LDFLAGS)` (donc `-ltbb`) avant les objets : TBB ignoré
   avec `--as-needed` (Ubuntu), symboles TBB indéfinis. Corrigé dans
   `MakefileGenerator.cxx` (master a le même défaut).
5. Test CUDA des propriétés matériau de master : un cas UO2 échoue chez nous
   (`errno` non nul après l'appel), sans lien avec nos modifications.
6. Avertissements nvcc 2810 (`nodiscard` sur `operator=` des types de
   TFEL/Math, faux positifs), 20208, 3189 : supprimés dans
   `config-cuda-nvcc.json`, avec `-arch=native` ajouté.

Mesures indicatives (RTX 5080, double précision bridée à 1/64), 100 003
points, appel complet (allocation, copie, synchronisation) : J2 0,3-0,4 ms,
Norton 6-7 ms ; `stlpar` 16 threads : J2 1,6-2,6 ms, Norton 17-21 ms.

### 2026-10-05, jalon 3 et balayage des lois de test

- Tests `tfel-check` : `compare.cxx` compilable par nvcc (mémoire managée)
  ou g++ ; les deux interfaces générées par le même appel de `mfront`
  (vérifie aussi leur coexistence). Composants `tfel-check` ajoutés dans
  `cmake/modules/tfel.cmake` (`mfront::behaviour::generic-parallel_...`).
- Balayage : chaque loi compilée avec le backend CUDA, classée OK si aucun
  appel de fonction hôte depuis le code GPU (avertissement nvcc 20011) et si
  aucun kernel n'a moins de 16 registres. **Un appel hôte n'est pas un
  simple avertissement : nvcc supprime le chemin de code correspondant**
  (kernel vide ou calcul faux sans erreur). 99 lois propres au premier
  passage, 120 après corrections.
- Corrections du générateur : DSL Runge-Kutta (`integrate`, dérivées,
  exceptions de divergence → `return FAILURE` sur GPU, pas de temps
  minimal), `MultipleIsotropicMisesFlows`, DSL isotropes (`computeFlow`,
  `NewtonIntegration`, opérateur tangent), dilatation libre
  (`computeStressFreeStrain/Expansion`), potentiel de Hooke (prédiction
  élastique), accesseurs de blocs du jacobien et de l'opérateur tangent,
  enveloppes des lois internes (`@BehaviourVariable`), drapeau
  `STANDARDTANGENTOPERATOR` qualifié (régression de ma part sur les lois
  générales), profileur désactivé sur GPU, tableaux dynamiques interdits si
  une interface ne les supporte pas (nouvelle méthode
  `areDynamicallyAllocatedVectorsAllowed` des interfaces), paramètres
  « static variables » : membres non statiques en code GPU (nvcc refuse
  l'ODR-use d'un membre `static constexpr`, or les expressions de TFEL/Math
  prennent les scalaires par référence). `FSAlgorithm` annoté.
- Code utilisateur non portable (19 lois) : `throw`, `std::cout`,
  `raise_if`, `static const` initialisé à l'exécution. Sur GPU, utiliser
  `return false;` ou `return FAILURE;`. À documenter.

### 2026-10-05, lois de matériau, valeurs propres et pile sur GPU

- Lois de matériau appelées par les lois (`@MaterialLaw`,
  `@ElasticMaterialProperties`, `@ComputeThermalExpansion`...) :
  l'interface `mfront` des propriétés matériau écrit en plus, dans son
  en-tête, une implémentation `TFEL_HOST_DEVICE inline` dans l'espace de
  noms `mfront::device` (même corps, paramètres constants, pas d'exception
  ni d'`errno`, `NaN` si une borne physique n'est pas respectée).
  `writeMaterialLaws` écrit `using mfront::device::Loi;` en code GPU ; les
  blocs de vérification des bornes de la loi (exceptions, `std::cerr`) ne
  sont compilés que pour l'hôte. Modèles (`@AxialGrowth`...) : classe déjà
  dans l'en-tête, annotée. Vérification : StrainHardeningCreep3 (Lamé
  calculés par `Inconel600_YoungModulus(T)`, grandeurs physiques) identique
  au bit près à l'interface `generic`.
- Valeurs propres : les 8 solveurs de `stensor` (TFEL, FSES ×5, GTE,
  Harari) annotés et vérifiés sur GPU (écart au CPU ≤ 3e-15, erreur de
  reconstruction ≤ 1e-13). **Nouveau piège silencieux : `std::fpclassify`
  en code GPU vide le kernel** (comme `std::string_view` depuis une chaîne
  C). FSES l'utilisait : remplacé par `tfel::math::ieee754::fpclassify`,
  qui fonctionne sur GPU.
- Pile des threads CUDA : 1024 octets par défaut, non augmentée
  automatiquement ; Hosford (puissance réelle, `pow` de libdevice) demande
  2160 octets (estimation statique, l'usage réel est un peu plus grand :
  seuil mesuré entre 2160 et 2220). Le backend GPU lit le besoin du kernel
  (`cudaFuncGetAttributes`, `localSizeBytes`) et augmente la limite si
  nécessaire (+512 octets, arrondi au ko, sous verrou car la limite est
  globale au contexte). Hosford identique à `generic` (≤ 6e-13).
- Tests `tfel-check` étendus : + Hosford, + StrainHardeningCreep3 (loi de
  matériau sur GPU, grandeurs). Les deux configurations passent (CUDA
  1 min 31, CPU 32 s). Le banc compare maintenant aussi les statuts et le
  bilan du premier pas (au lieu d'exiger son succès).
- Documentation : `docs/web/generic-parallel-behaviour-interface.md`.

### 2026-10-05, réorganisation, commits locaux, piège d'environnement

- Dossiers `tfel/` (ex `tfel-gpu-proto`), `mgis/`, `mfem-mgis/`, chacun
  avec `patches/`. 5 commits locaux dans `tfel/` (auteur t), série
  `tfel/patches/0001..0005`. Recompilation complète de TFEL après le
  renommage (chemins absolus dans le cache CMake).
- **Piège d'environnement** : le `LD_LIBRARY_PATH` et le `PATH` du shell
  contiennent la view spack `tfelgpu_local` (TFEL master plus ancien,
  `5696f406f`). `LD_LIBRARY_PATH` passe avant le RUNPATH des binaires :
  le `mfront` de l'arbre de compilation chargeait la libTFELMFront de
  spack. La première passe de la suite de tests de TFEL testait donc
  l'ancien MFront (erreurs BarBehaviour2/4 et TDLS sans rapport avec nos
  modifications). Règle : pour l'arbre de compilation,
  `env -u LD_LIBRARY_PATH PATH=<sans spack>` ; pour l'installation,
  `LD_LIBRARY_PATH=$P/install/tfel/lib:...` et `PATH=$P/install/tfel/bin:...`
  en tête (ce que font les scripts du proto).
- Vu au passage : `~/Documents/GitHub/mgis-gpu` (prototype
  de 01-02/2026) : noyaux CUDA écrits à la main sur les vues
  `StridedCoalescedMemoryAccess` de MGIS. Le format des champs à deux pas
  de l'interface correspond à ces vues (utile pour D1).
- Suite de tests de TFEL (commande validée, ~2 min 15 après
  compilation) : dans `build/tfel`,
  `CLEANPATH=$(echo "$PATH" | tr ':' '\n' | grep -v spack | paste -sd:)`,
  `touch mfront/src/mfront` (si MFront a changé),
  `env -u LD_LIBRARY_PATH PATH="$CLEANPATH" make -k -j16 build-tests`, puis
  `env -u LD_LIBRARY_PATH PATH="$CLEANPATH" ctest -j16`. Résultat :
  6963/6963.
- Question ouverte (interface `generic`, pas modifiée) : le post-traitement
  ne remplace que les forces thermodynamiques de l'état initial par celles
  de fin de pas, alors que le commentaire parle aussi des variables
  internes ; et il passe `*(nullptr)` comme état initial quand il n'est
  pas utilisé (comportement indéfini). La version GPU passera une
  référence valide et gardera la sémantique de `generic`.
- Erreur corrigée : TFEL a un dossier `patches/` suivi (patchs homebrew) ;
  ma série l'avait remplacé. Restauré ; la série du proto est dans
  `proto-patches/` (dans les trois bibliothèques).

### 2026-10-05, HIP sur GPU NVIDIA (C1, commit 0007)

- `hipcc` sur plateforme NVIDIA appelle `nvcc` et ne définit ni `__HIP__`
  ni `__HIPCC__` : `GPU.hxx` prenait silencieusement le chemin CUDA. Le
  backend `hip` inclut maintenant `hip/hip_runtime.h` avant `GPU.hxx`
  (comme le backend HIP des propriétés matériau), ce qui définit
  `__HIPCC__` ; sélection dans `GPU.hxx` : `__HIPCC__ || __HIP__`.
- Écarts CUDA/HIP corrigés : attribut de taille de warp
  (`cudaDevAttrWarpSize` / `hipDeviceAttributeWarpSize`), surcharge non
  typée seule pour `hipFuncGetAttributes`.
- `compare.cxx` : runtime HIP si `__HIP__` ou `MFRONT_GPB_COMPARE_USE_HIP`.
- Configuration `hip-hipcc` existante : drapeaux `nvcc` déplacés dans
  `TFEL_HIPCC_FLAGS` (plateforme nvidia seulement ; `-fPIC` sur amd) ;
  `-isystem @TFEL_HIP_INCLUDE_DIR@` retiré (HIP installé dans `/usr` :
  `-isystem /usr/include` casse les `#include_next` de libstdc++, et hipcc
  connaît déjà son répertoire). Reste le même défaut dans la configuration
  `hip-clang-nvidia-platform` (non activée ici).
- Lancement des tests installés : copier
  `install/tfel/share/tfel/test-suite/{behaviours,material-properties}/generic-parallel`
  dans `scratch/tc`, puis `tfel-check` dans le sous-dossier, avec
  `PATH=install/tfel/bin:...` et, pour HIP,
  `HIP_PLATFORM=nvidia HIP_COMPILER=nvcc CUDA_PATH=<SDK>/cuda/12.9`.

### 2026-10-05, `stlpar` sur GPU avec nvc++ 25.11 (C2, en cours)

- `nvc++ -stdpar=gpu` compile hôte et GPU en une passe, sans `__CUDA_ARCH__` :
  `TFEL_DEVICE_CODE` est maintenant aussi défini si `__NVCOMPILER_STDPAR_GPU`
  (toute fonction peut être déportée, donc le code « device » est choisi
  dans toute l'unité de traduction ; il est valide sur l'hôte).
- **Bogue nvc++ 25.11 n° 1 (prouvé, reproducteurs dans
  `results/nvhpc-bugs/`)** : dans le code GPU à `-O2`/`-O3`, une boucle
  `for (int k = 0; k != N; ++k)` (compteur **signé**, `!=`, borne
  constante) fait N+1 itérations (`int`, `long`, `short` ; `unsigned`,
  `size_t`, `unsigned short` et `<` corrects ; `-O1` correct ; aussi en
  CUDA via `nvc++ -cuda` ; nvcc 13 correct). Le PTX contient déjà la
  constante fausse (22,2 au lieu de 21,6) : bogue de l'optimiseur de nvc++,
  pas de ptxas. Variante : noyau entièrement vide (`empty-kernel-stdpar.cxx`).
  TFEL : seules mes boucles de `BoundsCheck` (surcharges stensor) et
  `fses::fill` avaient ce motif → compteurs `unsigned short` (cohérent avec
  le reste de TFEL). Les générateurs émettent des compteurs non signés. Le
  code utilisateur des `.mfront` reste exposé (à documenter).
- **Bogue nvc++ 25.11 n° 2** : en mode mémoire `unified` (défaut sur une
  machine HMM comme celle-ci), écriture `STG` (globale) vers une adresse de
  pile (`R212 = cvta.local(SP)` dans le SASS) dans le constructeur de la
  loi → `cudaErrorInvalidAddressSpace`, quel que soit `-O`. Pas réduit à un
  reproducteur minimal (objets auto-référents, constructeur non inliné :
  corrects isolément). **Contournement propre : `-gpu=mem:managed`**, qui
  est de toute façon le réglage portable de stdpar (ne dépend pas de la
  présence de HMM ; sur une machine HMM le GPU accède quand même à toute la
  mémoire).
- Pile : Hosford demande 6 à 8 Ko avec nvc++ (2,2 Ko avec nvcc) ;
  `NV_ACC_CUDA_STACKSIZE` (ou `NVCOMPILER_ACC_CUDA_STACKSIZE`) est pris en
  compte par stdpar.
- Mode managed, `-O2` : 5/6 lois OK d'emblée, Hosford OK avec
  `NV_ACC_CUDA_STACKSIZE=8192`. `-O1` : Hosford encore en échec à 8 Ko (à
  voir).

### 2026-10-05, backends Kokkos et SYCL, factorisation, corrections CMake

- Factorisation : `PointIntegrator` (objet fonction trivialement copiable
  intégrant la loi en un point et renvoyant les statistiques),
  `CombineStatistics` et `dispatch` (choix du chemin direct ou
  rassemblement/dispersion selon le rangement des champs) dans
  `Integrate.hxx`. `stlpar`, Kokkos et SYCL les utilisent ; plus de lambdas
  dupliquées, pas de lambdas étendues requises (nvcc, Kokkos). CUDA/HIP
  gardent leur noyau (réductions par warp).
- Kokkos (`Kokkos.hxx`) : réduction par foncteur (`value_type`, `init`,
  `join`), pas de `KOKKOS_LAMBDA` ; sur CUDA/HIP, noyau sonde
  (`stackSizeProbe`, jamais lancé) pour dimensionner la pile, Kokkos ne le
  faisant pas. Le code appelant initialise Kokkos ; `queue` peut pointer
  sur une instance de l'espace d'exécution.
- Configuration `kokkos` de TFEL : compilateur et options extraits des
  cibles importées de Kokkos (`find_package(Kokkos COMPONENTS
  separable_compilation)`, sinon Kokkos redirige toute la compilation de
  TFEL vers son compilateur) ; `-ccbin` = compilateur de TFEL ;
  `-cudart shared` (Kokkos utilise le runtime CUDA partagé : un seul
  runtime dans le processus). Pièges de `nvcc_wrapper` : `-I dir` avec
  espace pris pour un fichier source (utiliser `-Idir`) ; double `-std`
  → avertissement sur la sortie standard qui corrompt les fichiers `.d`
  (tfel-config ajoute déjà `-std=c++20`).
- Correction CMake de TFEL : la substitution d'un compilateur
  (`@TFEL_NVHPC_COMPILER@`, `@TFEL_HIPCC_COMPILER@`, ...) n'était
  enregistrée que si le compilateur était cherché par CMake (absente s'il
  est donné par l'utilisateur ou déjà en cache) ; les substitutions sont en
  cache : changer une valeur impose `cmake -U "tfel-substitutions*" .`.
- `compare.cxx` : couches d'exécution Kokkos (`ScopeGuard`,
  `kokkos_malloc<SharedSpace>`) et SYCL (file globale, USM partagée, file
  transmise par `BatchData::queue`).
- `scripts/env.sh` : environnement complet (TFEL, CUDA 12.9, LLVM 20,
  AdaptiveCpp, nvc++, Kokkos) pour reproduire tous les tests.
- SYCL / AdaptiveCpp : le compilateur générique (SSCP) est mono-passe →
  `TFEL_DEVICE_CODE` défini si `__ACPP_ENABLE_LLVM_SSCP_TARGET__`, comme
  pour nvc++. Sur GPU : `ptxas fatal: Unresolved extern variable
  _ZTVN10__cxxabiv117__class_type_infoE` : la loi est polymorphe (interface
  virtuelle pure de `MechanicalBehaviour`), sa vtable référence la RTTI ;
  AdaptiveCpp exige la RTTI, SYCL n'admet pas les fonctions virtuelles dans
  le code device. **Décision** : `MechanicalBehaviour` garde les types et
  constantes mais perd ses fonctions virtuelles (personne n'utilise les
  lois via cette base dans TFEL ; les lois générées sont `final`, même les
  « variables de comportement » dont les enveloppes n'ajoutent que l'amitié) ;
  l'interface est décrite par le concept C++20
  `tfel::material::MechanicalBehaviourConcept`, vérifié par `static_assert`
  dans les pilotes `generic` et `generic-parallel` ; les générateurs
  n'émettent plus `final`/`override` sur ces méthodes ni sur le
  destructeur. Sur CPU (OpenMP) : plantage dans le JIT d'AdaptiveCpp
  (`SubCfgFormationPass`), à revoir après ces changements.

### 2026-10-05, matrice complète des backends, corrections

- Une seule matrice de tests : `scripts/run-generic-parallel-tests.sh`
  (copie des tests installés dans `scratch/tc`, `tfel-check` par
  configuration, environnement `scripts/env.sh`).
- Piège : MFront lit la variable d'environnement `HIP_COMPILER` pour
  choisir le compilateur HIP ; `HIP_COMPILER=nvcc` (destiné à hipcc) lui
  faisait appeler le nvcc 13 du PATH. Avec `HIP_PLATFORM=nvidia`, hipcc
  utilise déjà nvcc : `HIP_COMPILER` retiré de `env.sh` et `build-tfel.sh`.
- CMake de TFEL : `HIP_PLATFORM` mis en cache (une reconfiguration
  automatique sans la variable d'environnement échouait).
- Pile des threads GPU : le pilote CUDA agrandit la pile au besoin statique
  du noyau, mais ce besoin ignore les routines internes du compilateur
  (`__cuda_sm20_div_rn_f64_full`, division double) → débordement de pile
  (vu avec compute-sanitizer pour Hosford en SYCL/AdaptiveCpp). Nouvel
  en-tête partagé `MFront/GenericParallel/DeviceStackSize.hxx` :
  `setMinimalCUDAStackSize` agit via le runtime CUDA déjà chargé dans le
  processus (`dlopen(..., RTLD_NOLOAD)`, aucune dépendance CUDA à la
  compilation, sans effet si CUDA n'est pas utilisé). Option
  `device_stack_size` des backends `stlpar` et `sycl`, pour les lois ET les
  propriétés matériau (les fichiers de configuration sont partagés).
- Rapports de violation de contrat : désactivés dans tout code device
  (`TFEL_DEVICE_CODE`), pas seulement CUDA/HIP ou après l'inclusion de
  SYCL (AdaptiveCpp sur CPU : `reportContractViolation` non défini).
- **Récursion dans TFEL** : `StensorComputeEigenVectors<3>::computeEigenVector`
  (valeur propre supposée simple) retombait, dans le cas dégénéré, sur la
  décomposition analytique complète qui rappelle `computeEigenVector` :
  cycle dans le graphe d'appel (détecté par une analyse de l'IR LLVM du
  noyau), qui fait boucler à l'infini une passe d'AdaptiveCpp sur CPU
  (`fillTransitiveSplitterCallers`). Repli remplacé par l'algorithme de
  Jacobi (`fses::syevj3`), sans récursion. Seul le cas dégénéré change
  (vecteurs propres non uniques).
- Générateurs : deux derniers `override` (`IsotropicMisesPlasticFlow`,
  `IsotropicStrainHardeningMisesCreep`) retirés ; la compilation des tests
  de TFEL les a trouvés.
- Commits locaux : `TFEL: behaviours without virtual functions`,
  `TFEL: device code selection for single pass compilers` ; les backends
  seront commités après la suite de tests de TFEL.
- **Défauts trouvés dans le code `generic` (à corriger dans la
  refonte B5)** :
  1. `logarithmic_strain::integrate`, opérateur `DTAU_DDF` en fin de pas :
     utilise `lgh0` (début de pas) avec `T1`, alors que la branche
     `DSIG_DF` utilise `lgh1` ;
  2. `executePostProcessing` : sauvegarde/restaure
     `d.s0.internal_state_variables` sans jamais y mettre les valeurs de fin
     de pas (le commentaire dit le contraire) : les post-traitements voient
     les variables internes initiales avec les contraintes finales ;
  3. `executePostProcessing` sans état initial : `*(nullptr)` passé par
     référence (comportement indéfini).
  4. (corrigés dans le commit 494929641, refonte des pilotes `generic` en
     cœurs appelables sur GPU) : statut `if (r)` des pilotes grandes
     déformations, fonctions d'initialisation GL/Hencky aux arguments
     inversés (erreur de compilation), décodage de la prédiction qui
     ignorait le décalage de la vitesse du son, pointeurs non restaurés si
     une fonction d'initialisation lève une exception.

### 2026-10-05 (soir), B5 grandes déformations, B2 initialisations et post-traitements

- **B5** (non commité au moment de l'écriture) : stratégies d'intégration
  `FiniteStrainIntegrator`, `GreenLagrangeStrainIntegrator`,
  `LogarithmicStrainIntegrator` (`FiniteStrainIntegrators.hxx`), tailles
  par point dépendant des options (forces : T pour PK1, S sinon ;
  opérateur : S·S pour DS_DEGL, T·T pour DPK1_DF, S·T sinon), traits
  génériques générés (`writeGenericBehaviourTraits`), outils grandes
  déformations annotés pour le device (`FiniteStrainBehaviourTangentOperator`,
  `LogarithmicStrainHandler` en en-tête seul, LU sans exception sur device,
  `TinyMatrixInvert` renvoie un booléen).
- **B2** : `_InitializeFunction_<nom>` et `_PostProcessing_<nom>` par lot,
  mêmes statuts que l'intégration (1 en cas de succès, alors que
  l'interface `generic` renvoie 0). Primitive commune des backends
  `apply(status, n, queue, rdt0, f, ...)` : intégration, fonctions
  d'initialisation et post-traitements = vérification + choix du noyau +
  `apply`. `PointView` (données `generic` d'un point, copies locales si les
  champs ne sont pas rangés point par point). Entrées nulles acceptées pour
  une fonction d'initialisation sans entrée (comme MGIS).
- Harnais de comparaison (`compare.cxx`) : options nommées
  `--finite-strain=SM,TO`, `--initialize-function=n`, `--post-processing=n`,
  états aléatoires pour initialisations et post-traitements, contrôle que la
  référence modifie bien quelque chose (sinon un test vide passerait),
  tolérance relative 1e-10 documentée (compilateurs différents pour les
  deux bibliothèques : contraction FMA différente, amplifiée par les Newton
  locaux ; SYCL/OpenMP donnait 1,1e-12 sur Hencky contraintes planes). Les
  comparaisons sont listées une seule fois dans `compare-all.sh` (motif de
  sélection optionnel), appelé par les 8 `test.check`.
- Défauts trouvés par ces tests et corrigés :
  1. rangement par composante : la vue locale ne recopiait pas l'état de
     fin de pas, lu par les post-traitements et seulement partiellement
     réécrit par certaines fonctions d'initialisation → `gatherState()`
     (test de mutation : le retirer fait échouer le nouveau test) ;
  2. **nvc++ 25.11** (`-stdpar=gpu`) : trois défauts de génération de code,
     reproducteurs dans `results/nvhpc-bugs/` :
     - visite de `std::variant` (opérateur tangent en grandes
       déformations) : « Unexpected branch type » puis alias LLVM
       invalides pour l'union interne de libstdc++ → `std::variant` exclu du
       code device : nouvelle classe `tfel::utilities::TrivialGenTypeBase`
       (union récursive, même interface que `GenTypeBase`, trivialement
       copiable, test unitaire `TrivialGenTypeBaseTest`), utilisée par
       `FiniteStrainBehaviourTangentOperator` ; constructeur par défaut de
       l'union vide (sinon nvc++ émet encore un alias) ;
     - initialiseur par défaut `= {}` d'un membre tableau :
       « parse integer constant must have integer type »
       (`aggregate-zero-init.cxx`) → membre initialisé dans le constructeur ;
     - **appel par pointeur sur fonction membre ignoré silencieusement sur
       GPU** (`member-function-pointer.cxx`, résultat faux à -O0..-O3, sans
       diagnostic) → les pilotes reçoivent un objet fonction ; MFront génère
       un foncteur nommé par fonction d'initialisation et post-traitement ;
       l'API `generic` à pointeur membre reste (enveloppe hôte, `MethodCaller`).
  3. `GenericFixedSizeArray` (base de tous les tenseurs de taille fixe) :
     copie et déplacement écrits à la main → types non trivialement
     copiables ; passés en `= default` (même sémantique sur les éléments
     logiques), ce qui rend `stensor`, `st2tost2`, etc. trivialement
     copiables (memcpy, `variant` sans visite, arguments SYCL).
  4. nvcc ne fait qu'avertir (#20011) quand du code device appelle une
     fonction hôte et produit une bibliothèque fausse → `-diag-error 20011`
     dans les configurations utilisant nvcc (`cuda-nvcc`, `hip-hipcc` sur
     NVIDIA).
- Constat amont non corrigé (sémantique) : en grandes déformations, la
  contrainte modifiée par une fonction d'initialisation est écrite dans un
  tampon temporaire et perdue (sauf mesure de Cauchy en grandes
  déformations standard) ; les lois de test FS initialisent donc des
  variables internes.
- Défauts amont dans le générateur des rotations `generic` (B3, à corriger
  dans la version par lots) : variante « tableau » d'un tenseur d'ordre 2
  qui n'est pas la première variable écrite sans le décalage du point ;
  scalaires non recopiés dans `dest` ; vecteurs non pris en charge ;
  `_dtau_ddF` utilise les arguments de `_dsig_dF` (sans effet).
- Balayage CUDA à refaire après installation (les 33 « UNSUPPORTED »
  étaient des lois grandes déformations).

### 2026-10-05 (soir, suite), B3 rotations des lois orthotropes

- Module partagé `MFront/BehaviourRotationFunctions.{hxx,cxx}` :
  description des fonctions de rotation d'une loi orthotrope (nom, objets
  tournés avec type et décalage, taille par point, sens de rotation,
  éventuelle raison de non-prise en charge) et écriture du code de rotation
  d'un point, utilisable sur GPU :
  `map<T>(dest + o) = change_basis(T{src + o}, r);`, scalaires recopiés.
- Interface `generic` : ses ~900 lignes de génération des rotations
  remplacées par ~60 lignes appuyées sur ce module (fonctions et noms
  inchangés). Défauts amont corrigés au passage : décalage du point oublié
  pour les tenseurs d'ordre 2 qui ne sont pas la première variable (variante
  « tableau »), scalaires non recopiés dans `dest`, vecteurs et blocs
  vecteur/vecteur non pris en charge (`std::exit(-1)` à l'exécution !) :
  les intégrateurs de thermique orthotrope de mfem-mgis en dépendent.
  Les blocs restant non pris en charge (tableaux, vecteur/tenseur) lèvent
  une exception (interface `generic`) ou renvoient le statut -2 (lots).
- TFEL/Math : `change_basis` pour `tvector<N>` et `tmatrix<N, N>`, même
  convention que les tenseurs (test `ChangeBasisTest` : cohérence avec
  `change_basis` des tenseurs symétriques en 2D et 3D).
- Interface `generic-parallel` : `mfront_gpb_RotationData` (destination,
  source, matrices par point ou uniformes si le pas de point est nul,
  file), fonctions `<loi>_<H>_rotate...` par lot, foncteur de rotation
  d'un point généré par hypothèse, `RotationExecutor`, `rotate` dans les 4
  backends ; SYCL : gestion des copies extraite dans `DeviceCopies`,
  réutilisée par `DeviceBatchData` et par les rotations.
- Tests : loi `OrthotropicStationaryHeatTransfer` (vecteurs, matrice de
  conductivité), copies de `OrthotropicElastic` et
  `OrthotropicSaintVenantKirchhoffElasticity` ; harnais `--rotation=n`
  (matrices aléatoires par point ou uniforme, deux rangements, référence =
  fonction `generic` sur place) ; 15 fonctions de rotation et 5 intégrations
  de lois orthotropes ajoutées à `compare-all.sh`.
- **Défaut trouvé par ces tests** : avec `@RequireStiffnessTensor` (ou la
  dilatation thermique), le code `generic` lit les propriétés élastiques en
  tête des propriétés matériau (MGIS les ajoute), mais les tailles par point
  du lot les excluaient → lecture de valeurs non initialisées en rangement
  par composante. Corrigé (taille totale lue par la loi).
- stlpar/g++ : 208/208 comparaisons.
- Matrice : `MAKEFLAGS=-j6` dans le script (MFront appelle `make`), option
  `GPB_OPTIMIZATION_LEVEL` (-O1, -O3 : configurations copiées dans chaque
  dossier de test, où MFront les cherche en premier) et
  `GPB_TESTS_DIRECTORY`.

### 2026-10-05 (nuit), B4 monocristaux, défauts trouvés par nvcc

- **Méthodes générées non annotées** : `execute<Nom>InitializeFunction` et
  `execute<Nom>PostProcessing` n'étaient pas `TFEL_HOST_DEVICE` (invisible
  en stlpar/g++, refusé par nvcc grâce à `-diag-error 20011,20014`) →
  annotées dans `GenericBehaviourInterface`. La matrice lancée avant ce
  correctif (CUDA, HIP, Kokkos en échec en 24 s) a été arrêtée.
- **Systèmes de glissement générés** (`@CrystalStructure`,
  `@SlidingSystem`) : singleton à `static` local initialisé dynamiquement,
  interdit sur GPU. Constructeur `constexpr` (toutes les valeurs sont des
  littéraux), méthodes `TFEL_HOST_DEVICE`, `areCoplanar` sans
  `std::min/max`, erreur d'`offset` sans exception sur device, accès par
  `tfel::utilities::getConstantObject<T>()` (nouvel en-tête
  `TFEL/Utilities/ConstantObject.hxx`) :
  - nvcc remplit de **zéros, sans diagnostic**, la copie device d'un objet
    `static constexpr` local dont le constructeur est appelé (reproducteur
    `results/nvhpc-bugs/nvcc-static-constexpr-object.cu`) → variable
    `__device__ constexpr` d'espace de noms en CUDA/HIP ;
  - SYCL (AdaptiveCpp, GPU et CPU) : variable `inline constexpr` lue
    correctement ;
  - nvc++ `-stdpar=gpu` refuse explicitement variables globales et
    statiques locales en code GPU → copie par valeur (à lier à une
    référence constante, comme le code généré le fait déjà).
- Classes historiques de TFEL : `MetallicCFCSlidingSystems` (même
  traitement) et `MetallicCFCGenericSlidingSystemsInteractionMatrix`
  (matrice renvoyée par valeur, calculée à chaque appel : supprime le cache
  du premier appel, piège documenté qui ignorait les coefficients des
  appels suivants ; impossible sur GPU).
- TFEL/Math : expressions `StensorProductExpr{1,2,3}D`,
  `BuildFromRotationMatrix`, accès matriciel des tenseurs, produit scalaire
  de vecteurs, dérivées de tableaux de taille fixe : annotés
  `TFEL_HOST_DEVICE`.
- nvcc : avertissement #20014 (même défaut que #20011, sans nom) aussi
  passé en erreur.
- Harnais : `MFRONT_GPB_COMPARE_NUMBER_OF_POINTS` (débogage),
  `--initial-internal-state-variables=v` ou liste complète, échec si tous
  les points échouent au premier pas (test vide).
- Tests ajoutés : `SingleCrystal_DD_FCC` (systèmes générés, densité de
  dislocations initiale) et `FiniteStrainSingleCrystal` (classes TFEL,
  `Fe` initial = identité). CUDA : écarts 1e-15 et 5e-12 ; le monocristal
  intègre 100 003 points en 65 ms sur GPU contre 640 ms pour la référence
  séquentielle.
- stlpar/g++ : 220/220 comparaisons. Réinstallation, régression et matrice
  relancées.

### 2026-10-05 (nuit, suite), B1 couverture numérique, nouveaux défauts nvcc et clang

- Harnais : énergies stockée et dissipée comparées à chaque intégration
  (`--energies=true` exige qu'elles soient non nulles), vitesse du son
  (`--speed-of-sound=<masse volumique>` : décalage +100 de `K[0]`, masse
  volumique uniforme), type d'opérateur tangent (`--tangent-operator`),
  perturbation aléatoire des variables internes initiales et politique de
  bornes dans les deux bibliothèques (`--out-of-bounds-policy=<loi>,Strict`)
  pour des lots à statuts mêlés.
- Lois ajoutées : `ImplicitNortonWithEnergies` (`@InternalEnergy`,
  `@DissipatedEnergy`, `@SpeedOfSound`, bornes standard), `Norton_rk42`
  (Runge-Kutta), `Tvergaard` (zone cohésive). stlpar : tout concorde ; avec
  `Strict`, 97 % d'échecs contre 50 % (bornes physiques seules), statuts
  identiques.
- **clang CUDA** (configuration `cuda-clang`) ne compilait aucune loi :
  `TFEL_HOST_DEVICE [[nodiscard]]` est refusé (un attribut C++11 doit
  précéder `__host__ __device__`) → 285 corrections dans 41 en-têtes et
  dans les générateurs de MFront ; définition de `fses::is_negligible` sans
  l'annotation de sa déclaration (refusé par clang).
- **nvcc 13.0, optimiseur NVVM** : la contrainte écrite par la fonction
  d'initialisation de SVK (grandes déformations, mesure de Cauchy) était
  perdue en rangement par composante (CUDA et HIP, qui utilise nvcc).
  Juste avec clang CUDA (même code, même GPU), avec `-G` et avec
  `-Xcicc -O0` ; le résultat change si l'on ajoute des `printf`.
  Reproducteur réduit (non minimal) dans `results/nvhpc-bugs/`. Contournement
  retenu : l'état de fin de pas est copié dans le constructeur de
  `PointView` (comme les autres entrées) au lieu d'une méthode
  `gatherState()` appelée après la construction ; les sorties de la vue sont
  relues par les pointeurs des données de la loi (`PointView` possède
  désormais `rdt`, accesseur `getTimeStepScalingFactor()`).

### 2026-10-05 (nuit, suite), B6 modèles autonomes

- Nouvelle interface de modèles `generic-parallel`
  (`GenericParallelModelInterface`, enregistrée à côté de `generic`) :
  mêmes backends et mêmes options que les lois, bibliothèque
  `GenericParallel<bibliothèque|matériau|Model>-<backend>`, en-tête
  `MFront/GenericParallelModel/<nom>-generic-parallel-<backend>.hxx`, une
  fonction par hypothèse (même nom que l'interface `generic`), symboles de
  l'interface `generic`, `<nom>_setParameter`.
- Un modèle est traité comme une loi sans gradients (sorties = variables
  internes, entrées = variables externes, propriétés matériau constantes =
  propriétés matériau) : `mfront::gpb::ModelIntegrator` (stratégie
  d'intégration, comme `StandardIntegrator`) construit le modèle à partir
  des données du point et des paramètres puis appelle `execute`. Aucun
  changement dans les noyaux des backends.
- Générateur de la classe du modèle factorisé dans `GenericModelInterface`
  (`writeGenericModelClass`, options `device` et `parameters_type`) :
  sortie hôte identique à l'ancienne sur 7 modèles des tests de TFEL ; en
  mode device, constructeur `TFEL_HOST_DEVICE` qui lit les paramètres dans
  une structure passée par valeur, fonctions `bool execute_<f>` (bornes :
  `return false;`), `execute()` qui les enchaîne, lois de matériau de
  l'espace `mfront::device` (`writeMaterialLaws`).
- Défaut amont corrigé : un modèle avec des propriétés matériau et des
  paramètres traités comme variables statiques
  (`parameters_as_static_variables`) ne compilait pas avec l'interface
  `generic` (les paramètres statiques étaient initialisés dans le
  constructeur) ; vérifié avec le mfront installé.
- Sémantique : échec → `rdt = 0.1` comme `generic` ; exceptions hôte
  (bornes physiques, lois de matériau) rattrapées point par point
  (`callAndReportExceptions`, sinon `std::terminate` dans l'algorithme
  parallèle) ; pas d'opérateur de prédiction : en prédiction, le modèle
  n'est pas exécuté et l'état n'est pas mis à jour (comme une loi). Le
  harnais saute la prédiction quand l'opérateur tangent est vide et
  n'écrit plus dans `K` vide.
- Tests : `GenericParallelModelTest` (paramètre, propriété matériau, loi
  de matériau, bornes, borne physique pour les échecs, deux fonctions
  enchaînées) et `GenericParallelModelWithQuantitiesTest` (`@UseQt`),
  ajoutés aux 8 configurations et à `compare-all.sh`. stlpar et CUDA :
  écarts nuls, statuts identiques (50 172 échecs sur 100 000 avec
  perturbation). Modèles de TFEL : 6/7 compilent pour CUDA ; le 7e
  (`SiC_IrradiationSwellingModel_GoFaster`) lance une exception dans son
  code (non supporté sur GPU, documenté).
- Second défaut amont corrigé : une variable statique d'un type avec
  unités (`@StaticVariable temperature Tref = 300;`) d'un modèle `@UseQt`
  ne compilait pas (`static constexpr temperature Tref = 300;`) ; elle est
  désormais écrite `static constexpr auto Tref = temperature{300};`
  (testé par `GenericParallelModelWithQuantitiesTest`, valeur correcte sur
  GPU).

### 2026-10-05 (nuit, suite), B7 solveur TDLS sur GPU

- `@LinearSystemSolver TDLS` (bibliothèque TDLS, construite dans
  `~/Documents/GitHub/tdls/build-install`, TFEL configuré avec
  `enable-tdls`) fonctionne tel quel avec `generic-parallel` : la
  bibliothèque est appelable sur GPU et le code généré ne passe que des
  types TFEL. `ImplicitNortonTDLS` (Norton implicite, opérateur tangent
  cohérent par l'inverse partiel du jacobien, donc `factorize`,
  `substitute_inplace` et `solve_inplace` de TDLS) : CUDA, écarts 1e-14 à
  2e-13 avec l'interface `generic` ; stlpar g++ OK.
- Tests : composant `tfel-check` `mfront::linear_system_solver::tdls`
  déclaré par CMake si `enable-tdls` ; sous-dossier `tdls/` dans chacune
  des 8 configurations (`@Requires` configuration + TDLS), script
  `compare-tdls.sh`. La fonction `run` et les variables communes sont
  passées dans `compare-functions.sh`, lu par `compare-all.sh` et
  `compare-tdls.sh`.
- `scripts/run-generic-parallel-tests.sh` : niveau d'optimisation appliqué
  aussi aux sous-tests ; résultat FAILED si un des tests échoue (avant, un
  seul SUCCESS suffisait) ; nombre de tests réussis affiché. Nouvelle
  version gardée dans `scratch/run-generic-parallel-tests.sh.new` le temps
  que la matrice en cours finisse (bash lit le script au fil de l'eau).
- Point ouvert (TDLS, hors GPU) : g++ 13 et 14 émettent un avertissement
  `-Wattributes` (« declared with greater visibility than the type of its
  field ») sur des lambdas de `tfel::math` en compilant
  `ImplicitNortonTDLS` avec l'interface `generic`, même sans les options
  de `tfel-config` ; absent sans TDLS et dans les tests TDLS de TFEL
  (brique `StandardElasticity`). Inoffensif ; cause non isolée (piste :
  paramètres de patron de type classe de TDLS et lambdas dans des
  initialiseurs de membres statiques).

### 2026-10-05 (nuit, suite), variables de comportement et balayage complet

- Constat : le balayage de compilation (`scripts/compile-sweep.sh`) ne
  couvrait que les 226 fichiers à la racine de `mfront/tests/behaviours`.
  Restaient hors champ : `StandardElastoViscoPlasticity` (140 fichiers, la
  brique des lois de plasticité réelles), `StandardElasticity` (21),
  `IsotropicDSLs` (10), `AuxiliaryModels` (7), `DDIF2` (6),
  `BehaviourVariables` (4), `FiniteStrainSingleCrystal` (2), `tdls` (12),
  `unit-tests` (117, dont 96 doivent échouer) et les modèles (7). Le script
  parcourt désormais tout `behaviours` et `models` ; un échec de MFront est
  classé `INVALID` si l'interface `generic` échoue aussi (tests unitaires
  qui doivent échouer).
- Variables de comportement (`@BehaviourVariable`,
  `@BehaviourVariableFactory`, `@Model` dans une loi) : appels hôte depuis
  le GPU (`ExternalModelCaller2/4` dans l'ancien balayage). Corrigé dans
  `BehaviourCodeGeneratorBase` : constructeur et `make()` de la fabrique,
  constructeurs par défaut de la loi interne et de ses classes
  `BehaviourData`/`IntegrationData` en `TFEL_HOST_DEVICE`. Les paramètres
  d'une loi interne (non statiques) étaient lus dans son singleton hôte
  (`ParametersInitializer::get()`, non exposé par la bibliothèque de la loi
  principale) : en code device, ils prennent leurs valeurs par défaut
  (méthode `writeBehaviourParameterInitialisationToDefaultValues`, extraite
  du code existant des paramètres traités comme variables statiques).
  CUDA et stlpar : `BehaviourVariable1`, `BehaviourVariableFactory1/2`,
  `ExternalModelCaller2/4` identiques à `generic`, et une variante à
  paramètres non statiques.
- Tests ajoutés aux 8 configurations : `BehaviourVariableTest` (fabrique,
  loi interne `IsotropicPlasticMisesFlow` avec unités et paramètres non
  statiques) et `ExternalModelCallerTest` (Runge-Kutta rk54 + modèle
  externe à paramètre).

### 2026-10-05 (nuit, suite), brique StandardElastoViscoPlasticity et en-têtes de TFEL

- Balayage CUDA (mfront de développement) des 140 lois de
  `StandardElastoViscoPlasticity` : 118 OK au premier passage, 140/140
  après corrections :
  - critères GTN et Rousselier, prédiction élastique de la déformation
    axiale (grandes déformations logarithmiques) : `throw` en cas de
    non-convergence → en code device, NaN (l'intégration échoue), comme les
    propriétés matériau ;
  - schéma étagé de la porosité (`tfel::raise<DivergenceException>` dans
    `@AdditionalConvergenceChecks`, qui ne peut pas renvoyer d'échec) : en
    code device, le compteur d'itérations du solveur est porté à `iterMax`
    pour arrêter la résolution sur un échec. `TinyNonLinearSolverBase`
    compare désormais avec `>=`/`<` au lieu de `==`/`!=` (identique en
    temps normal, bien défini si le compteur est forcé) ;
  - fonctions générées `mfront_computeYieldRadius*` (écrouissage isotrope
    défini par des données) sans `TFEL_HOST_DEVICE` ; **défaut amont** : le
    nom de `mfront_computeYieldRadiusAndDerivative` était doublé
    (`m2 + local_id`) dans la variante spline ;
  - `AitkenAccelerationAlgorithm`, interpolations linéaire et spline
    cubique : annotations ;
  - **défaut amont** : déclarations à 2 arguments des fonctions de
    nucléation Chu-Needleman (déformation) jamais définies, les définitions
    à 4 arguments utilisées par le code généré n'étant pas déclarées (donc
    hôte seulement) ;
  - **défaut nvcc 13.0** (cudafe++) : `auto i = decltype(s){};` avec
    `const auto s = c.size();` est réécrit pour le compilateur hôte en
    `const unsigned long{}`, invalide → `size_type{}` dans
    `LinearInterpolation.ixx` (plus lisible de toute façon).
- **Balayage clang CUDA des en-têtes** (`scratch/hdscan`) : chaque en-tête
  de `TFEL/Math` et `TFEL/Material` (335) compilé seul par clang en mode
  CUDA. clang refuse une définition sans `__host__ __device__` dont la
  déclaration en a (nvcc l'accepte) : 60 + 25 définitions corrigées
  (critères de Barlat, Cazacu, Drucker, Mohr-Coulomb, Michel-Suquet,
  Rousselier, GTN, plan pi, nucléations, `DDIF2Base`, solveurs de Broyden,
  Levenberg-Marquardt, Powell, bissection, Newton scalaire, schémas
  d'homogénéisation PCW) par un script (`fix_hd.py`) puis vérifiées ;
  attribut `[[gnu::visibility]]` placé après `__host__ __device__` dans
  `DDIF2Base` ; 19 en-têtes utilisant `TFEL_HOST_DEVICE` sans inclure
  `TFEL/Config/TFELConfig.hxx`. Plus aucun défaut de ce type.
- **Faux défaut GPU élucidé** : `OrthotropicNortonTest` (Hill) donnait
  66 % d'écart sur GPU (nvcc et clang) à 3e-3 de déformation, mais exact
  sur CPU, y compris avec `-DTFEL_DEVICE_CODE` et
  `-ftrivial-auto-var-init=pattern`. Le même écart apparaît sur CPU avec
  `-march=native -ffp-contract=fast` (FMA) : problème mal conditionné
  (`@Epsilon 1e-16`, quelques points en échec), Newton converge ailleurs
  selon les arrondis. Ce n'est pas un défaut du GPU. Critère retenu pour
  les lois de test : l'interface `generic` et `generic-parallel` compilée
  avec FMA sur CPU doivent concorder avant toute comparaison GPU.
- Tests ajoutés : `StandardElastoViscoPlasticityTest` (plastique, Mises,
  écrouissage isotrope par données, Armstrong-Frederick + Chaboche 2012),
  `OrthotropicStandardElastoViscoPlasticityTest` (Norton, Hill, Voce +
  linéaire), `PorousStandardElastoViscoPlasticityTest` (GTN, nucléation
  Chu-Needleman, schéma étagé, porosité initiale 0,005) : robustes (FMA)
  et identiques à 1e-15 près sur CUDA.

### 2026-10-05 (nuit, suite), commits 0016 à 0018

- Régression TFEL : 6965/6965 avec tout ce qui précède (installation
  fraîche). Commits locaux : 7b1ed663c (TFEL : annotations cohérentes,
  critères et solveurs), 4aa867593 (MFront : variables de comportement,
  modèles, écoulements poreux), bd5528869 (generic-parallel : modèles et
  tests). Paquet `backups/tfel-generic-parallel-behaviours-2026-10-05.bundle`
  vérifié, `proto-patches/` régénéré (0001 à 0018).
- Matrice O2 sur l'installation fraîche en cours
  (`results/generic-parallel-tests-2026-10-05-c.txt`).
- Script `scripts/clang-cuda-header-scan.sh` : contrôle des annotations
  hôte/device de tous les en-têtes par clang CUDA (0 incohérence).

### 2026-10-05 (soir), C5 backend OpenMP, points d'entrée factorisés

- Points d'entrée communs dans `Integrate.hxx` : `integrateBatch`,
  `executeInitializeFunctionOnBatch`, `executePostProcessingOnBatch`,
  `rotateBatch` (contrôle des données, choix de l'objet fonction, appel de
  la primitive `apply` du backend passée comme objet appelable). Les
  backends stlpar, CUDA/HIP et Kokkos ne font plus que lier leur primitive
  à leurs options (`internals::makeApply`) ; SYCL garde ses points d'entrée
  (copies vers le device avant le choix). stlpar et CUDA : 46/46
  comparaisons identiques après la factorisation.
- Backend `openmp` (`OpenMP.hxx`, `OpenMPBackend`) : `target teams
  distribute parallel for map(to: f)` avec réductions membre à membre des
  statistiques (`+` pour les nombres d'échecs et de résultats peu fiables,
  `min` pour `rdt`) ou `parallel for` sur l'hôte (option `offload`).
  `TFEL_DEVICE_CODE` défini pour `nvc++ -mp=gpu` (`__NVCOMPILER_OPENMP_GPU`,
  compilateur mono-passe comme `-stdpar=gpu`), `getConstantObject` renvoie
  une copie comme avec `-stdpar=gpu`.
- **Défaut nvc++ 25.11 (OpenMP, -O2 et -O3)** : les lois des DSL implicites
  échouent dès qu'un opérateur tangent est demandé (même élastique) ; le
  Newton converge sur GPU, mais `computeConsistentTangentOperator` renvoie
  faux (décomposition LU du jacobien). Juste à -O0 et -O1 ; un `printf`
  dans la branche d'échec suffit à corriger le résultat (défaut de
  l'optimiseur). Pas reproduit hors contexte (LU 7×7 seule dans une région
  `target` : juste), ni corrigé par `-Mnounroll`, `-Mnovect`, `-Mnoinline`,
  `-Mnofma`, ni par `map(to:)` au lieu de `firstprivate`. Le même code est
  juste avec `nvc++ -stdpar=gpu` à -O2. Contournement documenté : la
  configuration `openmp-offload-nvhpc-gpu` est compilée à -O1. clang ne
  peut pas servir de compilateur de référence ici : le LLVM 20 de spack
  est construit sans le support du déport OpenMP (`offload=False`).

### 2026-10-05 (soir, suite), propriétés matériau : backends Kokkos et OpenMP

- L'interface `generic-parallel` des propriétés matériau
  n'avait que CUDA, HIP, stlpar et SYCL. Ajout de `kokkos` (option
  `execution_space`, `Kokkos::parallel_for`, statuts des bornes dans une
  `Kokkos::View` en `SharedSpace`) et `openmp` (options `offload`,
  `device_stack_size`, mêmes constructions que pour les lois).
- Factorisation dans `CxxProgrammingModelBackendBase` : choix du noyau
  selon les pas (`writeKernelSelection`, auparavant recopié dans stlpar et
  SYCL), écriture des noyaux (`writeKernels`), spécificateurs des lambdas
  (`getKernelSpecifiers`, `KOKKOS_FUNCTION` pour Kokkos). nvcc interdit les
  initialiseurs dans les captures des lambdas étendues (`__host__
  __device__`) : le pointeur vers les statuts des bornes et les paramètres
  sont désormais déclarés avant les noyaux et capturés par nom (conteneur
  `mfront_bounds_statuses_storage` dans chaque backend).
- Test commun `test-cxx.cxx` (ancien `test-stlpar.cxx`) pour stlpar,
  openmp et kokkos : backend choisi par `MFRONT_GPMP_TEST_BACKEND`, tableaux
  en mémoire partagée Kokkos si `MFRONT_GPMP_TEST_USE_KOKKOS`. Résultats :
  stlpar g++, OpenMP hôte, OpenMP déporté (nvc++ -O2, 19 noyaux GPU),
  Kokkos CUDA et Kokkos hôte : tous les tests réussis. Dossiers de tests
  `kokkos`, `kokkos-host`, `openmp-host-gcc`, `openmp-offload-nvhpc-gpu`
  ajoutés ; plus aucune configuration réservée aux lois.

### 2026-10-05 (soir, suite), matrice c, OpenMP, modèles compilés par nvc++

- Matrice O2 (`results/generic-parallel-tests-2026-10-05-c.txt`, état des
  commits 0016-0018) : CUDA, HIP, Kokkos, Kokkos hôte, stlpar g++ réussis
  avec les sous-tests TDLS ; `stlpar-parunseq-nvhpc-gpu` : 280 comparaisons
  justes et 4 fausses, toutes sur le cas d'échec des modèles
  (`GenericParallelModelTest_PlaneStress` avec perturbation). Cause : la
  bibliothèque de référence (`generic`) est elle aussi compilée par
  `nvc++ -stdpar=gpu`, où `TFEL_DEVICE_CODE` est défini dans toute l'unité
  de compilation : les contrôles de bornes renvoient alors un booléen au
  lieu de lever une exception, et les fonctions hôte des modèles (`void`)
  l'ignoraient (bornes physiques ignorées en silence). Corrigé :
  `writeGenericModelClass` écrit toujours des fonctions `bool` qui
  renvoient faux si un contrôle échoue et une méthode `execute`, que
  l'implémentation hôte utilise (échec signalé, `rdt = 0.1`). Sur un hôte
  normal, rien ne change (les exceptions sont levées). 12/12 avec nvc++.
- OpenMP hôte (g++ `-fopenmp`) : `offload` devient un paramètre de patron
  (`if constexpr`) : le code hôte ne dépend plus des points d'entrée de
  déport du support d'exécution OpenMP (`GOMP_teams4` absent des
  bibliothèques `libgomp.so.1` de NVHPC et de LLVM présentes dans
  `LD_LIBRARY_PATH`). Suite complète : 284/284.
- OpenMP déporté (nvc++ -O1) : chaque loi est juste séparément, mais la
  bibliothèque de toute la suite ne se charge pas sur GPU (« Could not find
  symbol _st_... on the device »). Bissection par réédition des liens des
  objets existants : `HosfordPerfectPlasticity` seule en cause ; son code
  device (IR `.gpu` gardé par `-gpu=keep`) définit un bloc statique interne
  de 317 octets jamais utilisé, supprimé par l'optimiseur, que le support
  d'exécution cherche au chargement. Défaut de nvc++ 25.11 ; avec celui de
  -O2, la configuration `openmp-offload-nvhpc-gpu` reste fournie mais n'est
  pas activée (documenté dans la page de l'interface et dans
  `results/nvhpc-bugs/README.md`).
- Fin de la matrice c : SYCL GPU et SYCL CPU (AdaptiveCpp, compilateur
  générique mono-passe) échouent sur le même cas des modèles que nvc++
  (même cause, corrigée par les fonctions `bool`) ; toutes les propriétés
  matériau réussissent ; régression 6965/6965.
- `scripts/run-generic-parallel-tests.sh` : tfel-check signale comme réussi
  un test dont les exigences (`@Requires`) ne sont pas remplies (toutes ses
  commandes sautées) ; le script compte désormais séparément les tests
  réussis et sautés (configurations non activées, TDLS absent).
- Chaîne de validation lancée (`scratch/validate2.sh`) : reconfiguration
  avec `openmp-host-gcc`, installation, régression, matrices -O2, -O1, -O3
  (`results/generic-parallel-tests-2026-10-05-d-O*.txt`), balayage CUDA de
  tous les tests (`scratch/sweep-config-cuda-nvcc/summary.tsv`).
- Régression après installation (validate2) : 6965/6965. Commits locaux
  4da85b2c5 (TFEL : `TFEL_DEVICE_CODE` avec nvc++ -mp=gpu), b0290dce3
  (MFront : contrôles de bornes des modèles), b401968ce (generic-parallel :
  OpenMP, Kokkos pour les propriétés matériau, points d'entrée communs) ;
  paquet de sauvegarde et `proto-patches/` (0001 à 0021) à jour.
- Étude MGIS (D1, lecture seule) : `scratch/d1-mgis-design.md`. Constat
  important : MGIS charge aujourd'hui sans erreur une bibliothèque
  `generic-parallel` isotrope (mêmes noms de fonctions, signatures
  différentes) et l'appel est alors un comportement indéfini ; il faut
  détecter l'interface (`_mfront_interface` == "GenericParallel").
  Recommandation : option B (MGIS assemble un champ par groupe, sans copie
  dans les cas courants), mémoire via `std::pmr::memory_resource`, nouvelle
  surcharge `integrate(BatchExecutionPolicy, MaterialDataManager&, ...)`.

### 2026-10-06 (nuit), versions de NVHPC et de CUDA, défauts d'OpenMP déporté

- Notes de version de NVHPC 26.1 à 26.9 : aucune liste de défauts corrigés.
  Reproducteurs rejoués avec NVHPC 26.3 et 26.9 installés dans
  `deps/nvhpc/` (archives vérifiées par les empreintes du paquet spack) :
  seul `aggregate-zero-init` est corrigé (26.9) ; `run-nvcxx-reproducers.sh`
  et le nouveau reproducteur `std-variant-visit.cxx` dans
  `results/nvhpc-bugs/`. Défauts de nvcc présents de 13.0.88 à 13.4.92
  (compilateurs des paquets pip `nvidia-cuda-nvcc`).
- **Hosford sous OpenMP déporté** (symbole `_st_…__T1_0` introuvable) : le
  bloc de 317 octets est le texte `__PRETTY_FUNCTION__` des `assert` de
  `StensorComputeEigenVectors<3>::computeEigenVectors` ; nvc++ supprime
  l'`assert` du code device mais garde l'enregistrement du texte. Les trois
  `assert` étaient redondants (impliqués par les conditions des branches),
  l'un était faux à l'égalité (`<` au lieu de `<=`) : remplacés par des
  commentaires énonçant les invariants. Ce sont les seuls `assert` des
  en-têtes de TFEL utilisés dans du code device.
- **Opérateur tangent faux à -O2/-O3 sous OpenMP déporté** : ni la LU, ni
  la pile (cadres identiques à -O1 et -O2, ptxas), ni une variable non
  initialisée (`-ftrivial-auto-var-init`). Marqueurs distincts de `rdt` sur
  chaque chemin d'échec : `computeConsistentTangentOperator` renvoie faux
  sans passer par aucun de ses `return false`. Déclencheur : l'objet local
  `GetPartialJacobianInvert` généré par MFront, qui garde une référence
  `bool&` vers le drapeau local `mfront_success` (sa seule construction
  suffit). Générateur corrigé (`ImplicitCodeGeneratorBase`) : chaque objet
  (`GetPartialJacobianInvert`, `GetIntegrationVariablesDerivatives_*`)
  possède son drapeau, combiné dans `mfront_success` après le code de
  l'utilisateur ; la substitution passe par une variable locale (un échec
  n'est plus effacé par un appel suivant réussi).
- icpx (oneAPI 2026.1.1, apt) installé : SYCL sur CPU seulement (greffon
  NVIDIA retiré depuis 2025.3) ; `libintelocl` exige le TBB de oneAPI
  (2023.1) : chemins `compiler/2026.1/lib`, `umf/1.1/lib`,
  `tbb/2023.1/lib/intel64/gcc4.8` à ajouter en fin de `LD_LIBRARY_PATH`.
- clang CUDA : clang 20 échoue avec CUDA 13.0 (`texture_fetch_functions.h`
  supprimé) et avec le CUDA 12.9 de NVHPC (en-têtes curand dans
  `math_libs`). Configurations clang sur GPU reportées au nouvel
  environnement (LLVM 23, CUDA 13) ; préfixe CUDA regroupant `math_libs`
  nécessaire si le CUDA vient de NVHPC.
- Tests des lois ajoutés pour `stlpar-parunseq-clang`,
  `stlpar-parunseq-icpx` et `sycl-default-icpx` (à valider).

### 2026-10-06 (nuit, suite), champs par variable, symboles, rôles

- `BatchData` : les propriétés matériau et les variables externes peuvent
  être données par groupe (un champ) ou par variable (un champ par
  variable, membres optionnels `material_properties_by_variable` et
  `external_state_variables_by_variable`). Une structure initialisée
  comme avant garde son sens. Les deux formes sont converties côté hôte
  dans `InputFields` (champs et nombres de composantes), passé par valeur
  aux noyaux : un seul chemin de code par point. Les propriétés matériau
  et les variables externes d'un point sont toujours copiées, le mode
  direct ne dépend plus que des autres champs (une propriété uniforme ne
  le désactive plus). SYCL copie chaque variable si nécessaire.
  `PointTraits` donne le nombre de composantes de chaque variable (lois
  et modèles), lu seulement côté hôte.
- `compare.cxx` : troisième rangement « par variable » (valeurs
  uniformes, pas nul) pour toutes les comparaisons : 396 au lieu de 284
  avec stdpar g++.
- Symboles : `<loi>_generic_parallel_backend` (lois et modèles).
  Défaut corrigé : les modèles `generic-parallel` exportaient
  `mfront_interface = "Generic"`.
- Rôles de TFEL/MFront, MGIS et mfem-mgis : section de `DESIGN.md`.
  Conception MGIS (`scratch/d1-mgis-design.md`, 2.3) : MGIS passe un champ
  par variable, sans espace de travail.
- Développement dans l'arbre de travail `tfel-dev` (branche
  `batchdata-by-variable`) pendant la validation de `tfel/` ; MFront de
  développement dans `scratch/build-dev2`, tests par
  `scratch/dev-test.sh <configuration>`.

### 2026-10-06 (nuit, suite), pilotes de l'interface generic, OpenMP en mode direct

- Test de développement OpenMP sur GPU (-O2) des champs par variable : un
  seul échec, la fonction d'initialisation de SaintVenantKirchhoffElasticity
  en mode direct (champs point par point). La contrainte écrite est perdue
  (le lot garde la contrainte initiale). Juste à -O1. Déplacer une copie
  dans `PointView` fait échouer les trois rangements : défaut de
  l'optimiseur sensible à la forme du code, comme celui de nvcc 13.0 sur la
  même fonction.
- Déclencheur : les pilotes de l'interface `generic` (intégration,
  fonctions d'initialisation, post-traitements des lois en petites
  déformations, en grandes déformations, en déformation de Green-Lagrange
  et logarithmique) modifiaient temporairement les pointeurs des données
  de la loi données par l'appelant, puis les restauraient. Ils travaillent
  désormais sur une copie locale des données dont les pointeurs désignent
  les valeurs locales ; pour la contrainte de Cauchy, les pilotes en
  grandes déformations appellent directement ceux des petites
  déformations. Résultats inchangés ; OpenMP sur GPU -O2 : 396/396.
- Défaut latent corrigé au passage : la contrainte calculée par une
  fonction d'initialisation n'était exportée que pour la contrainte de
  Cauchy des lois en grandes déformations standard. Elle est maintenant
  convertie dans la mesure de l'appelant (gradient de la transformation
  du début du pas) pour toutes les mesures et les lois écrites en
  déformation de Green-Lagrange ou logarithmique. Conversion commune dans
  `MFront/GenericBehaviour/StressMeasures.hxx`. Tests ajoutés (PK2, PK1).
- Structures C de `BatchData.h` nommées (`typedef struct nom {...} nom;`) :
  cudafe (nvcc 13.0) plantait sur le destructeur d'un `std::vector` de
  champs anonymes.
- Matrice -O2 de la nuit (commits du 05/10 soir et des corrections OpenMP) :
  OpenMP sur GPU et stdpar icpx réussissent ; stdpar clang échoue à cause
  de l'environnement : le clang du système charge la bibliothèque du LLVM
  de spack (`LD_LIBRARY_PATH`), construite sans la cible AArch64. Le clang
  de ce LLVM est désormais utilisé.

### 2026-10-06 (matin), configurations clang CUDA et HIP

- Les lois sont maintenant testées avec `cuda-clang` et `hip-clang`
  (plateforme NVIDIA). Seules les propriétés matériau l'étaient.
- clang 20 ne connaît que les CUDA sortis avant lui (jusqu'à 12.8). La
  variable CMake `TFEL_CLANG_CUDA_PATH` choisit le CUDA de clang (par
  défaut celui de nvcc). Ici : CUDA 12.9 du HPC SDK (préfixe fusionné de
  `deps/spack/make-nvhpc-cuda-prefix.sh`).
- clang refuse tout appel d'une fonction hôte depuis le code GPU. nvcc
  accepte les fonctions hôte `constexpr` avec `--expt-relaxed-constexpr`.
  Deux trous trouvés ainsi : `pull_back` et `invert(st2tost2)`, appelés par
  les conversions des lois en grandes déformations. Ils sont désormais
  `TFEL_HOST_DEVICE` (`invert` s'appuie sur `TinyMatrixInvert`, déjà
  utilisable sur GPU). Il ne reste qu'une fonction hôte `constexpr` dans
  TFEL : `tensor::operator()(i, j)`, qui lève une exception.
- HIP du système dans `/usr` : `-isystem /usr/include` plaçait ce dossier
  avant les en-têtes de libstdc++ et cassait leurs `#include_next`. Les
  configurations clang passent le dossier de HIP avec `-I`, que les
  compilateurs ignorent pour un dossier système.
- Sous-tests TDLS ajoutés aux cinq configurations qui n'en avaient pas
  (cuda-clang, hip-clang, stlpar clang, stlpar icpx, SYCL icpx).
- Résultats (arbre de développement, -O2) : cuda-clang et hip-clang
  402/402 chacun, dont 114 avec les champs par variable ; TDLS 18/18 pour
  les cinq configurations.
- `scripts/compile-sweep.sh` reconnaît les erreurs de clang. Balayage de
  toutes les lois, modèles et propriétés de MFront avec cuda-clang lancé.
- Branche `clang-cuda` de `tfel-dev` (3 commits), à intégrer après la
  validation en cours.

### 2026-10-06 (matin, suite), balayage clang CUDA de tous les tests de MFront

- Balayage des 608 lois, modèles et propriétés de `mfront/tests` avec
  cuda-clang (MFront et en-têtes de développement, 21 min) : 458 OK au
  premier passage. clang refuse tout appel hôte depuis le code GPU, il
  complète donc le balayage nvcc.
- Appels hôte (38) : 37 dans le code écrit par l'utilisateur des lois de
  test (`throw`, `std::cout`, `tfel::raise`). Un seul trou de TFEL : la
  vue `tensor` d'un `stensor`, utilisée par tout produit ou somme mêlant
  `tensor` et `stensor` (par exemple `F * S`), et les accesseurs `(i, j)`
  des tenseurs étaient réservés à l'hôte, car ils lèvent une exception
  sur un indice invalide. Nouvelle fonction
  `internals::reportInvalidTensorIndex` : exception sur l'hôte, NaN sur
  GPU, comme `TinyMatrixInvert`. Plus aucune fonction `constexpr`
  réservée à l'hôte dans `TFEL/Math`.
- Loi de test ajoutée à toutes les configurations :
  `DefaultGenericSaintVenantKirchhoffElasticity` (gradient F, force PK1,
  produit `tensor` × `stensor`, conversion de l'opérateur tangent).
  nvcc et clang : 18/18.
- Erreurs de compilation (32) :
  - 26 tests unitaires `WarningTest*`, qui ne doivent pas compiler. Le
    balayage compile désormais aussi l'interface `generic` avant de
    conclure : ils sont classés `INVALID` ;
  - propriétés matériau définies par interpolation de données (`@Data`,
    5 fichiers) : le noyau CUDA/HIP nommait ses champs d'entrée
    `mfront_<entrée>_values`, nom des abscisses des données. Renommés
    `mfront_<entrée>_field`. Échouait aussi avec nvcc ;
  - propriété matériau appelant une autre loi (`@MaterialLaw`) : clang
    analyse aussi le corps d'un noyau `__global__` pendant la passe hôte,
    où seules les versions hôte des lois appelées étaient visibles. Les
    noyaux importent désormais toujours les implémentations
    `mfront::device`.
- Tests ajoutés aux 16 configurations des propriétés matériau :
  `LinearDataInterpolationTest` et `MaterialLawCallTest`, vérifiés par une
  fonction commune dans chaque programme de test. Défaut latent corrigé :
  `hipMallocManaged` reçoit une taille en octets (4 octets alloués pour
  4 réels). Six configurations essayées (nvcc, clang, hipcc, hip-clang,
  OpenMP, stdpar) : toutes réussissent.
- « Suspects » (39) : tous des propriétés matériau, dont les noyaux sont
  petits avec clang. Le critère des registres ne s'applique plus qu'aux
  lois et aux modèles.
- Branche `clang-cuda` : 8 commits au-dessus de la branche validée.

### 2026-10-06 (matin, suite), DSL non exécutés sur GPU et deux défauts de nvc++ OpenMP

- Couverture à l'exécution : quatre DSL de MFront n'étaient que compilés
  par le balayage, jamais exécutés sur GPU. Lois de test ajoutées à toutes
  les configurations : `IsotropicMisesCreepTest`,
  `MultipleIsotropicMisesFlowsTest` (écrouissage de Voce et fluage),
  `ImplicitGenericBehaviourTest` et `ImplicitModelTest`. CPU, nvcc, clang,
  hipcc : 60/60 ; SYCL AdaptiveCpp : 84/84.
- Limite amont, sans lien avec le GPU : le DSL
  `MultipleIsotropicMisesFlows` ne fournit pas d'opérateur tangent
  cohérent, et le pilote de l'interface `generic` refuse alors toute
  demande d'opérateur, même élastique. Aucun test amont n'exécute une loi
  de ce DSL (elles sont seulement compilées). La comparaison porte sur les
  intégrations sans opérateur et sur la prédiction.
- nvc++ 25.11 OpenMP sur GPU, -O2 et -O3, deux défauts trouvés par la
  loi `DefaultGenericSaintVenantKirchhoffElasticity` (produit `F * S`) :
  - erreur interne « Computed goto control expression must be a 32-bit
    integer » : le `switch` de la vue `tensor` d'un `stensor`, inliné dans
    celui du produit. Contournement : vue écrite avec une table d'indices
    et une comparaison, sans `switch` ;
  - code du noyau supprimé, aucun résultat écrit : la mise à zéro d'un
    membre tableau (`v()`, `v{}`, `= {}`) est réalisée par la copie d'un
    tableau global `.common` non initialisé. Reproducteur autonome de 30
    lignes. Contournement : mise à zéro par une boucle dans le
    constructeur par défaut de `GenericFixedSizeArray`.
  Reproducteurs et contournements dans `results/nvhpc-bugs/README.md`.
  Méthode : la sentinelle compte. Avec `map(from:)`, les cases non
  écrites rendent la mémoire du GPU d'un essai précédent ; seul
  `map(tofrom:)` avec une valeur sentinelle distingue un noyau vide.
- Après contournements, les cinq nouvelles lois sur nvc++ : OpenMP sur GPU
  78/78 à -O1, -O2 et -O3, stdpar 78/78. Tests unitaires `Math` de TFEL
  (hôte) : 80/80.
- `scratch/dev-test.sh` : variable `DEVTEST_OPTIMIZATION_LEVEL`, comme la
  matrice.
- Branche `clang-cuda` : 10 commits au-dessus de la branche validée.

### 2026-10-06 (matin, fin), validation de la branche fusionnée

- Les matrices ne tournent plus qu'à -O3 : quand -O3 passe, -O1 et -O2
  passent dans la très grande majorité des cas. `VALIDATE_LEVELS` choisit
  d'autres niveaux si besoin.
- Branche `clang-cuda` fusionnée dans `generic-parallel-behaviours`
  (9cae5c235). `scratch/validate3.sh 2026-10-06c` :
  - régression de TFEL : 6965/6965 ;
  - matrice -O3 : les 31 tests réussissent (15 configurations de lois,
    toutes avec leur sous-test TDLS, 16 de propriétés matériau), aucun
    sauté (`results/generic-parallel-tests-2026-10-06c-O3.txt`) ;
  - balayages de compilation des 608 lois, modèles et propriétés,
    identiques avec nvcc et clang : 504 OK, 36 code utilisateur, 60
    invalides (dont 26 tests unitaires qui ne doivent pas compiler), 7
    fragments, 1 appel hôte (`results/compile-sweep-*-2026-10-06.tsv`).
- Le dernier appel hôte venait de MFront : la comparaison à la jacobienne
  numérique (`@CompareToNumericalJacobian`) affiche ses résultats avec
  `std::cout`. Elle n'est plus générée que pour l'hôte (53d30ff6c).
- Les 36 cas de code utilisateur sont tous dans les lois de test :
  `throw`, `tfel::raise`, `std::cout`, statiques initialisés à
  l'exécution, systèmes de glissement des polycristaux lus dans un
  fichier (limite documentée).

### 2026-10-06 (après-midi), corrections de la relecture (R1)

Branche `review-fixes` de TFEL, au-dessus de 53d30ff6c :

- A4 : la classe de loi a la même définition sur l'hôte et le GPU
  (propriété d'interface `requiresDeviceCompatibleBehaviourClass`, option
  de DSL `device_compatible_class` propagée aux variables de comportement).
- A1 : la référence des tests est compilée par le compilateur par défaut
  de MFront (`build-libraries.sh`), indépendamment de la configuration
  testée.
- A2 : les paramètres des variables de comportement sont comparés sur
  l'hôte à leur valeur par défaut, utilisée sur le GPU (statut -2 sinon).
  Les paramètres des lois matériau appelées sont des constantes sur l'hôte
  et le GPU : ce point de la relecture était faux.
- B2 : une seule implémentation des points d'entrée. Le code généré
  appelle `integrateBatch` (et les trois autres) avec la primitive
  `makeApply` du backend ; SYCL ajoute ses transferts. 677 lignes de moins.
- A3 : options d'exécution (`mfront_gpb_ExecutionOptions`) : file,
  statistiques asynchrones (CUDA, HIP, Kokkos, SYCL sans copie),
  `preserve_stack_size`. Un seul mutex de pile par bibliothèque. Le second
  pas des comparaisons est asynchrone.
- B4 : `MFRONT_GPB_API_VERSION`, exporté par les bibliothèques. Les
  messages par point ne sont pas transmis (pas d'exceptions ni de tampon
  par thread sur GPU) : la loi `generic` donne le message d'un point.
- B1 : options communes aux deux interfaces (`BackendOptions`), espace
  d'exécution Kokkos validé, conventions de statut documentées.
- B5 : `compare` lit les tailles dans la bibliothèque de référence et
  déduit la fonction testée de son nom. Les tailles lues sont identiques
  aux tailles données à la main pour les 75 comparaisons.
- A5 : vecteurs propres d'une valeur propre multiple par Jacobi sur GPU
  seulement (résultats CPU inchangés), incompatibilités dans les notes de
  version, effet des compilateurs à passe unique documenté, principe 4 de
  DESIGN.md corrigé, plan de découpage en PR (`UPSTREAM.md`).
- B3 et C : déduplication du générateur documentée, statut initialisé
  pour un lot nul, politique hors bornes des modèles, `GenTypeBase`
  restauré, surcharge de `callAndReportExceptions`, utilitaire hôte/GPU
  dans le DSL Runge-Kutta, liste CMake vide retirée.
  `getConstantObject` garde son type de retour : copie imposée par nvc++,
  référence ailleurs, API inchangée pour les autres compilateurs.

- B7 : accès directs aux champs rangés point par point contre copies
  locales seules (`scratch/b7-measure.sh`, couche point par point) :

  | | compilation nvcc | compilation gcc | exécution |
  |---|---|---|---|
  | SingleCrystal_DD_FCC | 7,9 s contre 5,8 s | 2,6 s contre 2,6 s | identique |
  | ImplicitNorton | 16,8 s contre 10,6 s | 5,3 s contre 4,4 s | identique |

  Le temps d'exécution est le même à 1 % près sur GPU (RTX 5080) ; sur
  CPU, le bruit d'OpenMP (16 fils) masque toute différence. Le mode direct
  est retiré : une seule variante des noyaux.

Tests à -O3 (`scratch/dev-test.sh`) : 489/489 pour cuda-nvcc, cuda-clang,
hip-hipcc, kokkos, kokkos-host, openmp-host-gcc,
openmp-offload-nvhpc-gpu, stlpar-parunseq-gcc,
stlpar-parunseq-nvhpc-gpu ; 513/513 pour sycl-default-acpp (avec les
copies depuis la mémoire hôte) ; TDLS 18/18 ; propriétés matériau OK sur
cinq configurations.

### 2026-10-06 (soir), lois de test portables (R2)

Branche `portable-test-laws` de TFEL, au-dessus de 7c1427f57.

- `tfel::fail<Exception>(...)` : lance l'exception sur l'hôte comme avant,
  renvoie `false` sur GPU. Le bloc de code renvoie cette valeur.
- Renvoyer `false` ne fait pas échouer l'intégration partout : dans
  `@Integrator` des DSL implicites, la correction de Newton est divisée
  par deux ; dans `@Derivative` des DSL Runge-Kutta, le pas de temps est
  réduit. Une exception, elle, arrête tout. Sur GPU, le même effet
  s'obtient par `this->abortResolution()` (solveur non linéaire de TFEL,
  compteur d'itérations au maximum) et `this->abortIntegration()` (lois
  Runge-Kutta, générée par MFront).
- `@SharedData Type nom = expression;` : donnée construite une fois sur
  l'hôte (orientations des grains d'un polycristal lues dans un fichier).
  La loi y accède par un membre référence. L'interface generic-parallel
  la copie au premier appel dans une mémoire accessible au GPU (mémoire du
  GPU pour CUDA, HIP, SYCL et Kokkos, tas sinon) et passe son adresse aux
  noyaux par la structure des paramètres. `PolyCrystalsSlidingSystems` en
  tableaux de taille fixe. Interdite dans les variables de comportement.
- 36 lois de `mfront/tests` converties sans changer leurs résultats sur
  l'hôte (mêmes exceptions, mêmes messages) : 30 à exceptions ou `cout`,
  6 polycristaux. Les messages construits à l'exécution et les `cout`
  restent sur l'hôte (`#ifndef TFEL_DEVICE_CODE`).
- Nouveau test de comparaison : polycristal de Berveiller-Zaoui à quatre
  grains (`PolyCrystalSharedDataTest`).
- Défauts d'environnement de NVHPC 25.11 trouvés en chemin (dans
  `results/nvhpc-bugs/README.md`) : `-mp=gpu -gpu=mem:managed` ne lie pas
  `libnvhpcmanaux` (contrairement à `-stdpar=gpu`), d'où un
  `free(): invalid pointer` quand du code g++ libère de la mémoire gérée ;
  `libnvcpumath.so` ne dépend pas de `libm` : un programme qui ne charge
  pas `libm` ne peut pas charger une bibliothèque nvc++ (sans effet sur
  MGIS).
- Balayages de compilation des 608 fichiers avec nvcc et clang : 541 OK,
  0 code utilisateur, 0 appel hôte (contre 504 OK, 36 code utilisateur et
  1 appel hôte le 06/10 au matin), 60 invalides et 7 fragments inchangés.
- Le pilote de l'interface `generic` construisait la loi hors de la zone
  qui convertit les exceptions en échecs. Avec `@SharedData`, la lecture
  du fichier d'orientations a lieu dans le constructeur : un fichier
  manquant aurait fait sortir l'exception de la fonction générée au lieu
  de renvoyer un échec comme avant. La loi est maintenant construite dans
  la zone protégée (vérifié : statut -1 et message de l'exception).
- `compare` accepte le pas de temps et les valeurs des variables externes
  (pour R3).
- Validation à -O3 : 498/498 sur cuda-nvcc, cuda-clang, hip-hipcc, kokkos,
  kokkos-host, openmp-host-gcc, openmp-offload-nvhpc-gpu,
  stlpar-parunseq-gcc, stlpar-parunseq-nvhpc-gpu ; 522/522 sur
  sycl-default-acpp ; TDLS 18/18 ; propriétés matériau OK sur cinq
  configurations. Après le correctif du pilote `generic` : openmp-host-gcc,
  cuda-nvcc et openmp-offload-nvhpc-gpu à nouveau 498/498,
  sycl-default-acpp 522/522.

### 2026-10-07 (nuit), essais d'exécution de toutes les lois de test (R3)

Les lois de `mfront/tests` qui compilent avec l'interface `generic-parallel`
sont maintenant comparées à l'exécution, et plus seulement compilées :

- `scripts/runtime-sweep.py` (lois de comportement et modèles) : pour
  chaque loi, bibliothèques `generic` (référence) et `generic-parallel`,
  données lues dans le fichier MTest de la loi (propriétés matériau,
  variables externes, valeurs initiales des variables internes,
  paramètres, pas de temps, amplitude du chargement ; définitions évaluées
  dans n'importe quel ordre, évolutions comprises), puis `compare`. Sans
  MTest, les propriétés matériau prennent des valeurs synthétiques
  choisies d'après leur nom. Si la référence échoue en tous les points,
  des pas de temps et des amplitudes plus petits sont essayés ;
- `scripts/property-sweep.py` (propriétés matériau) : fonction `generic`
  appelée point par point contre les deux prototypes `generic-parallel`,
  avec un argument uniforme, et avec des valeurs hors bornes et la
  politique stricte ;
- `scripts/roundoff-check.sh` : confirme qu'un petit écart vient des
  arrondis (bibliothèque du lot compilée sur l'hôte avec contraction FMA,
  comme nvcc sur le GPU, contre une référence sans FMA).

Configurations : openmp-host-gcc et cuda-nvcc (comparaisons une par une
sur le GPU).

Défauts de TFEL trouvés et corrigés :

- tableaux de gradients et de forces thermodynamiques (`@Gradient Stensor
  a[2]`) : la taille d'un point ignorait la taille des tableaux, seule la
  première composante était recopiée (`getPointDataSizes`) ;
- conversion de l'opérateur tangent des lois en grandes transformations :
  le code généré lit l'opérateur calculé par l'utilisateur avec
  `Dt.template get<T>()`, qui ne vérifie le type que sur l'hôte. Si le
  bloc de l'utilisateur ne calcule pas le bon opérateur, l'hôte échoue
  (exception GenTypeCastError), mais le GPU renvoyait un succès avec une
  matrice fausse. Le type est maintenant testé avant la lecture : même
  exception sur l'hôte, échec sur le GPU ;
- paramètres des modèles partagés entre bibliothèques (code amont de
  l'interface `generic`) : le singleton des paramètres est la variable
  statique d'une fonction inline, que GCC émet comme symbole unique
  (STB_GNU_UNIQUE) et que le chargeur unifie entre bibliothèques, même
  chargées avec RTLD_LOCAL. Les bibliothèques `generic` et
  `generic-parallel` d'un même modèle partageaient leurs paramètres. Le
  code généré des modèles est maintenant dans un espace de noms anonyme ;
- bornes physiques des arguments des propriétés matériau (code amont,
  toutes les interfaces) : non exportées quand l'argument n'avait pas de
  bornes standard ;
- interpolation par spline cubique sur GPU (nvcc) : NaN entre et après les
  points de collocation, pour les propriétés matériau comme pour les
  écrouissages définis par des données (dès qu'il y a de la plasticité).
  Cause : nvcc 13.0 compile sans diagnostic, en fonction vide, toute
  fonction `constexpr` utilisant `__builtin_constant_p`, que libstdc++ 14
  appelle dans `std::advance` ; `std::advance`, `std::next` et `std::prev`
  ne font donc rien sur GPU (reproducteur
  `results/nvhpc-bugs/nvcc-builtin-constant-p.cu`). Fonctions internes
  dans `CubicSpline.ixx`, seul code device de TFEL concerné ;
- DigitSeparatorTest : la règle d'écoulement ne calculait pas `df_dseq`,
  que le Newton lisait sans initialisation. Le résultat dépendait des
  restes de la pile (écart en OpenMP à plusieurs fils, avec compare compilé
  en -O3). Démontré avec `-ftrivial-auto-var-init=zero` et `pattern`.

Défauts du programme de comparaison, corrigés : variables externes
vectorielles ou tensorielles données variable par variable (un champ par
composante au lieu d'un champ par variable). Nouveau : option
`--parameters` (paramètres changés dans les deux bibliothèques).

Couverture ajoutée à la suite `generic-parallel` :

- `TensorialVariablesTest` : tableaux de gradients et de forces, variables
  externes et internes de tous les types (scalaires, vecteurs, tenseurs
  symétriques, tenseurs, tableaux) ;
- paramètres changés à l'exécution pour un comportement et un modèle ;
- `CubicSplineHardeningTest` : écrouissage interpolé par une spline
  cubique, chargement plastique.

Chaque nouveau test a été vérifié dans les deux sens : il passe avec les
correctifs et échoue sans (ancien programme de comparaison, anciennes
tailles, paramètres fixés dans la seule référence).

Les écarts restants viennent des arrondis : `scripts/roundoff-triage.py`
rejoue sur l'hôte la comparaison la plus défavorable de chaque écart, la
bibliothèque du lot compilée avec et sans contraction FMA. Pour chacun,
aucun écart sans FMA, un écart du même ordre avec. Aux grandes amplitudes,
certaines lois sont à la limite de convergence de leur solveur local :
une différence d'arrondi y fait échouer un point d'un côté et pas de
l'autre (par exemple Chaboche2 à l'amplitude 1e-2 : 101 échecs sur 1009
dans la référence, 2 de plus sur GPU, 4 statuts différents sur l'hôte avec
FMA). Même sur l'hôte, g++ contracte différemment les deux bibliothèques
d'une loi : quelques écarts apparaissent aussi en OpenMP, et disparaissent
sans FMA.

Cas non concluants restants : échecs voulus (MaterialPropertiesBoundsCheck,
ExternalStateVariablesBoundsCheck, StateVariablesBoundsCheck,
InitializationFailureTest : même statut dans la référence et le lot), et
lois dont les valeurs synthétiques font diverger la référence (certaines
sont couvertes par la suite avec des valeurs physiques). Les tests des
avertissements de MFront WarningTest80 à 82 nomment `f` et `df_dseq` sans
les calculer : résultats indéfinis par construction.

Partage d'un GPU entre processus : PolyCrystals demande 113 à 121 Ko de
pile par fil, soit 14,8 à 15,9 Go réservés sur la RTX 5080 quel que soit le
nombre de points. Un autre processus (même son seul contexte CUDA) suffit
à manquer de mémoire : l'interface renvoie proprement -3 « out of
memory ». Six processus simultanés sur une loi ordinaire (pile de 3,4 Ko,
100 003 points chacun) : aucune erreur. Question ouverte pour MGIS : 7
« illegal memory access » sur environ 1000 lancements, uniquement quand 6
processus réservaient chacun environ 15 Go de pile sur 16 Go ; jamais
seul, jamais sous compute-sanitizer. Une loi dont la pile occupe presque
tout le GPU doit l'avoir pour elle seule.

Résultats :

- passe d'exécution finale, amplitudes 1e-4, 1e-3 et 1e-2 (105 min) :
  openmp-host-gcc 467 OK sur 485, cuda-nvcc 401 OK sur 485 ; tous les
  autres résultats expliqués (76 écarts d'arrondi confirmés par le tri,
  72 sur GPU et 4 sur l'hôte,
  WarningTest80 à 82, 10 cas non concluants identiques dans les deux
  configurations, un délai dépassé à l'amplitude 1e-2, un manque de
  mémoire GPU passager) ;
- propriétés matériau : 56/56 sur openmp-host-gcc et cuda-nvcc ;
- suite `generic-parallel` à -O3 : openmp-host-gcc, stlpar-parunseq-gcc
  et kokkos-host 531/531 (avant l'ajout de CubicSplineHardeningTest),
  cuda-nvcc 540/540 (avec) ;
- balayages de compilation des 608 fichiers avec AdaptiveCpp (28 min),
  nvc++ -stdpar=gpu (329 min) et nvc++ -mp=gpu (335 min) : 541 OK, 60
  fichiers invalides et 7 fragments, comme nvcc et clang, fichier par
  fichier. Avec nvc++ -mp=gpu, les polycristaux de 30 grains
  (BerveillerZaouiPolyCrystal, PolyCrystals) demandent environ 4 h et
  5 Go au front-end : `compile-sweep.sh` a maintenant un délai maximal
  par loi (`COMPILE_TIMEOUT`, 6 h par défaut, statut TIMEOUT).

### 2026-10-07 (matin), grandes validations à -O3 (R4, première partie)

TFEL de la branche (tous les correctifs de R3) reconfiguré, installé et
validé par `scratch/validate3.sh` :

- ctest complet de TFEL : 6965 tests sur 6965 ;
- matrice des tests `generic-parallel` à -O3
  (`results/generic-parallel-tests-2026-10-07-r4-O3.txt`) : 15
  configurations de lois (cuda-clang, cuda-nvcc, hip-clang-nvidia-platform,
  hip-hipcc, kokkos, kokkos-host, openmp-host-gcc,
  openmp-offload-nvhpc-gpu, stlpar-parunseq-clang, stlpar-parunseq-gcc,
  stlpar-parunseq-icpx, stlpar-parunseq-nvhpc-gpu, sycl-default-acpp,
  sycl-default-acpp-omp, sycl-default-icpx) et 16 de propriétés matériau,
  toutes réussies, TDLS compris ;
- balayages de compilation nvcc et clang : 541 OK, 60 fichiers invalides
  et 7 fragments, identiques fichier par fichier au balayage de référence.

Reste de R4 : passage à l'environnement spack mmm-gpu (tout y est
installé sauf MFEM, reporté), puis nouvelle validation.

### 2026-10-07 (après-midi), passage à l'environnement spack mmm-gpu (R4, en cours)

`scripts/env.sh` choisit l'environnement selon `PROTO_ENVIRONMENT` :
`tfelgpu_local` (défaut, validé) ou `mmm-gpu` (deps/spack/mmm-gpu :
Kokkos 5.2.2, AdaptiveCpp 25.10, HIP 7.2.3 et LLVM 20.1.8 de spack, CUDA
13.0.2 du HPC SDK pour nvcc et HIP). TFEL est construit dans
`build/tfel-mmm-gpu` et installé dans `install/tfel-mmm-gpu` ;
`build-tfel.sh`, `scripts/run-generic-parallel-tests.sh`,
`scratch/dev-test.sh` et `scratch/validate3.sh` suivent cette variable. Vue
de mmm-gpu régénérée (`spack env view regenerate` : `acpp` manquait). MFEM
reste reporté.

Validation dans mmm-gpu (interrompue à 14:40, arrêt de la machine) :

- ctest complet de TFEL : 6965 tests sur 6965 ;
- matrice à -O3 (`results/generic-parallel-tests-2026-10-07-mmm-gpu-O3.txt`,
  partielle) : cuda-nvcc, kokkos, kokkos-host, openmp-host-gcc,
  openmp-offload-nvhpc-gpu, stlpar-parunseq-clang, stlpar-parunseq-gcc et
  stlpar-parunseq-icpx réussis ;
- cuda-clang et hip-clang : `libcudart.so.12` introuvable (clang 20 compile
  avec CUDA 12.9). Corrigé dans `scripts/env.sh` (bibliothèques de CUDA
  12.9 après celles de CUDA 13) ;
- hip-hipcc : en-têtes de HIP 5.7 du système pris à la place de ceux de
  HIP 7.2 (hipcc ne donne pas les en-têtes de HIP à nvcc sur la plateforme
  NVIDIA). Défaut de TFEL corrigé (`cmake/modules/tfel.cmake`, patch 0066) :
  dossier des en-têtes de HIP ajouté aux options de hipcc s'il n'est pas un
  dossier implicite du compilateur C++ (inchangé pour HIP dans /usr).

Reprise :

```bash
PROTO_ENVIRONMENT=mmm-gpu ./build-tfel.sh install -j12   # reconfiguration (patch 0066)
PROTO_ENVIRONMENT=mmm-gpu scratch/validate3.sh 2026-10-07-mmm-gpu
```

Puis, si tout passe : mmm-gpu par défaut dans `scripts/env.sh`, R4 cochée.
Ensuite : campagne de performance (grands lots, CPU multicœur contre GPU,
plusieurs familles de lois), clang-tidy et relecture indépendante, puis D1
(MGIS).

### 2026-10-07 (soir), environnement mmm-gpu validé (fin de R4)

Validation complète dans l'environnement spack mmm-gpu, après les
correctifs de l'après-midi (bibliothèques de CUDA 12.9 dans
`LD_LIBRARY_PATH`, en-têtes de HIP donnés à hipcc) :

- ctest complet de TFEL : 6965 tests sur 6965 ;
- matrice des tests `generic-parallel` à -O3
  (`results/generic-parallel-tests-2026-10-07-mmm-gpu-O3.txt`) : 15
  configurations de lois et 16 de propriétés matériau, toutes réussies,
  dont hip-hipcc et hip-clang avec HIP 7.2.3, kokkos avec Kokkos 5.2.2 et
  sycl-default-acpp avec l'AdaptiveCpp de spack ;
- balayages de compilation nvcc et clang : 541 OK, identiques fichier par
  fichier au balayage de référence.

Correctif amont supplémentaire de TFEL (patch 0067) : argument du `endif`
de la configuration hip-clang dans `tfel.cmake` (avertissement de CMake à
chaque configuration). Remarque : les substitutions de TFEL sont dans le
cache de CMake ; une reconfiguration qui en change une valeur demande
`cmake -U "tfel-substitutions*" .` (fait par `scratch/validate3.sh`).

mmm-gpu est maintenant l'environnement par défaut de `scripts/env.sh` et
de `build-tfel.sh` ; `PROTO_ENVIRONMENT=tfelgpu_local` garde l'ancien
(validé lui aussi le 07/10 au matin). TFEL est donc validé dans les deux
environnements.

Suite proposée : campagne de performance (grands lots, CPU multicœur
contre GPU, plusieurs familles de lois), clang-tidy et relecture
indépendante, puis D1 (MGIS).

### 2026-10-08 (nuit), clang-tidy et relecture indépendante (R5)

clang-tidy 20 (vérifications bugprone, clang-analyzer, performance et
quelques modernize, readability et cppcoreguidelines) sur les sources de
l'interface dans MFront et sur le code généré de cinq lois (comportements
en petites et grandes déformations, modèle, propriétés matériau), avec les
en-têtes de l'interface. Aucun défaut dans les en-têtes exécutés au
calcul. Dans le générateur : variables mortes (code amont des propriétés
matériau) et deux copies inutiles, corrigées. Le reste est cosmétique.

Deux relectures indépendantes, l'une des en-têtes exécutés au calcul,
l'autre du générateur. Défauts confirmés et corrigés :

- lois de comportement et modèles :
  - un lot vide dont les champs sont nuls (données d'un `std::vector`
    vide, processus sans point d'un matériau) était refusé (`-2`) ;
  - des options de calcul invalides d'une loi en grandes déformations
    donnaient un échec en chaque point (`-1`), ce qui pousse un code
    éléments finis à réduire le pas sans fin. Elles sont maintenant
    refusées pour le lot (`-2`) ;
  - les fonctions d'initialisation et les post-traitements écrasaient les
    facteurs de pas `rdts` proposés par l'intégration ;
  - le statut de sortie n'était pas complet après une exception levée
    avant la vérification des données ;
  - SYCL, politique `automatic` : l'étendue copiée d'un champ utilisait
    les tailles maximales et non celles des options. Pour une loi en
    grandes déformations avec la contrainte de Cauchy, la copie lisait
    puis réécrivait des valeurs après la fin du tableau de l'appelant ;
  - SYCL, politique `automatic` : deux sorties rangées dans la même zone
    (contraintes et variables internes d'un point stockées à la suite)
    s'écrasaient au retour. Seules les valeurs d'un champ sont recopiées ;
  - SYCL : le champ d'un groupe donné variable par variable était copié
    alors que le contrat le dit ignoré ;
  - SYCL : aucun ordre entre le noyau et le travail déjà soumis à une file
    hors ordre. Les calculs asynchrones demandent maintenant une file
    ordonnée (`-2` sinon) ; en synchrone, une file hors ordre est
    synchronisée avant le noyau ;
- propriétés matériau (code amont) :
  - avec `disable_runtime_checks`, la variante avec pas ne compilait pour
    aucun backend (`mfront_errno_old` non déclaré) ;
  - SYCL : la mémoire du périphérique n'était jamais libérée (fuite à
    chaque appel) ; étendue des champs avec pas trop grande ; valeurs
    situées entre les valeurs d'une sortie avec pas écrasées au retour ;
  - le tableau des pointeurs des arguments, souvent sur la pile de
    l'appelant, était lu dans les noyaux (Kokkos, OpenMP, stdpar, SYCL).
    Les pointeurs sont maintenant copiés dans le noyau. Ce défaut était
    masqué par la carte de développement : le pilote y active HMM, qui
    rend toute la mémoire de l'hôte lisible par le GPU ;
  - transferts SYCL faits hors du bloc `try` ; garde d'inclusion sans le
    nom du backend ; bornes et valeurs par défaut écrites avec 6 chiffres ;
    un bloc CUDA de trop ; message d'erreur faux ; fichier parasite.

Nouveau défaut de nvc++ 25.11 trouvé en corrigeant : une lambda qui capture
un `std::array` de pointeurs, passée par `map(to : ...)` à une région
OpenMP, lit des pointeurs faux sur GPU (-O2, -O3), ou ne compile pas (-O0,
-O1). Reproducteur `results/nvhpc-bugs/omp-lambda-pointer-array.cxx`. Les
pointeurs sont capturés un par un.

Tests ajoutés : lot vide à champs nuls et options invalides (programme de
comparaison, tous les backends), compilation d'une propriété matériau avec
les vérifications désactivées. Chacun échoue avec les anciennes
bibliothèques et passe avec les nouvelles.

Documenté sans correction : les `@SharedData` sont copiées une fois, sur le
périphérique du premier appel (un processus à plusieurs GPU ne peut pas
partager la bibliothèque) ; hors CUDA, HIP, SYCL et Kokkos, le tas doit être
accessible par le périphérique.

Restent ouverts, de moindre gravité : symboles `generic_parallel_backend`
et version d'API absents pour les propriétés matériau, duplication entre
les interfaces comportement et modèle et entre les deux noyaux des
propriétés matériau, encodage de `rdt` écrit deux fois (C et C++),
`device_stack_size` ignoré sans message si le runtime CUDA est absent.
Les relectures sont dans `scratch/review-r5` (non versionné).

Validation à -O3 avec l'arbre de développement :

- comportements : openmp-host-gcc, cuda-nvcc et kokkos 601/601 (540
  comparaisons, 50 lots vides, 11 options invalides), sycl-default-acpp
  627/627 (RTX 5080, mémoire de l'hôte : copies `automatic` exercées) ;
- propriétés matériau : stlpar-parunseq-gcc, openmp-host-gcc, kokkos-host,
  kokkos, cuda-nvcc, hip-hipcc, sycl-default-acpp, sycl-gpu-acpp,
  stlpar-parunseq-nvhpc-gpu et openmp-offload-nvhpc-gpu réussis.

`scratch/dev-test.sh` produit maintenant tous les fichiers de
configuration nommés par un test (sycl-gpu-acpp utilise celui de
sycl-default-acpp).

Limite de la validation sur cette machine : le pilote NVIDIA y active HMM,
qui rend toute la mémoire de l'hôte accessible au GPU. Un accès du GPU à
la pile ou au tas de l'hôte y fonctionne, alors qu'il plante sur une
machine sans HMM. Pour s'en assurer, il faudrait relancer les tests GPU
avec le module `nvidia_uvm` chargé avec `uvm_disable_hmm=1` (droits
administrateur).

### 2026-10-08 (nuit, suite), données partagées, symboles, duplication (R6)

Trois points laissés ouverts par R5.

- `@SharedData` : une copie par GPU (GPU courant en CUDA et HIP), par
  périphérique et par contexte en SYCL. Kokkos garde une copie (un GPU par
  processus). La copie est rangée sous l'adresse de la fonction qui fournit
  la donnée (`getSharedDataCopy<&Behaviour::mfront_get_shared_data_x>`,
  cache `SharedDataCopies<T, Device>`). L'intégration, les fonctions
  d'initialisation et les post-traitements partagent donc la même copie.
  Coût : une requête du GPU courant par appel. Nouveau test SYCL : la loi
  polycristalline est intégrée sur le GPU, puis sur le CPU.
- propriétés matériau : symboles `<nom>_generic_parallel_backend` et
  `<nom>_generic_parallel_api_version`, comme pour les lois. La version est
  définie dans le nouvel en-tête C
  `MFront/GenericParallelMaterialProperty/MaterialProperty.h`. Les tests de
  tous les backends vérifient ces symboles.
- duplication : un seul code des noyaux des propriétés matériau pour tous
  les backends (1385 lignes avant, 1013 après). Les deux copies avaient
  divergé : ordre des vérifications, messages, restauration de `errno`.
  Les lois et les modèles partagent le début et la fin des fonctions
  générées, les noms de fichiers et de bibliothèque, et le choix du
  backend. Un backend inconnu donne la liste des backends disponibles. Côté
  hôte, le facteur de pas est encodé par les fonctions C du contrat.

Le code généré de toutes les lois de test a été comparé avant et après,
sur 6 backends : seuls les changements voulus apparaissent.

Validation complète à -O3, avec HMM désactivé (`uvm_disable_hmm=Y`) :

- installation, tests de TFEL : 6965/6965 ;
- lois : 15 configurations sur 15 ;
- propriétés matériau : 14 sur 16. cuda-nvcc et cuda-clang échouaient à
  cause du test amont `test-cuda.cu` : le second appel de `test9` écrivait
  le résultat dans la mémoire de l'hôte depuis un noyau (compute-sanitizer :
  écriture invalide). HMM masquait ce défaut. Le test utilise maintenant
  la mémoire du périphérique : cuda-nvcc et cuda-clang réussis,
  compute-sanitizer sans erreur. Le test HIP utilise déjà de la mémoire
  managée ;
- balayages de compilation : nvcc 541 OK, 60 INVALID, 7 FRAGMENT, comme
  avant ; clang identique (541 OK).

Sans HMM, aucun accès caché du GPU à la mémoire de l'hôte n'apparaît dans
le code de l'interface.

### 2026-10-08 (nuit, fin), paramètres des variables de comportement, rotations (R7)

Trois limites restantes.

- Paramètres chaîne des modèles : rien à corriger pour le GPU. Un modèle
  avec un paramètre chaîne échoue déjà à l'analyse sur CPU, avec
  l'interface `generic`. Le reproducteur `reproducers/device-strings`
  montre au passage que `std::string_view` ne fonctionne pas dans le code
  GPU : résultats faux sans message avec nvcc et hipcc, erreurs avec nvc++
  (appels à `std::char_traits`). Un petit type sans `char_traits` fonctionne
  avec les six compilateurs.
- Paramètres des variables de comportement : chaque classe compatible GPU
  a une structure `mfront_parameters_values`, avec ses paramètres et ceux
  de ses variables de comportement (récursivement). Elle est remplie sur
  CPU depuis les initialiseurs (`mfront_get_parameters_values`), puis
  appliquée sur le périphérique (`mfront_set_parameters_values`). Les
  fabriques de variables de comportement la reçoivent aussi. La structure
  des paramètres de l'interface `generic-parallel` en dérive : plus de
  types dupliqués, plus de refus des valeurs changées. Nouveau test : le
  paramètre `hardening_factor` de `BehaviourVariableTestPlasticity` est
  changé par fichier, puis les résultats CPU et GPU sont comparés. Une
  version qui ignore la valeur échoue sur GPU (écart de 1,5e-3).
- Rotation de l'opérateur tangent : tous les blocs sont tournés, sur CPU
  (`generic`) comme sur GPU. Les blocs de tableaux sont tournés sous-bloc
  par sous-bloc, car MFront range les sous-blocs l'un après l'autre. Les
  paires sans type TFEL (vecteur et tenseur, tenseur et tenseur
  symétrique : TFEL n'a pas de `change_basis` pour `st2tot2`) passent par
  un produit générique `Q1 K Q2ᵀ`, utilisable sur GPU
  (`GenericBehaviour/RotateTangentOperatorBlock.hxx`). Le code généré des
  cas existants ne change pas. Nouvelle loi de test
  `OrthotropicTangentOperatorBlocksTest` : 9 blocs de tous les types.

Vérifications limitées au code touché :

- reproducteur `reproducers/rotations`, indépendant du code qui tourne les
  blocs : `rotateTangentOperatorBlocks(K) x =
  rotateThermodynamicForces(K rotateGradients(x))` sur les 5 hypothèses
  (écart de l'ordre de 1e-15), rotation sur place identique. Une version
  volontairement fausse est détectée ;
- comparaisons ciblées des lois orthotropes (`scratch/targeted-test.sh`),
  CPU et GPU : 137 OK ;
- suites complètes des lois, qui contiennent aussi les comparaisons des
  rotations de la nouvelle loi : openmp-host-gcc et
  cuda-nvcc 637/637, sycl-default-acpp 664/664, stlpar-parunseq-nvhpc-gpu
  637/637 avec le contournement ci-dessous.

Nouveau défaut de nvc++ 25.11 (`-stdpar=gpu`), révélé par la suite
stlpar-parunseq-nvhpc-gpu : l'édition des liens du code device de la
bibliothèque des lois échoue (« nvlink error : Multiple definition of
'__staticinit_110045_…' », entre `HosfordPerfectPlasticity` et
`InitializeFunctionTest`). nvc++ range certaines constantes locales du code
device (ici la lambda vide de `assign` et `addAndAssign` des tableaux de
TFEL) dans des données globales visibles, nommées par un compteur propre à
chaque fichier. Deux fichiers qui donnent le même compteur à la même
fonction inline entrent en conflit. Le code des paramètres n'a fait que
décaler les compteurs : la validation R6 passait par chance. Sans l'en-tête
des rotations, le conflit reste (compteurs décalés de 8). Le reproducteur
`results/nvhpc-bugs/stdpar-staticinit-link.cxx` (deux copies d'un fichier de
quatre lignes) échoue de -O0 à -O3. Contournement : code device de chaque
fichier lié seul (`-gpu=nordc` dans `config-stlpar-parunseq-nvhpc-gpu.json`),
possible car les sources générées n'appellent pas de fonction device
définie dans un autre fichier. OpenMP déporté n'émet pas ces symboles.

### 2026-10-08 (matin), gestion des erreurs homogène (R8)

Principe, le même partout : le chemin sans erreur ne coûte rien. Les
noyaux n'écrivent dans les statistiques ou les statuts que si un événement
survient en un point (échec, résultat peu fiable, baisse du facteur de pas,
borne violée, résultat invalide). Pas de chaîne de caractères sur le GPU :
le message détaillé d'un point s'obtient en le rejouant avec l'interface
`generic`.

État avant (relevé dans `scratch/r8-map.md`, non versionné) :

- lois, chemin synchrone CUDA/HIP : allocation, deux copies et libération
  à chaque appel ; SYCL : allocation, trois attentes et libération à
  chaque appel ;
- propriétés matériau : tableau de statuts alloué à chaque appel par tous
  les backends (`cudaMallocManaged` puis `cudaFree`, qui synchronise le
  GPU ; `std::vector`, donc `cudaMallocManaged` avec nvc++ ;
  `Kokkos::View` partagée ; vecteur USM partagé), puis
  `cudaDeviceSynchronize` ;
- aucun indice du point en échec : il fallait parcourir le tableau des
  statuts par point.

Après :

- contrat C des lois en version 2 : `first_failure_index` dans
  `mfront_gpb_BatchStatistics` et `mfront_gpb_OutputStatus`. C'est le plus
  petit indice en échec, quel que soit l'ordre des calculs (minimum
  atomique ou réduction). Le message donne le nombre d'échecs et ce point ;
- CUDA/HIP : un vote par warp ; les warps sans événement ne font plus
  aucune réduction ni opération atomique ;
- nouvel en-tête partagé `MFront/GenericParallel/StatusBlockPool.hxx` :
  réserve de blocs accessibles par l'hôte et le GPU, alloués au premier
  appel sur un GPU puis réutilisés, un bloc par appel en cours (appels
  concurrents possibles). Il sert au chemin synchrone des lois (CUDA/HIP :
  mémoire managée, SYCL : USM partagée, une seule attente) et aux statuts
  des propriétés matériau de tous les backends (managée, tas, USM
  partagée, `Kokkos::SharedSpace` libérée au `finalize`). L'hôte
  initialise le bloc et le lit après la synchronisation : sans événement,
  aucune donnée ne passe entre l'hôte et le GPU ;
- propriétés matériau CUDA/HIP : synchronisation du flux au lieu de tout
  le GPU ;
- Kokkos synchrone, OpenMP et algorithmes parallèles : réductions
  inchangées, avec l'indice du premier échec en plus.

Défaut de test révélé : le test HIP des propriétés matériau supposait la
mémoire de `hipMallocManaged` nulle. Les allocations libérées à chaque
appel le masquaient. La mémoire est maintenant mise à zéro, comme les
`thrust::device_vector` du test CUDA.

Le contrat des propriétés matériau (`mfront_gmp_OutputStatus`, commun avec
l'interface `generic`) n'a pas de champ pour un indice de point : elles
disent quelle borne est violée, pas où.

Mesure (`reproducers/sync-call-overhead`, RTX 5080, loi élastique, données
en mémoire managée) : un appel synchrone CUDA passe de 135 à 29 µs pour
1000 points, et de 546 à 442 µs pour 100000 points. Le coût fixe par appel
baisse d'environ 105 µs.

Propriétés matériau, dans la foulée : avec la politique `Warning`, une borne
violée donne maintenant le statut 1 (résultat calculé, à utiliser avec
précaution), comme l'interface `generic`. Le statut restait nul. Les codes
propres à l'interface (-2, -5, -6, -7) sont documentés. Nouveau test
`test14` (tests C++ et CUDA).

Vérifications :

- tests ciblés des lois avec échecs (`ImplicitNortonWithEnergies` en
  politique stricte, modèle `GenericParallelModelTest`, fonctions
  d'initialisation, post-traitements) : 50/50 sur cuda-nvcc, hip-hipcc,
  sycl-default-acpp, kokkos, openmp-host-gcc, stlpar-parunseq-gcc. Une
  version sans le minimum atomique de l'indice est détectée (12 échecs) ;
- propriétés matériau : cuda-nvcc, hip-hipcc, stlpar-parunseq-gcc,
  openmp-host-gcc, kokkos, sycl-default-acpp ;
- validation complète sans HMM (`validate3.sh 2026-10-08-r8`) : tests de
  TFEL 6965/6965, matrice -O3 31/31 (stlpar-parunseq-nvhpc-gpu compris,
  avec `-gpu=nordc`), balayages nvcc et clang 541 OK, comme avant.

`scratch/dev-test.sh` saute en mode ciblé les comparaisons appelées
directement par un test (lois non construites).

### 2026-10-08 (matin, suite), MGIS : chargement et calcul de lots (D1, étapes 1 à 3)

Branche `gpu-proto` de `mgis/` (base `818f0f9`). Les patchs sont dans
`patches/mgis`, régénérés par `scripts/sync-mgis.sh`. MGIS est construit
dans `build/mgis` contre le TFEL du proto. Le TFEL de spack doit être
écarté : `CMAKE_PREFIX_PATH` le trouvait en premier.

- Étape 1, chargement. Une bibliothèque `generic-parallel` était chargée
  sans erreur et ses fonctions de lot appelées avec les arguments d'un
  point (comportement indéfini). `BehaviourDescription` lit maintenant
  l'interface, le backend et la version des structures (2). `Behaviour`
  porte les pointeurs de lot (`batch_b`, `batch_f`, trois rotations).
  `LibrariesManager` refuse les fonctions d'un point pour une bibliothèque
  de lots, et inversement. MGIS garde sa copie des structures C
  (`BatchDataView.hxx`), vérifiée contre `BatchData.h` par un test de
  disposition mémoire.
- Étape 2, intégration. `integrate(BatchExecutionPolicy, MaterialDataManager,
  options, dt[, b, e])` construit la vue du lot sur les tableaux du
  gestionnaire, sans copie. Les propriétés matériau et les variables
  externes sont données variable par variable : valeur uniforme ou valeur
  par point. Les fonctions existantes basculent seules, avec la politique
  du gestionnaire (`batch_execution_policy`). Les statuts -2 et -3 sont
  levés comme exceptions, jamais pris pour un échec de convergence. Un
  échec donne l'indice du premier point (`n`). `integrate_debug` est
  refusé : la bibliothèque `generic` sert au diagnostic.
- Étape 3 : fonctions d'initialisation, post-traitements et rotations de
  lots, avec le même basculement (pool de fils compris).

GPU : avec le backend `cuda`, l'appelant fournit toute la mémoire en
mémoire managée (initialiseurs de MGIS, `std::span`). Les valeurs uniformes
gardées par MGIS sont en mémoire hôte : elles sont refusées avec un message
clair pour `cuda` et `hip`. Test `GenericParallelCUDATest` (option
`enable-generic-parallel-cuda-tests`) : Norton sur 10 pas et 1000 points,
élasticité avec module d'Young par point, HMM désactivé. Résultats égaux
à ceux du calcul par point aux arrondis près, compute-sanitizer sans
erreur, noyaux vus par Nsight (0,5 ms pour 1000 points).

Défauts de MGIS trouvés en chemin, corrigés dans des commits séparés :

- les axes des matrices de rotation étaient refusés quand leur nombre de
  valeurs était un multiple de la dimension (test inversé) : aucune
  `RotationMatrix2D` ni `RotationMatrix3D` ne pouvait être construite. Le
  nombre de points était aussi comparé à un reste au lieu d'un quotient ;
- plantage intermittent en fin de processus : les bibliothèques étaient
  déchargées alors que les fils de TBB (algorithmes parallèles de
  libstdc++) vivaient encore. Elles ne sont plus déchargées
  (`RTLD_NODELETE`) ;
- une propriété matériau ne pouvait pas être donnée par un `std::span`
  d'une seule valeur, alors que l'évaluation le traitait déjà comme
  uniforme. C'est la seule façon de placer une valeur uniforme en mémoire
  accessible par le GPU.

Noté sans correction : `executeInitializeFunction` par point ne met pas
`exit_status` à -1 en cas d'échec (seul `n` est rempli). Le chemin des lots
signale l'échec.

Vérifications : 44 tests de MGIS (41 d'origine, chargement, comparaison
lots/points sur CPU, comparaison sur GPU), chaque cas de comparaison
vérifié par une version volontairement fausse (décalage de plage,
rotations uniformes forcées, suffixe de mesure de contrainte faux). Une
loi orthotrope en grandes déformations a été ajoutée aux tests de MGIS
pour les rotations avec mesure de contrainte et opérateur tangent.

Étape 4a : l'emplacement des tableaux (`MemoryLocation` : hôte, GPU,
partagé) est déclaré par l'initialiseur du gestionnaire et vérifié. Un
backend `cuda` ou `hip` refuse des données en mémoire hôte. Hors de la
mémoire hôte, tous les tableaux des états doivent être donnés, ainsi que
ceux des opérateurs tangents et de la vitesse du son s'ils servent : MGIS
les allouerait en mémoire hôte. `update` et `revert`, qui copient sur
l'hôte, refusent des données sur le GPU. Chaque refus est testé par son
message. L'en-tête `BatchExecutionPolicy.hxx` n'était pas installé :
corrigé.

Reste du plan : ressource mémoire (étape 4b, qui rendra automatique ce que
le test GPU fait à la main), asynchrone et `update`/`revert` sur le GPU
(étape 5), propriétés matériau `generic-parallel` et documentation.

### 2026-10-08 (après-midi), espace mémoire (R9) et relecture de MGIS (D1, étape 4b)

TFEL (R9), deux commits plus un de dates :

- chaque bibliothèque exporte l'espace mémoire de ses champs,
  `<f>_generic_parallel_memory_space` : `host`, `device` ou `any`. La
  valeur est fixée à la compilation. `cuda` et `hip` : `device`. `stlpar`
  et `openmp` sur l'hôte : `host`. nvc++ (`-stdpar=gpu`, `-mp=gpu`) :
  `device`, ou `any` en mémoire unifiée. `kokkos` : selon l'accès de
  l'espace d'exécution à la mémoire hôte. `sycl` : `any` avec les
  transferts automatiques, `device` sinon. En-tête
  `MFront/GenericParallel/MemorySpace.hxx`. Les propriétés matériau
  exportent le même symbole. `compare` et les tests des propriétés
  vérifient sa valeur ;
- les configurations `stlpar` liées à TBB (gcc, clang, icpx) ajoutent
  `-Wl,-z,nodelete` : la bibliothèque reste chargée jusqu'à la fin du
  processus. Les fils de TBB tournaient encore quand elle était déchargée
  à la sortie. MGIS n'a plus besoin de `RTLD_NODELETE` pour toutes les
  bibliothèques.

MGIS : relecture indépendante de la branche (19 points) et trois points
relevés avant. Corrigé :

- la liste des backends écrite dans MGIS (`cuda`, `hip`) est remplacée par
  le symbole de TFEL : `isMemoryLocationAccessible(comportement, lieu)`.
  Les cas manqués avant (Kokkos sur GPU, OpenMP déporté, SYCL sans
  transfert, nvc++ sans mémoire unifiée, données sur GPU avec un backend
  hôte) sont couverts ;
- les fonctions d'un point (intégration, `integrate_debug`, initialisation,
  post-traitement, rotations, liaisons C) appelaient un pointeur nul pour
  une bibliothèque de lots. Elles lèvent une erreur explicite ;
- les valeurs gardées par MGIS (valeur uniforme, `LOCAL_STORAGE`) sont
  refusées quand la mémoire hôte n'est pas accessible aux calculs. Avant,
  seul le scalaire était vérifié : un `std::vector` passait et le noyau
  CUDA plantait, contexte CUDA perdu ;
- la masse volumique n'est passée que pour la vitesse du son. Elle peut
  être donnée par un `std::span` d'une valeur. Une propriété uniforme en
  `LOCAL_STORAGE` est stockée comme un scalaire ;
- tableaux de variables : la relecture les croyait refusés. MFront exporte
  chaque élément comme une variable (`v2_esv[0]`, `v2_esv[1]`) et refuse
  les types tableau pour les propriétés matériau et les variables
  externes. Un champ par variable, comme TFEL : le découpage par élément
  ajouté puis testé était faux pour un type tableau, il est retiré ;
- initialisation, post-traitements et rotations de lots prennent une
  politique d'exécution (file), comme l'intégration. Les rotations `sycl`
  échouaient toujours : file nulle. Les fonctions de rotation de lots sont
  chargées par des accesseurs typés, et les suffixes de grandes
  déformations écrits une fois (une valeur inconnue lève une erreur au lieu
  de prendre PK1) ;
- les gestionnaires d'état connaissent l'emplacement de leurs tableaux.
  `update`, `revert`, `updateValues`, l'extraction, `save`, `restore` et
  le calcul point par point refusent des tableaux sur le GPU. Un
  gestionnaire sans point, ou un modèle sans opérateur tangent, est
  accepté hors de la mémoire hôte ;
- un échec donne l'indice absolu du premier point dans le message, et un
  statut inattendu lève une exception. Les fonctions de statistiques
  asynchrones inutilisées sont retirées : elles reviendront avec l'étape 5 ;
- CMake : architecture CUDA `native` seulement avec CMake 3.24 ; sources
  des lots régénérées quand `mfront` change (les anciennes sources, sans le
  nouveau symbole, restaient).

Étape 4b : l'initialiseur du gestionnaire prend une ressource mémoire
(`std::pmr::memory_resource*`). MGIS alloue avec elle les tableaux non
donnés. Leurs valeurs sont initialisées par l'hôte : la ressource doit
allouer de la mémoire accessible par l'hôte (managée ou partagée). Sur le
GPU seul, tous les tableaux restent à donner. Ressources fournies en
en-têtes, compilées par le code appelant : mémoire managée CUDA et HIP,
mémoire partagée SYCL et Kokkos. Elles vérifient l'alignement demandé. Coût :
un pointeur par tableau, un appel virtuel par allocation, rien par calcul.

Tests ajoutés : plage avec gradients modifiés, post-traitement sur une
plage, échecs partiels donnés par une variable externe hors bornes
physiques (indice absolu), toutes les variables de
`TensorialExternalStateVariableTest` vérifiées valeur par valeur, grandes
déformations avec quatre couples mesure de contrainte et opérateur
tangent, vitesse du son, rotations 2D et 3D par point, fonctions avec
politique, ressource mémoire et emplacements refusés. Sur GPU : ressource
mémoire managée sans aucun tableau donné, puis intégration,
initialisation, post-traitement et rotation sur un flux CUDA qui n'est pas
le flux par défaut. Deux versions volontairement fausses (pas des points
des variables, indice relatif) sont détectées. Les 44 tests de MGIS
passent, HMM désactivé. Les en-têtes HIP, SYCL et Kokkos compilent.

Une sauvegarde plus ancienne de `backups/` a été écrasée par un
renommage. Ses commits restent dans `tfel-dev` (branches et reflog).

Reste pour MGIS : asynchrone et `update`/`revert` sur le GPU (étape 5,
après mfem-mgis), propriétés matériau `generic-parallel` (étape 6),
documentation (étape 7). Statuts par point non exposés. Les rotations
données par `RotationMatrix2D`/`3D` construisent les matrices sur l'hôte à
chaque appel : pour des rotations répétées, les donner en tableau.

### 2026-10-08 (après-midi, suite), MGIS : propriétés matériau et modèles de lots (D1, étape 6)

Tout ce que l'interface `generic-parallel` de TFEL produit est maintenant
utilisable par MGIS, comme pour les lois.

- Propriétés matériau. Au chargement, MGIS lit l'interface, la version des
  fonctions (1), le backend et l'espace mémoire. Il charge la fonction
  d'un ensemble de points (`batch_fct`) au lieu de celle d'un point. Avant,
  une bibliothèque `generic-parallel` donnait un pointeur de mauvaise
  signature. Chaque fonction refuse l'autre type de bibliothèque.
- `evaluate(valeurs, propriété, entrées[, emplacement, politique])` évalue
  une propriété sur un ensemble de points. Les entrées sont uniformes ou
  données par point. Un seul appel pour une propriété `generic-parallel`,
  calcul point par point sur l'hôte sinon : même basculement que
  `integrate`. L'emplacement est vérifié contre l'espace mémoire. Les
  statuts gardent le sens de l'interface `generic` (0, 1 hors bornes avec
  la politique Warning, négatif en cas d'échec).
- Modèles : l'interface `generic-parallel` des modèles exporte les mêmes
  symboles que les lois. MGIS les chargeait déjà par `load` : un test le
  vérifie, avec un paramètre modifié par `setParameter`.
- `MemoryLocation` passe dans l'espace de noms `mgis`
  (`MGIS/MemoryLocation.hxx`) : la règle d'accès est écrite une fois pour
  les lois et les propriétés. Les fonctions CMake des bibliothèques de lots
  (lois, modèles, propriétés) partagent leur code
  (`cmake/modules/generic-parallel.cmake`).

Asymétrie qui vient de TFEL : les fonctions des propriétés matériau ne
prennent pas d'options d'exécution (file, statistiques asynchrones). Elles
calculent sur la file par défaut de leur backend, de façon synchrone.
`evaluate` n'a donc pas de surcharge avec politique. L'ajouter demanderait
un changement du contrat des propriétés matériau dans TFEL (version 2).

Tests (minimaux) : propriétés sur l'hôte (`stlpar`) et sur le GPU
(`cuda`, mémoire managée), valeurs par point et uniformes, bornes
physiques, statuts des trois politiques de bornes, emplacements, arguments
invalides ; version et disposition du statut vérifiées contre TFEL ;
modèle `ode_rk54` en lot sur l'hôte et sur le GPU. Deux versions
volontairement fausses de `evaluate` et un paramètre ignoré sont détectés.
45 tests de MGIS, tous verts, HMM désactivé.

Le message `Error at line ... of file 'src/targets.lst'` pendant la
construction vient d'appels parallèles de `mfront` dans un même dossier.
Il existe aussi dans la construction de MGIS d'origine. Il est sans effet
sur les sources.

Reste pour MGIS : asynchrone et `update`/`revert` sur le GPU (étape 5,
après mfem-mgis), documentation (étape 7).

### 2026-10-08 (soir), options d'exécution des propriétés matériau (R10)

Les fonctions des propriétés matériau `generic-parallel` n'avaient pas
d'options d'exécution. Le backend SYCL ajoutait la file comme argument
spécifique, au milieu des autres : la signature dépendait du backend.
`evaluate` de MGIS, écrit pour la signature commune, aurait appelé une
bibliothèque SYCL avec des arguments décalés. CUDA lançait sur le flux par
défaut, Kokkos sur l'instance par défaut avec un `Kokkos::fence()` global.

TFEL, version 2 des fonctions (`MFRONT_GPMP_API_VERSION`) :

- un dernier argument `const mfront_gpmp_ExecutionOptions*` (nul pour les
  options par défaut) : file et préservation de la taille de pile, même
  sens que pour les lois. Flux CUDA et HIP, `sycl::queue*` obligatoire
  (statut -5 sinon), instance de l'espace d'exécution Kokkos (attente de
  cette instance seulement). Les algorithmes parallèles et OpenMP n'ont
  pas de file ;
- l'argument spécifique de SYCL et le mécanisme d'arguments spécifiques
  disparaissent : une seule signature pour tous les backends ;
- le code de la taille de pile (`stlpar`, OpenMP déporté, SYCL) est écrit
  une fois et respecte `preserve_stack_size` ;
- l'évaluation reste synchrone. L'asynchrone demanderait un autre compte
  rendu des bornes (le message nomme la variable, il est construit après
  le calcul) : non fait.

Tests : les six fichiers de test passent les options. Un flux CUDA et HIP
non bloquant est vérifié par un noyau qui écrit les arguments sur ce flux
après un délai, sans attente de l'hôte : sur le flux par défaut, le test
échoue (vérifié). Instance Kokkos (`kokkos-host` reçoit l'espace hôte de sa
configuration), file SYCL absente refusée. Dix configurations, une par
backend : `stlpar-parunseq-gcc`, `openmp-host-gcc`, `kokkos-host`,
`kokkos`, `cuda-nvcc`, `hip-hipcc`, `sycl-default-acpp`, `sycl-gpu-acpp`,
`stlpar-parunseq-nvhpc-gpu`, `openmp-offload-nvhpc-gpu`. Toutes passent.
La matrice et la suite complète de TFEL ne sont pas relancées.

MGIS : `BatchExecutionPolicy` passe dans l'espace de noms `mgis`, partagé
par les lois, les modèles et les propriétés. `evaluate(politique, ...)`
existe comme `integrate(politique, ...)`. Les statuts d'arguments
invalides et d'erreur du backend (-5, -6, -7) sont levés comme
exceptions, comme pour les lois. 45 tests de MGIS, tous verts, dont
l'évaluation sur un flux CUDA non par défaut.

Lois, modèles et propriétés matériau ont maintenant la même interface
d'exécution dans TFEL et dans MGIS. Seule différence : les lois savent
calculer en asynchrone dans TFEL (statistiques sur le GPU), pas les
propriétés matériau. MGIS n'utilise l'asynchrone pour aucun des deux
(étape 5).

### 2026-10-08 (soir, suite), évaluation asynchrone des propriétés matériau (R10, suite)

Dans TFEL, les lois savaient calculer en asynchrone (statistiques données
par les options d'exécution), pas les propriétés matériau. C'est réglé :
synchrone ou asynchrone selon le choix de l'appelant, sur tous les
backends, pour les deux.

- Les options d'exécution des propriétés prennent les statuts d'une
  évaluation asynchrone (`statuses`, tableau de
  `<name>_generic_parallel_number_of_statuses` entiers mis à zéro par
  l'appelant). La fonction rend alors la main dès le calcul soumis à la
  file (CUDA, HIP, Kokkos, SYCL si aucune valeur n'est copiée). Le statut
  de sortie ne rapporte que les arguments invalides et les erreurs du
  backend.
- L'analyse des statuts est générée une fois, dans une fonction exportée
  (`<name>_generic_parallel_report`). L'évaluation synchrone l'appelle
  elle-même ; après une évaluation asynchrone, l'appelant l'appelle une
  fois la file synchronisée. Le compte rendu est le même dans les deux cas,
  message nommant la variable compris.
- Les noyaux ne changent pas : ils n'écrivent les statuts qu'en cas de
  borne violée ou de résultat invalide. Coût nul en cas de succès, comme
  avant.
- La structure n'est publiée que dans le proto : la version reste 2.

Tests limités à ce qui change, une configuration par backend (dix, toutes
vertes) : évaluation asynchrone avec une température hors de ses bornes
physiques, compte rendu après synchronisation, valeurs des autres points.
Avec CUDA et HIP, un noyau lent écrit les arguments sur le flux avant
l'évaluation : le flux est encore occupé quand la fonction rend la main
(asynchrone vérifié) et les valeurs sont justes (ordre vérifié).

MGIS : la copie de la structure d'options suit TFEL (vérifiée par le test
de disposition). MGIS n'utilise l'asynchrone ni pour les lois ni pour les
propriétés (étape 5) : les deux sont symétriques.

### 2026-10-08 (soir, fin), files SYCL et blocs de fils des propriétés matériau (R10, fin)

Deux écarts restaient entre les lois et les propriétés matériau. Ils sont
alignés sur les lois.

- **File SYCL non ordonnée.** Le noyau doit s'exécuter après les calculs
  déjà soumis à la file, qui peuvent produire ses arguments. Une file
  ordonnée le garantit. Les lois traitaient déjà la file non ordonnée :
  calcul asynchrone refusé, file vidée avant le noyau en synchrone. Les
  propriétés matériau ne le faisaient pas : leur noyau pouvait passer avant
  les calculs de ses arguments. Même règle maintenant (statut `-5` pour
  une évaluation asynchrone).
- **Fils par bloc de CUDA et HIP.** 64 fixes pour les propriétés matériau,
  128 arrondis au multiple de la taille du warp pour les lois. Même règle
  maintenant : 128 par défaut, taille du warp lue sur le GPU à l'exécution.

Tests : `test-sycl.cxx` vérifie le refus de l'asynchrone sur une file non
ordonnée (test 16). Le test de l'évaluation asynchrone utilise une file
ordonnée. Tests installés des propriétés matériau, après réinstallation de
TFEL : `cuda-nvcc`, `hip-hipcc`, `sycl-default-acpp`, `sycl-gpu-acpp`, tous
verts. MGIS : 45 tests sur 45 contre ce TFEL.

### 2026-10-08 (nuit), documentation de TFEL et de MGIS (R11)

La documentation datait en partie des versions officielles. Elle a été
relue en entier et mise à jour, en trois parties.

**Documentation du prototype** (`docs/`). Restructurée : vue d'ensemble,
détails de TFEL, détails de MGIS, validation, découpage pour les dépôts
officiels, plan pour mfem-mgis. Les anciens documents de travail sont
dans `docs/historique/`. Ce journal a reçu un index et la table des
phases.

**TFEL** (`tfel-dev`). Trois relectures indépendantes, puis corrections :

- pages web des deux interfaces `generic-parallel`, notes de version,
  pages de l'interface `generic` et de `tfel-config`, menu, liste des
  documentations ;
- commentaires de code : espace mémoire des lots, files SYCL, listes de
  backends, conventions des backends des propriétés matériau,
  documentation écrite dans le code généré (pas des tableaux, paramètres
  du constructeur, ordre des arguments) ;
- erreurs anciennes dans les fichiers touchés : commentaires déplacés,
  noms de paramètres faux, paramètres de patron documentés comme des
  paramètres ;
- section `@SharedData` ajoutée aux pages des mots-clés des DSL de lois.
  Ces pages ne sont pas régénérées en entier : TFEL y a corrigé des fautes
  directement, pas dans les fichiers d'aide dont elles sont tirées. La page
  du DSL `IsotropicStrainHardeningMisesCreep`, qui contenait un message
  d'erreur, est remplacée par la page générée.

Chaque affirmation de la documentation a été vérifiée dans le code. Neuf
constats des relectures étaient faux ou imprécis : ils ont été corrigés
avant d'être appliqués. Un défaut de CMake est noté, pas corrigé : les
tests `sycl-gpu-acpp` utilisent les fichiers de configuration de
`sycl-default-acpp`, installés seulement si cette configuration est aussi
choisie.

**MGIS** (`mgis`). Nouvelle page `generic-parallel.md`, reliée depuis
l'accueil, le menu et les pages qui en parlent. Notes de version,
installation et README mis à jour. Chaque fonction publique dit ce
qu'elle fait pour une bibliothèque `generic-parallel` : bascule, refus,
emplacement des tableaux. Erreurs anciennes corrigées : `update` décrit à
l'envers, noms des rotations des lois en grandes déformations, paramètres
mal annotés.

Seuls des commentaires ont changé dans le code : vérifié en comparant les
sources sans leurs commentaires (MGIS), en lisant le diff et en compilant
(TFEL et MGIS).

