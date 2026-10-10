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
| R12 | statistiques des propriétés matériau calquées sur les lois |
| R13 | audit final de TFEL et de MGIS, corrections, seconde relecture |
| D1 | MGIS : lots de lois, de modèles et de propriétés matériau |
| D2 | mfem-mgis : choix retenus, travail à venir ([docs/plan-mfem-mgis.md](docs/plan-mfem-mgis.md)) |

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
| [2026-10-08 (nuit, suite), calcul asynchrone dans MGIS (D1, étape 5)](#2026-10-08-nuit-suite-calcul-asynchrone-dans-mgis-d1-étape-5) |
| [2026-10-08 (nuit, fin), premier point en échec des propriétés matériau (R12)](#2026-10-08-nuit-fin-premier-point-en-échec-des-propriétés-matériau-r12) |
| [2026-10-08 (nuit, clôture), audit final de TFEL et de MGIS (R13)](#2026-10-08-nuit-clôture-audit-final-de-tfel-et-de-mgis-r13) |
| [2026-10-08 (fin), choix pour mfem-mgis (D2)](#2026-10-08-fin-choix-pour-mfem-mgis-d2) |

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

### 2026-10-08 (nuit, suite), calcul asynchrone dans MGIS (D1, étape 5)

TFEL offrait le calcul asynchrone pour tout. MGIS ne l'exposait pas.
C'est fait, le synchrone restant le défaut.

- Une politique asynchrone (`mgis::AsynchronousBatchExecutionPolicy`)
  reprend la file et la pile de la politique synchrone. Elle ajoute la
  ressource mémoire qui alloue le statut du calcul, et son emplacement.
- Chaque fonction de lot l'accepte : intégration, fonctions
  d'initialisation, post-traitements, rotations, évaluation des
  propriétés matériau. Elle rend un calcul en attente
  (`PendingBatchComputation`, `PendingBatchEvaluation`).
- Après la synchronisation de la file par l'appelant, `getResult` rend le
  résultat de la fonction synchrone. Le code est commun aux deux modes,
  sauf l'obtention du statut. MGIS ne synchronise rien : il ne dépend
  d'aucune bibliothèque GPU.
- Le calcul en attente garde la mémoire que le GPU écrit ou lit : le
  statut, et les matrices de rotation construites à partir des axes du
  matériau. Ces matrices sont allouées par la ressource de la politique :
  en asynchrone, elles marchent aussi pour une loi qui calcule sur le GPU.
- MGIS vérifie que le statut est accessible par l'hôte (initialisation,
  lecture) et par les calculs.
- Les propriétés matériau chargent la fonction de compte rendu et le
  nombre de statuts de leur bibliothèque.

Tests : sur CPU (`stlpar`), mêmes résultats qu'en synchrone, indice du
premier point en échec compris, refus d'une loi `generic` et d'un statut
mal placé. Sur GPU (CUDA, flux non par défaut, mémoire managée, sans HMM) :
intégration, rotation à partir des axes du matériau, bornes d'une
propriété matériau. MGIS : 45 tests sur 45.

TFEL : les tests `sycl-gpu-acpp` utilisent les fichiers de configuration
de `sycl-default-acpp`. Ils sont maintenant installés si l'une des deux
configurations est choisie (vérifié par une configuration avec
`sycl-gpu-acpp` seule).

Reste à discuter : l'indice du premier point en échec pour les propriétés
matériau (trois options dans [docs/tfel.md](docs/tfel.md#comparaison-entre-lois-et-propriétés-matériau)).

### 2026-10-08 (nuit, fin), premier point en échec des propriétés matériau (R12)

Les propriétés matériau `generic-parallel` reprenaient le statut de
l'interface `generic`, conçu pour un seul point. Elles ne disaient ni
combien de points avaient échoué, ni lequel en premier. C'était la
dernière différence avec les lois. Option retenue : un statut propre,
calqué sur celui des lois.

TFEL :

- statistiques (`mfront_gpmp_Statistics`) et statut de sortie
  (`mfront_gpmp_OutputStatus`) propres à l'interface : nombre de points
  en échec, premier d'entre eux, nombre de points hors bornes avec la
  politique `Warning`, première vérification en échec (pour nommer la
  variable dans le message) ;
- options d'exécution avec les membres de celles des lois : file,
  statistiques (asynchrone), pile. Le tableau de statuts et le symbole
  `<name>_generic_parallel_number_of_statuses` disparaissent ;
- chaque noyau calcule le statut de son point dans une lambda, et ne le
  fusionne qu'en cas d'événement : opérations atomiques pour CUDA, HIP,
  SYCL et Kokkos, réduction pour les algorithmes parallèles et OpenMP,
  comme pour les lois. Coût nul en cas de succès ;
- nouveau défaut de `nvc++` 25.11 (`-mp=gpu`, -O0 à -O3) : une lambda
  imbriquée qui capture par référence un pointeur capturé par le noyau
  n'écrit pas par ce pointeur. Contournement : capture par copie.
  Reproducteur : `results/nvhpc-bugs/omp-lambda-reference-capture.cxx`.

Le signe du rang de la variable hors bornes a été examiné. La doc de
l'interface `generic` le dit négatif pour une borne physique. Le code le
rend négatif dès que l'évaluation échoue, borne avec la politique
`Strict` comprise, et un test du master de MGIS l'attend. Le comportement
est gardé, l'interface `generic-parallel` le suit, et la doc est corrigée.
Deux défauts de l'interface `generic` sont corrigés : message de la borne
physique supérieure, `errno` non restauré.

MGIS : copies des nouvelles structures (dispositions vérifiées à la
compilation), `evaluate` rend un `MaterialPropertyEvaluationResult`
(statut, rang, indice du premier point en échec, message), le même pour
les deux interfaces, comme `BehaviourIntegrationResult`. L'indice est le
même en `generic` et en `generic-parallel` (testé).

Tests : propriétés matériau de TFEL sur les 16 configurations (tous les
backends et compilateurs), avec le nombre d'échecs, le premier point en
échec et les résultats peu fiables ; tests `generic` de TFEL (6412) ;
MGIS 45 sur 45, dont CUDA.

Il reste une différence, écartée : le tableau optionnel de statuts par
point des lois n'existe pas pour les propriétés matériau. MGIS ne l'expose
pour aucune des deux.

### 2026-10-08 (nuit, clôture), audit final de TFEL et de MGIS (R13)

Deux relectures indépendantes, en lecture seule, ont parcouru tout le code
`generic-parallel` de TFEL (lois, modèles, propriétés matériau, six
backends) et de MGIS. Chaque constat a été vérifié dans le code avant
correction. Les constats et leur traitement :

TFEL :

- Kokkos, lois en synchrone : réduction, alors que les propriétés
  matériau utilisaient des atomiques. Les lois passent aux atomiques sur
  un bloc de `Kokkos::SharedSpace`, comme les propriétés. Les deux
  interfaces exigent cet espace (assertion statique).
- SYCL, propriétés matériau : classe de copies générée à part, différente
  de celle des lois. Statistiques jamais copiées (écriture du GPU dans la
  mémoire de l'hôte au premier échec), allocation non testée, sortie
  indéterminée aux points en échec, `wait()` au lieu de
  `wait_and_throw()`, notion d'accessibilité différente. Les deux
  interfaces partagent maintenant `SYCLDeviceCopies.hxx`.
- Pile des fils : CUDA, HIP et Kokkos ne la vérifiaient pas pour les
  propriétés matériau, et `preserve_stack_size` y était sans effet. Les
  deux interfaces partagent `GPUStackSize.hxx` et `KokkosStackSize.hxx`.
- Le noyau devient un paramètre template de la fonction qui lit ses
  attributs. Avec un argument, la fonction était instanciée à
  l'identique dans plusieurs bibliothèques. L'éditeur de liens dynamique
  la résolvait vers une seule. Le runtime CUDA, lié statiquement, de
  cette bibliothèque ne connaissait pas le noyau des autres :
  `invalid resource handle`. Trouvé par le test CUDA des propriétés
  matériau, qui charge trois bibliothèques.
- Exceptions des runtimes (SYCL, Kokkos, algorithmes parallèles,
  OpenMP) : statut `-2` pour les propriétés matériau, `-6` pour les
  erreurs CUDA et HIP. Toutes en `-6` maintenant.
- Compte rendu asynchrone : fonction C `mfront_gpb_report` pour les lois,
  pendant de `<nom>_generic_parallel_report` des propriétés. Le compte
  rendu synchrone des lois l'utilise : un seul texte de message.
- Lot vide : accepté sans aucune vérification dans les deux interfaces
  (avant : `-2` avec SYCL sans file, `-7` pour une sortie uniforme).
  `errno` restauré à chaque retour (garde RAII).
- Lois SYCL avec `@SharedData` et sans file : `-3` au lieu de `-2`.
- Nombre de points en `size_t` dans les noyaux CUDA, HIP, Kokkos et avec
  AdaptiveCpp (`int` avant, limite à 2^31 points). AdaptiveCpp déporte
  bien la boucle avec un index `size_t` (vérifié).
- Doc : coût des réductions (quelques entiers copiés par appel),
  statistiques de `stlpar` et `openmp` accessibles par l'hôte.
- Tests ajoutés : `preserve_stack_size` (propriétés CUDA, HIP, Kokkos ;
  lois CUDA et HIP), statistiques en mémoire hôte avec SYCL, lois SYCL
  avec file nulle et file non ordonnée, lois SYCL avec la politique de
  transfert `none` (configuration `sycl-default-acpp-2`), exécution d'une
  propriété compilée sans vérifications, ensemble vide.

MGIS :

- Affectation par déplacement des objets en attente supprimée : avec des
  ressources différentes, elle recopiait les statistiques pendant que le
  GPU les écrivait (écriture après libération). Le déplacement reste.
- Propriétés matériau : une évaluation `generic-parallel` ne rend que
  `1`, `0`, `-1`, `-4`. Tout autre statut lève une exception, dès la
  soumission en asynchrone (le `-2` était perdu).
- Conversions des grandes transformations sur un gestionnaire : contrôle
  de l'accès hôte.
- Rotations par lots : emplacement des tableaux en argument (`HOST` par
  défaut), comme `evaluate`. `buildRotationMatrices` construit une fois
  les matrices dans une ressource donnée. Rotation en place documentée.
- Emplacement des valeurs données par `std::span` au gestionnaire et des
  tableaux des fonctions de lot : le même que les tableaux du
  gestionnaire (contrat documenté). `allocatePostProcessingVariables`
  avec une ressource.
- Tests ajoutés : échecs partiels sur GPU (premier échec par atomiques,
  synchrone, asynchrone, plage), résultat non fini (`-4`, hôte et CUDA),
  emplacements refusés, matrices de rotation construites une fois,
  rotation sans politique sur GPU, déplacement des objets en attente.

Limites gardées, avec leur raison (voir [docs/mgis.md](docs/mgis.md#limites-actuelles)) :
agencement point par point du `MaterialDataManager` (une API sur des vues
est à définir avec mfem-mgis), pas de variante avec `Context`,
`MaterialFunctionManager` sur l'hôte, `update`/`revert` sur l'hôte, pas
de tableau de statuts par point, tests `generic-parallel` absents de la
CI de MGIS tant que TFEL master n'a pas l'interface.

Seconde relecture, sur le diff des corrections, avec le rapport du premier
audit en entrée (statut de chaque constat, défauts introduits, asymétries
restantes) :

- TFEL : tous les constats corrigés, aucun défaut introduit de classe A ou
  B. Points mineurs traités : commentaires sur le coût des réductions,
  statut de succès d'un ensemble vide (`1` ou `0` selon la convention),
  remontée des erreurs SYCL par le gestionnaire asynchrone de la file,
  fuite possible si la description d'une copie SYCL ne peut être rangée,
  fonctions C des propriétés CUDA et HIP sans `try`.
- MGIS : deux trous trouvés et corrigés. Les rotations sans politique
  passaient l'emplacement `HOST` aux fonctions de lot : elles refusaient
  toute bibliothèque `device`, régression de l'étape précédente. Elles
  prennent maintenant l'emplacement, comme les autres.
  `allocatePostProcessingVariables` sans ressource donnait de la mémoire
  hôte à une loi qui ne peut pas la lire : refusé maintenant. Une
  ressource nulle a le même sens partout (ressource par défaut).

Tests :

- propriétés matériau : les 16 configurations ;
- lois : 3 configurations en entier (`openmp-host-gcc`, `stlpar-parunseq-gcc`,
  `sycl-default-acpp-2`), les 13 autres en passe ciblée (intégration,
  échecs, fonction d'initialisation, post-traitement, rotations, données
  partagées pour SYCL, files SYCL, pile préservée sur CUDA et HIP) ;
- MGIS : 45 tests sur 45 contre le TFEL installé, CUDA compris.

`scratch/dev-test.sh` ne copiait pas le fichier de configuration quand
`test.check` le passe à `build-libraries.sh` (`sycl-default-acpp-omp`) :
corrigé.

En reconstruisant les bibliothèques de test de MGIS avec `-Wshadow`, un
avertissement est apparu dans le code des post-traitements de lots (une
variable locale masquait un membre) : corrigé.

### 2026-10-08 (fin), choix pour mfem-mgis (D2)

Le plan pour mfem-mgis est tranché. Il est décrit dans
[docs/plan-mfem-mgis.md](docs/plan-mfem-mgis.md). En bref :

- deux buts : explorer (mesures, variantes) et proposer un code propre aux
  développeurs de mfem-mgis. Tout est choisissable par des options, le
  chemin actuel reste le chemin par défaut ;
- deux points d'entrée dans MFEM : assemblage partiel de `NonlinearForm`
  et opérateur propre à mfem-mgis, sur les mêmes noyaux par lot ;
- deux rangements des données aux points de Gauss : point par point et
  composante par composante. MGIS n'accepte aujourd'hui que le premier :
  une API sur des vues est la première étape ;
- trois options de l'opérateur tangent (sans matrice, point fixe élastique,
  tangente assemblée) et un solveur non linéaire de point fixe ;
- Ginkgo dans la fabrique de solveurs, solveurs directs gardés sur CPU ;
- MFEM 4.9, la version utilisée jusqu'ici par mfem-mgis, à compiler avec
  CUDA. La spécification de `deps/spack/mmm-gpu/spack.yaml` (4.10) sera
  changée avant la compilation.

Les sources de mfem-mgis (master `a011d97d6f`), de mfem-mgis-examples
(master `a242567`) et de mm-opera-hpc (main `c31c4bc`) sont suivies dans
ce dépôt depuis leur base, sans modification (`scripts/sync-repository.sh`).

### 2026-10-09, API de vues de MGIS (D2, étape 1)

Première étape du plan de mfem-mgis : MGIS accepte des tableaux rangés
composante par composante, ou appartenant à l'appelant, décrits par des
vues à pas quelconques. La proposition (version 4, après trois relectures)
est versionnée dans
[docs/historique/conception-vues-mgis-2026-10-08.md](docs/historique/conception-vues-mgis-2026-10-08.md).
L'état est décrit dans [docs/mgis.md](docs/mgis.md#vues-et-rangements).

Ce qui est fait, dans `mgis/` (patchs `0024` à `0033`) :

- `FieldView` (deux pas, type régulier, contrat vérifié à la
  construction), `MaterialStateView` et `MaterialDataView`, `make_view` ;
- toutes les fonctions du gestionnaire sur une vue, sans politique (espace
  de travail explicite), avec politique synchrone ou asynchrone. Les
  fonctions du gestionnaire appellent celles de la vue : un seul chemin ;
- rangement d'un gestionnaire (`layout` de l'initialiseur), suivi par tout
  ce qui est donné par point. `save` et `restore` restent point par point ;
- `MaterialFunctionManager` avec initialiseur et vues *strided coalesced* ;
  rangement exposé en C, Fortran et Python ;
- `FieldOperations` : `update`, `revert` et `updateValues` sur les
  tableaux de l'appelant, sur GPU compris. En-têtes CUDA, HIP, SYCL, Kokkos
  et OpenMP, compilés par l'appelant ;
- rotations et `evaluate` sur vues ; les surcharges sur tableaux passent
  par elles ;
- plus d'allocation sur le tas par appel dans le cas courant (tampons sur
  la pile), pool recommandé pour le statut des calculs asynchrones.

Défauts trouvés en chemin et corrigés : restauration HDF5 d'une variable
externe à plusieurs composantes (écriture hors du tableau), ressources
mémoire CUDA, HIP et Kokkos inutilisables sous un pool (alignement refusé),
initialiseur sans constructeur en Python, message d'erreur des rotations
de la tangente. Trouvé et non corrigé : `MaterialAxisStorage` non exporté
(copier un `RotationMatrix3D` hors de MGIS ne se lie pas).

Vérifications : 48 tests de la construction principale (CUDA compris),
94 de la construction avec HDF5 et les liaisons C, Fortran et Python, et
les en-têtes GPU compilés et exécutés avec les cinq chaînes
(`scripts/check-mgis-gpu-headers.sh`). Les résultats sur vues sont
identiques bit à bit à ceux du gestionnaire, dans tous les rangements. Des
mutations volontaires du code (recopie omise, pas ignoré, matrices mal
lues, mise à jour incomplète) sont détectées par les tests.

Documentation de MGIS : page `docs/web/views.md`, page des lots et notes de
version mises à jour.

### 2026-10-09 (suite), relectures de l'API de vues et corrections

Avant de déclarer l'étape finie, un audit de couverture et trois
relectures indépendantes (justesse, conception, puis relecture du diff des
corrections). Patchs `0035` à `0053` (le commit `a2f8b21` annonçait à tort
`0035` à `0057`).

- **Couverture** : l'audit a relevé les cas non testés (prédictions,
  vitesse du son, modèle, énergies, pool de threads, plages, gestionnaire
  composante par composante hors cas simple, vues sur GPU, refus attendus).
  Tous ont un test, dont une nouvelle loi de test qui calcule les énergies.
  Chaque calcul sur une vue part maintenant d'une copie de l'état avant
  calcul : un calcul qui ne fait rien est détecté.
- **Coût** : les contrôles construisaient à chaque appel les messages
  d'erreur, jamais affichés (une trentaine d'allocations par appel). Un
  test compte les appels à `operator new` : plus aucune allocation dans le
  thread appelant pour l'intégration, l'initialisation, le
  post-traitement, la fin de pas, les rotations et `evaluate`.
- **Défauts corrigés** : pool de threads (une tâche en échec laissait les
  autres lire des objets détruits), sorties d'un point en échec recopiées
  depuis le point précédent, message d'erreur lu dans une exception
  détruite, conversions des grandes transformations sans contrôle de la
  sortie, copies de vues vides, masse volumique (règle écrite), noyaux GPU
  sans division par valeur, file SYCL gardée par référence.
- **API complétée** : `getPoints`, blocs de la tangente d'une vue,
  `updateValues` avec contexte, extraction par opérations sur le
  gestionnaire, conversions sur vues du gestionnaire, en-tête propre de
  l'espace de travail, refus d'une loi temporaire.
- **Défauts anciens de MGIS corrigés** : copies des matrices de rotation,
  forme des opérateurs tangents en Fortran, tableaux numpy à pas lus comme
  contigus, contexte alloué par `update` et `revert`.
- **Laissé de côté** : les vues dans les liaisons C, Fortran et Python (le
  prototype ne travaille qu'en C++), le contrôle du recouvrement entre deux
  sorties différentes (documenté).

Leçon : les tests des liaisons sont construits par la cible `check`, pas
par `build-tests`. Une mutation du code Fortran passait inaperçue avec un
exécutable périmé ; elle est détectée une fois la cible `check` utilisée.

Vérifications finales : 50 tests (construction principale, CUDA compris),
96 tests (HDF5 et liaisons C, Fortran, Python, par `check`), les cinq
chaînes GPU (`scripts/check-mgis-gpu-headers.sh`), et des mutations
volontaires détectées par les nouveaux tests.

### 2026-10-09 (nuit), verrouillage de MGIS, liste des solveurs, environnement spack

- **MGIS verrouillé.** Nouveau script `build-mgis.sh` : construction
  neuve de `build/mgis` (50 tests, CUDA compris) et installation dans
  `install/mgis`, où mfem-mgis le trouvera ; `build-mgis.sh full` pour
  `build/mgis-full` (HDF5, liaisons C, Fortran, Python : 96 tests). Les
  deux reconstruites depuis des dossiers vides, vertes.
- **Fuite d'environnement corrigée.** La première construction neuve a
  échoué : l'environnement spack `tfelgpu_local`, activé dans le shell,
  fournissait HDF5, et sa vue contient un ancien TFEL dont les en-têtes
  passaient avant ceux du TFEL du prototype. `scripts/env.sh` et
  `build-tfel.sh` neutralisent maintenant cet environnement pour
  `mmm-gpu`. Une configuration neuve de TFEL donne le même cache que
  `build/tfel-mmm-gpu`, aux bibliothèques CUDA inutilisées près (cuFile,
  nvfatbin) et au hachage git affiché par TFEL, figé à la configuration
  (`9f10d65a4` au lieu de `631019468`, cosmétique).
- **MGIS, patch 0054** : les lois CUDA des tests sont compilées avec les
  options de la configuration `cuda-nvcc` de TFEL (avertissements de nvcc
  sur les en-têtes de TFEL masqués, appels de fonctions hôte depuis le GPU
  en erreur). Plus aucun avertissement de nvcc.
- **Liste des solveurs** : [docs/solveurs-gpu.md](docs/solveurs-gpu.md).
  Points notables, vérifiés dans les sources : la stratégie
  « Elasticity » de BoomerAMG dans mfem-mgis arrête le programme avec
  hypre GPU (`SetElasticityOptions`) ; MFEM 4.10 enveloppe cuDSS et les
  solveurs distribués de Ginkgo ; AmgX 2.5.0 précède les corrections pour
  CUDA 13 ; l'enveloppe AmgX de MFEM est en double seulement.
- **Environnement spack** : MFEM 4.9 avec CUDA, Ginkgo, PETSc (Kokkos),
  SUNDIALS, AmgX, SuperLU_DIST, STRUMPACK, MUMPS, UMFPACK ; cuDSS seul.
  Problèmes rencontrés et réglés : cache de spack partagé entre deux
  installations (index des patchs) ; BLAS du HPC SDK choisi par défaut
  (OpenBLAS imposé) ; Kokkos 5 sans variante `cuda_lambda` exigée par la
  recette de PETSc (recette modifiée dans `mmm_gpu`) ; AmgX 2.5.0 trop
  ancien pour CUDA 13 (commit de la branche principale) ; OpenMPI avec
  CUDA impossible sans reconstruire LLVM et TBB (OpenMPI sans CUDA,
  PETSc lancé avec `-use_gpu_aware_mpi 0`) ; préfixe CUDA dont `bin` était
  un lien, si bien que nvcc prenait les en-têtes du SDK sans cuBLAS
  (`make-nvhpc-cuda-prefix.sh` corrigé, préfixes 13.0 et 12.9 refaits).

### 2026-10-09 (nuit, suite), plan détaillé de mfem-mgis et relecture

- **Cartographie de mfem-mgis** (lecture seule) : boucle de résolution,
  intégrateurs, données aux points, solveurs, paramètres, tests, avec les
  références `fichier:ligne`. Elle sert de base au plan.
- **Plan détaillé** : [docs/plan-mfem-mgis-detaille.md](docs/plan-mfem-mgis-detaille.md),
  onze lots (L0 à L10) avec fichiers touchés, tests et critères de fin,
  paramètres proposés, stratégie de test, campagne de mesures.
- **Relecture indépendante** de la liste des solveurs et du plan : 50
  constats, tous intégrés. Les plus importants :
  - le `Mult` de mfem-mgis n'intègre pas la loi (c'est une action du
    Newton de mfem-mgis) : les solveurs non linéaires de MFEM, SUNDIALS
    et PETSc demandent un opérateur qui l'intègre ; en mode PETSc, un
    échec de la loi arrête le programme ;
  - le Newton de mfem-mgis n'est pas celui de MFEM (Newton inexact et
    recherche linéaire à y écrire) ; L-BFGS de MFEM ignore le
    préconditionneur ;
  - `mfem::forall` ne lance de noyau CUDA que dans un fichier compilé par
    nvcc ;
  - règle d'intégration à 27 points sur l'hexaèdre d'ordre 1 (option de
    règle ajoutée : environ 3,4 fois moins de calcul avec 8 points) ;
  - les forces thermodynamiques ne sont pas tournées à l'intégration ;
  - le chemin actuel suppose partout le rangement point par point ;
  - validité des copies hôte et GPU de `mfem::Memory` autour des appels à
    MGIS ;
  - `SystemHypreBoomerAMG` écrase deux réglages de BoomerAMG à chaque
    `SetOperator` ; `AmgXSolver::UpdateOperator` est une réutilisation
    complète, pas partielle.

### 2026-10-09 (nuit, fin), compilation de l'environnement

- Deuxième passe : tout s'installe (hypre, Ginkgo, AmgX, SuperLU_DIST,
  STRUMPACK, SUNDIALS, MUMPS, cuDSS, kokkos-kernels) sauf PETSc : ses
  fichiers Kokkos étaient compilés sans OpenMP, que les en-têtes d'un
  Kokkos compilé avec OpenMP exigent. Corrigé par `petsc +openmp`.
- Troisième passe : PETSc s'installe ; MFEM échoue dès la configuration,
  compilé en C++17 alors qu'Umpire 2026 (par camp) exige C++20. Corrigé
  par `mfem cxxstd=20`.
- Kokkos a été reconstruit (`+pic`, exigé par kokkos-kernels partagé) : les
  tests Kokkos de TFEL passent (lois et propriétés matériau, `kokkos` et
  `kokkos-host`, 4/4), et les en-têtes GPU de MGIS avec Kokkos aussi.
- Leçon d'outillage : dans les surveillances de journaux, `grep` du shell
  est une fonction (ugrep) qui ne remontait rien ; attendre plutôt la fin
  d'un processus par son PID.

### 2026-10-09 (matin), MFEM 4.10, OpenMPI compatible GPU, constructions surveillées

- **Tout le prototype passe en MFEM 4.10** (cuDSS, Ginkgo distribué,
  `GroupCommunicator` sur le périphérique), avec les backends OpenMP et
  libCEED. Les faits vérifiés sur 4.9 ont été revérifiés sur 4.10.
- **OpenMPI avec CUDA.** Son hwloc interne activait le greffon OpenCL avec
  un chemin invalide ; la recette `openmpi` de `mmm_gpu` le construit avec
  `--disable-opencl`. Le hwloc de spack avec CUDA aurait reconstruit LLVM
  et TBB.
- **Constructions robustes.** Le régulateur de mémoire suspend tous les
  processus de la construction sauf le plus gros (qui tient la mémoire et
  doit finir) ; chaque phase de `install-mmm-gpu.sh` est relancée jusqu'à
  trois fois ; une recette `mfem` garde la date de `config.mk` pour qu'une
  construction interrompue reprenne sans tout recompiler. Un superviseur
  hors de la session relance l'installation jusqu'au succès.
- **Deux fichiers de MFEM 4.10 démesurés avec CUDA 13** (sm_120) :
  `fem/lor/lor_batched.cpp` (`ptxas` à 23 Go) et
  `fem/integ/bilininteg_convection_ea.cpp` (`cicc` pendant plus de deux
  heures : noyau 3D spécialisé jusqu'aux ordres 8×9, boucles déroulées).
  Avec deux tâches, le noyau a tué `ptxas` faute de mémoire malgré la
  pause du régulateur ; la phase repart avec une tâche, en gardant les
  objets compilés.
- **Plan** : confronté à tout le code de mfem-mgis (82 fonctionnalités) et
  complété ; techniques de point fixe de l'étude préalable ajoutées (point
  fixe sur `K0`, accélérateurs de MTest portés, Newton modifié, AMG de `K0`
  gelé, simple précision, prédiction).

### 2026-10-09 (midi), MFEM 4.10 installé, lot L0 commencé

- **MFEM 4.10 installé** (12 h 32), avec CUDA, OpenMP, libCEED, Ginkgo,
  AmgX, cuDSS, PETSc, SUNDIALS, SuperLU_DIST, STRUMPACK, MUMPS, SuiteSparse,
  Umpire et gslib. L'exemple `ex1p` tourne sur CPU, sur GPU en assemblage
  partiel, avec BoomerAMG sur GPU à deux processus (OpenMPI compatible
  CUDA) et avec le backend `ceed-cuda`.
- **Dernières difficultés de construction** :
  - `fem/integ/bilininteg_convection_ea.cpp` : cicc grossissait de 8 Go
    toutes les 20 minutes (17 Go après 2 h 30) ; patch qui garde les
    spécialisations 3D jusqu'à l'ordre 3 (le noyau générique prend la
    suite), 1 min 42 ;
  - édition de liens : `-lnvToolsExt` (absent de CUDA 13), demandé par les
    options d'AmgX de MFEM, retiré ;
  - `libgs.a` de gslib non relogeable : recompilée avec `-fPIC`.
  L'installation actuelle a été finie à la main (même contenu que les
  recettes) ; une reconstruction entièrement par les recettes reste à
  faire (`deps/spack/README.md`).
- **Lot L0** : les trois dépôts se construisent et passent leurs tests
  contre MFEM 4.10 avec CUDA (périphérique CPU) : mfem-mgis 128/128,
  mfem-mgis-examples 20/20, mm-opera-hpc 6/6. Corrections :
  - mfem-mgis : `FindMFEM.cmake` tronquait les valeurs de `config.mk`
    contenant un signe égal ; options de nvcc traduites pour le compilateur
    C++ ; toolkit CUDA ajouté ; `tfel.cmake` lit plusieurs dossiers rendus
    par `tfel-config` (TDLS) ; tolérance absolue des Krylov de hypre
    (code mort et mauvais appel pour FGMRES) avec un test qui détecte la
    régression ; option `--device` avec un test ;
  - bases des exemples et de mm-opera-hpc : dépôts de référence
    (`latug0/mfem-mgis-examples` master `8397eec`,
    `rprat-pro/mm-opera-hpc` main `8bd56eb`) ; les anciennes bases, prises
    sur des miroirs plus anciens, ne compilaient plus avec mfem-mgis ;
  - mm-opera-hpc : la boucle des bulles, sur des données de l'hôte,
    utilisait `mfem::forall`, qui exige depuis MFEM 4.10 une compilation
    CUDA quand MFEM a CUDA (et lancerait la boucle sur le GPU) :
    `mfem::hypre_forall_cpu`.

### 2026-10-09 (après-midi), relecture du lot L0 et corrections

Une relecture indépendante des commits du lot L0 a donné 25 constats.
Tous sont traités :

- **Option `--device`** : le nom court `-dev` n'était pas reconnu par
  `initialize` (le calcul tournait sur l'hôte sans message), comme les
  noms courts des options de PETSc. Le périphérique est maintenant créé
  avant l'initialisation de PETSc. Chaque processus d'un nœud prend un GPU
  du nœud à tour de rôle. Le nombre de GPU n'est demandé que pour une
  configuration `cuda` ou `hip` : MFEM s'arrête sinon quand il n'a pas de
  GPU, ou sur une machine sans GPU.
- **Sources contenant des noyaux de MFEM** : module CMake installé
  `MFEMMGISDeviceCode.cmake`, chargé aussi quand un projet trouve MFEM
  avant mfem-mgis. nvcc et son compilateur hôte sont par défaut ceux de
  MFEM (`mpic++` d'OpenMPI, qui trouve `mpi.h`). Architectures lues sous
  toutes les formes des options de nvcc. Erreurs explicites pour HIP et
  quand le langage CUDA n'est pas activé. Runtime CUDA partagé
  (`CUDA::cudart`, comme le reste du code de mfem-mgis ; MFEM embarque sa
  propre copie statique du runtime). Définitions et options de visibilité du C++ données aux sources
  CUDA. Tests du GPU hors du bloc TFEL, option `enable-gpu-tests`, test
  du GPU à deux processus.
- **TFEL** : seules les options `-I` et `-L` de `tfel-config` sont lues ;
  `env.sh` ajoute chaque dossier de bibliothèques (mfem-mgis et
  mm-opera-hpc).
- **Documentation** : `AbsoluteTolerance` des Krylov de hypre ; entiers
  32 bits de hypre et absence de `+gpu-aware-mpi` (`deps/spack/README.md`).
- **Scripts du prototype** : LLVM et les autres racines indépendantes de
  MFEM sont construits avec 4 tâches (et non 1 ou 2) ; le régulateur de
  mémoire reprend les processus qu'il a suspendus, même s'ils ont quitté
  l'arbre et même quand il s'arrête ; les dossiers de construction sont
  reconfigurés à chaque appel ; patch de la convection limité à CUDA 13 ;
  recette PETSc décrite comme copie complète.
- **Boucle des bulles de mm-opera-hpc** : `hypre_forall_cpu` ajoute une
  boucle OpenMP quand hypre est compilé avec OpenMP, ce que l'ancien
  `mfem::forall` sans périphérique ne faisait pas. Sans effet ici (hypre
  sans OpenMP), à surveiller si hypre passe en OpenMP (sursouscription
  avec plusieurs processus par nœud).

Reconstruction complète dans un dossier neuf : mfem-mgis 132/132,
mfem-mgis-examples 20/20, mm-opera-hpc 6/6, dont le noyau de MFEM sur le
GPU à un et deux processus.

### 2026-10-09 (après-midi), profilage et GPU

- MGIS : `mgis::setProfilingSynchronization` fixe une fonction appelée au
  début et à la fin de chaque section de profilage active, avant la
  lecture du temps (aucune par défaut). Sans elle, la durée des noyaux
  lancés de façon asynchrone sur le GPU tombait dans la première section
  suivante qui les attend. Test exact du nombre d'appels.
- mfem-mgis : avec un périphérique CUDA ou HIP, `initialize` y met une
  synchronisation du GPU.
- Exemples ssna303 en 3D (ex4) : option `--profiling`. Sur un pas, à deux
  processus, l'intégration des lois (dans
  `NS::processNewUnknownsEstimate`) prend déjà 25 % du temps, la tangente
  et la résolution 24 % chacune.
- Suites : MGIS 50/50, mfem-mgis 132/132, exemples 20/20, mm-opera-hpc
  6/6.

### 2026-10-09 (après-midi), seconde relecture des corrections du lot L0

Une seconde relecture indépendante des corrections a confirmé 18
constats corrigés et en a trouvé 6 partiels et 1 non corrigé, plus 9
nouveaux défauts mineurs (gravité 1 ou 2). Tous sont traités :

- choix du GPU de chaque processus aussi pour l'alias `gpu` de MFEM (et
  `raja-gpu`, `ceed-gpu`...), qui mettait tous les processus sur le
  premier GPU ; erreur si le périphérique de MFEM est déjà configuré ;
  le test vérifie le GPU de chaque processus ;
- toolkit CUDA pris à côté du nvcc de MFEM ; CMake 3.18 vérifié dès que
  MFEM utilise CUDA ; avertissement si le compilateur CUDA n'est pas celui
  de MFEM ; langage CUDA vérifié dans le dossier de la cible ; C++20 au
  moins pour les sources CUDA ; module documenté dans le guide
  d'installation ; avec un MFEM HIP, la configuration n'échoue plus ;
- `env.sh` : une ligne par dossier de bibliothèques de TFEL, écrite à la
  configuration (compatible zsh) ;
- régulateur de mémoire : un processus suspendu qui reçoit un signal (par
  exemple SIGTERM de spack) est repris aussitôt pour le recevoir (essayé) ;
- `build-mfem-mgis.sh` retire du cache les compilateurs CUDA d'un ancien
  dossier.

Nouveau dossier `cas/` : cas de la campagne propres au prototype. `poutre`
(plasticité J2, Newton de mfem-mgis) ; ses tests comparent la norme de la
solution finale, identique pour PCG et BoomerAMG, MUMPS et UMFPACK, et un
test vérifie que la comparaison détecte une solution fausse.
`cas/mesures-reference.sh` mesure les temps de référence du lot L0.

Suites : mfem-mgis 132/132, exemples 20/20, mm-opera-hpc 6/6, cas 4/4.

### 2026-10-09 (après-midi), temps de référence, fin du lot L0

Temps de référence du chemin actuel sur CPU pour les cas de la campagne :
`results/reference-L0-2026-10-09/README.md` (poutre, ssna303 en 3D,
polycristal, thermomécanique). La résolution linéaire (jusqu'à 61 %) et
l'intégration des lois (jusqu'à 55 %) dominent ; avec MUMPS, la prédiction
refactorise l'opérateur élastique à chaque pas (40 % du temps de la poutre
raffinée).

Comparaison avec l'ancienne pile (MFEM 4.9), même cas et même Newton :
PCG et BoomerAMG un peu plus rapides (15,4 s contre 17,5 s), MUMPS 1,55
fois plus lent (65,6 s contre 42,2 s), dont environ 10 % dus à OpenMPI
compilé avec CUDA (`--mca accelerator null`). Le reste est à éclaircir
(MUMPS 5.9.1, OpenBLAS avec OpenMP, MFEM 4.10) ; MUMPS avec METIS est à
essayer à la reconstruction de MFEM.

Le lot L0 est terminé : trois dépôts construits et testés contre MFEM
4.10 avec CUDA, option `--device`, compilation des noyaux de MFEM,
sections de profilage synchronisées avec le GPU, temps de référence.

### 2026-10-09 (après-midi), lot L1 : données aux points dans la mémoire de MFEM

Étapes faites (mfem-mgis, patchs 0016 à 0020) :

- **Tableaux des matériaux dans des `mfem::Vector`** (`MaterialArrays`,
  base privée de `Material` déclarée avant le gestionnaire de MGIS) :
  gradients, forces, variables internes et énergies des deux états, blocs
  de la tangente. MGIS les utilise sans copie (emplacement hôte). MGIS ne
  met plus F à l'identité quand on lui donne les gradients : mfem-mgis le
  fait, dans les deux rangements.
- **Validité des copies** : MGIS écrit par des pointeurs bruts, que MFEM
  ne voit pas. Chaque accès hôte de mfem-mgis passe d'abord par
  `synchronizeOnHost` (intégration, fin de pas, retour arrière,
  post-traitements, accesseurs ; lecture seule par
  `synchronizeOnHostForReading`). Sans GPU, rien ne change.
- **Éléments d'un matériau lisibles sur le GPU** : tableaux des éléments
  et des décalages de leurs points dans `PartialQuadratureSpace`.
- **Copies des états par noyaux de MFEM** (`batch::MFEMFieldOperations`,
  première source de la bibliothèque compilée comme du CUDA). La fin de
  pas et le retour arrière passent par `Material`, y compris pour les
  modèles ponctuels, qui ne synchronisaient pas.
- **Rangement composante par composante** : choisi à la construction du
  matériau ; les accesseurs fondés sur un pas entre points le refusent.

Reportés : propriétés matériau et variables externes en mémoire de
l'utilisateur (sans objet tant que MGIS reçoit les pointeurs de l'hôte),
paramètre `QuadratureDataLayout` du problème (pas encore de chemin par
lots), variante des tableaux alloués par MGIS.

Tests : `MaterialArraysTest` et `MFEMFieldOperationsTest` sur les
périphériques `cpu`, `debug` et `cuda`. Suites : mfem-mgis 138/138,
exemples 20/20, mm-opera-hpc 6/6, cas 4/4. (Correction après relecture :
cette entrée disait que chaque synchronisation était vérifiée par
mutation ; c'était faux pour l'intégration, les post-traitements et le
retour arrière. Voir l'entrée suivante.)

### 2026-10-09 (après-midi), relecture du lot L1 et corrections

Une relecture indépendante du lot L1 a donné 27 constats, aucun de
gravité 3 ou plus sur le chemin actuel (rien n'écrit encore les tableaux
sur le GPU), mais quatre qui deviendraient des résultats faux dès qu'un
noyau les écrirait. Tous sont traités :

- **Points d'accès hôte manquants** : intégration par PETSc (dans
  l'assemblage du résidu, sans passer par `integrate`) : synchronisation
  au début de `solve` ; modèles ponctuels (`computeNextState`) ; fonctions
  publiques de calcul des forces moyennes et des résultantes ;
  évaluateurs de contraintes et de gradients tournés ; boucles de
  mm-opera-hpc (bulles) et de l'exemple ex8 qui lisent les tableaux
  directement.
- **Post-traitements en lecture seule** : ils gardaient valides les
  copies du GPU, au lieu de les invalider à chaque pas.
- **Règle de validité complète** : une phase d'accès hôte commence par la
  synchronisation, et aucun accès au GPU n'a lieu pendant cette phase ;
  les vues rendues par les accesseurs ne valent que jusque-là.
- **Tableaux privés** : plus de réallocation possible derrière MGIS ;
  accès par pointeurs (`Read`, `Write`, `ReadWrite`).
- **Fin de pas sur l'hôte** par `mgis::HostFieldOperations`, qui garde les
  copies denses (les noyaux génériques étaient 3,5 à 12 fois plus lents) ;
  les opérations par noyaux de MFEM restent pour des tableaux sur le GPU.
  Avec un MFEM compilé avec HIP, leur source n'est plus compilée (la
  configuration échouait).
- **CMake** : runtime CUDA porté par la bibliothèque elle-même, contrôle
  de l'API de MGIS à la configuration.
- **Types** : tailles en `std::size_t` (les contrôles « trop grand pour
  MFEM » ne pouvaient pas se déclencher), exceptions des synchronisations
  rapportées dans le contexte au lieu d'arrêter le programme.
- **Rangement composante par composante** refusé aussi par les
  évaluateurs tournés.
- **Tests** : chaque point d'entrée (accesseurs, intégration,
  post-traitements, retour arrière, fin de pas) est appelé après une
  écriture de tous les tableaux sur le GPU ; fin de pas et retour arrière
  dans les deux rangements ; six accesseurs refusés ; loi qui calcule une
  énergie ; matériau sans point ; copies et remplissages dans les deux
  sens. Vérifié par mutation : les synchronisations de l'intégration, du
  retour arrière, de la fin de pas, des blocs de la tangente et des
  accesseurs font échouer le test quand on les retire (périphérique
  `debug`, et `cuda` quand des valeurs sont comparées). Exception : la
  synchronisation globale avant les post-traitements double celle que font
  déjà les fonctions qu'ils appellent ; la retirer ne fait rien échouer.

### 2026-10-09 (fin d'après-midi), seconde relecture du lot L1

Une seconde relecture des corrections a confirmé 22 constats corrigés, 5
partiels, et trouvé 10 défauts nouveaux mineurs (gravité 1 ou 2) et un
oubli plus ancien. Traités :

- synchronisation ciblée sur un tableau : les accesseurs, les évaluateurs
  et le calcul des forces moyennes ne recopient plus que le tableau lu
  (avant : les onze tableaux, environ 15 fois trop pour lire les
  contraintes) ;
- opérateurs linéarisés publics (`getLinearizedOperators`) : synchronisés
  à leur création ;
- exceptions : `integrate` et les évaluateurs les rapportent au lieu
  d'arrêter le programme ; le modèle ponctuel synchronise son statut entre
  processus même après une exception ;
- les évaluateurs tournés utilisent le contrôle de rangement commun ;
- CMake : l'API de MGIS est vérifiée par la compilation d'un test (et non
  par la présence d'un en-tête) ; la configuration installée cherche le
  toolkit CUDA utilisé à la construction ; avec HIP, l'en-tête des
  opérations par noyaux n'est plus installé ;
- noms des méthodes selon les règles de mfem-mgis (`readWrite`, `write`) ;
- tests : modèle ponctuel, rotation des forces d'une loi orthotrope (et
  son refus en composante par composante), forces moyennes appelées
  directement, quatre accesseurs de plus, énergies, blocs de la tangente
  après la fin de pas, matériau sans point, fichiers de sortie distincts
  par périphérique.

Mutations (chaque synchronisation retirée tour à tour) : toutes font
échouer le test sous `debug`, et presque toutes sous `cuda`. Deux
exceptions connues : la synchronisation globale avant les
post-traitements (redondante) et celle du début de `solve`, qui ne sert
qu'au chemin PETSc (non testé sur le périphérique `debug`).

Limite connue, comme le reste de mfem-mgis : une erreur locale à un
processus (ici une exception de synchronisation, possible seulement avec
MFEM compilé avec les exceptions) est rapportée par ce processus seul ;
les autres peuvent attendre dans une opération collective.

### 2026-10-09 (soir), lot L2 : lois intégrées par lots sur CPU

Paramètre `BehaviourIntegration` des problèmes non linéaires :
`PointByPoint` (défaut, chemin actuel) ou `Batch`. Valeur par défaut
modifiable par la variable d'environnement
`MFEM_MGIS_BEHAVIOUR_INTEGRATION`, pour passer une suite de tests ou une
application dans ce mode sans toucher au code.

Par lots :

1. les gradients sont calculés élément par élément (code de
   `implementIntegrate` sans l'appel de la loi), sauf pour une
   prédiction ;
2. chaque intégrateur intègre sa loi en un seul appel de MGIS
   (`mgis::behaviour::integrate` sur le gestionnaire de données du
   matériau) : point par point dans MGIS pour une bibliothèque `generic`,
   en un lot pour une bibliothèque `generic-parallel` (backends `openmp`,
   `stlpar`...) ;
3. pour une loi orthotrope, les gradients sont tournés dans le repère du
   matériau avant l'appel, les blocs de la tangente dans le repère global
   après, en une fois pour tous les points (fonctions de rotation de MGIS
   sur des vues ; fonction de rotation des lots pour une loi
   `generic-parallel`, qui n'a pas de fonction de rotation d'un point).

Contraintes dans le repère global pour le résidu : tournées point par
point, comme avant, quand la loi fournit la rotation d'un point ; pour une
loi orthotrope `generic-parallel`, tournées pour tous les points par
`integrateBatch`, juste après l'intégration. Les post-traitements (forces
moyennes, évaluateurs tournés) les calculent pour tous les points à chaque
appel.

Propriétés matériau de début de pas : pour l'appel, celles de fin de pas
leur sont substituées sous forme de vues (comme le chemin point par
point), puis les leurs sont remises en place, sans copie.

Refus explicites, avec message : `Batch` avec PETSc, intégrateurs sans
support des lots (endommagement micromorphe), loi `generic-parallel` point
par point.

**Relecture indépendante** : un défaut grave trouvé et reproduit. La
première version gardait les contraintes tournées dans un cache du
matériau, invalidé par des versions des tableaux ; MGIS écrit par ses vues
sans changer les versions, et l'assemblage du résidu avec PETSc (intégration
puis résidu, élément par élément) relisait des contraintes périmées pour
toute loi orthotrope. Cache et versions retirés (conception ci-dessus) ; la
séquence de la relecture est devenue un test. Corrigés aussi : copie des
propriétés matériau à chaque appel (vues), documentation, message de la
variable d'environnement, tests en 2D, FBar orthotrope, prédiction.
Documenté sans changement : MGIS donne à la loi un facteur de pas initial
de 10 (1 point par point), sans effet pour mfem-mgis.

Tests :

- mêmes résultats, au bit près, point par point et par lots : plasticité,
  grandes transformations, FBar, propriétés matériau de fin de pas
  variables en espace, lois orthotropes en petites et grandes
  transformations (nouvelle loi de test
  `OrthotropicSaintVenantKirchhoffElasticity`), rotation uniforme ou axes
  variables d'un point à l'autre, déformations planes, prédiction ;
- bibliothèques `generic-parallel` (`openmp`, `stlpar`) comparées au calcul
  point par point avec les rotations d'un point : au bit près (1e-13 avec
  prédiction en grandes transformations), ce qui vérifie les rotations de
  tous les points indépendamment ;
- échec de la loi en certains points ;
- résidu élément par élément juste après l'intégration de chaque élément ;
- suite de mfem-mgis : 139/139 ; dans le mode `Batch` (variable
  d'environnement) : 133/139, les 6 échecs sont les refus attendus
  (intégrateurs micromorphes) ; exemples : 18/20 en `Batch` (refus de
  PETSc) ; mm-opera-hpc : 6/6 en `Batch` ;
- cas du prototype : poutre et polycristal de mm-opera-hpc par lots avec
  leurs lois `generic-parallel`, comparés au calcul point par point.

**Deux défauts trouvés par les mesures, hors de mfem-mgis :**

- `scripts/env.sh` exporte `OMP_PROC_BIND=spread OMP_PLACES=threads`
  (pour Kokkos et AdaptiveCpp). libgomp, chargée par les bibliothèques des
  cas, attache alors le fil principal au premier processeur : le backend
  `stlpar` (TBB) ne voyait qu'un processeur et restait séquentiel, et
  `openmp` à 8 fils, répartis sur les fils matériels, n'occupait que 4
  cœurs. Le script de mesures retire ces variables.
- Les fichiers de configuration de `g++` du fork de TFEL compilaient les
  lois `generic-parallel` en `-O2` seul : à `-O2`, g++ 14 ne vectorise
  qu'avec un modèle de coût très restrictif. La loi du polycristal était
  17 % plus lente qu'avec l'interface `generic` (78,4 s contre 67,1 s
  d'intégration sur 2 pas). Corrigé (patch 0093 de TFEL : les options de
  `tfel-config --oflags`, dont `-ftree-vectorize`). Les autres fichiers de
  configuration (GPU) n'ont pas non plus `-DNDEBUG` : à vérifier en L6.

Mesures sur machine calme (`results/L2-polycristal-2026-10-09/`,
`results/L2-poutre-2026-10-09/`) : sur le polycristal, la loi pèse 70 % du
calcul sur un processus ; par lots, elle passe de 104,9 s à 17,8 s avec
`openmp` à 8 fils et à 13,9 s avec `stlpar` (calcul complet : 150 s puis
59 s). À 8 processus d'un fil, 14,5 s dans tous les modes. Sur la poutre,
la loi pèse moins de 5 % : gain invisible.

### 2026-10-09 (nuit), lot L3 commencé : gradients et forces intérieures par noyaux

Premières briques du chemin par lots sans passe élément par élément sur
l'hôte (`src/Batch/`, sources compilées comme du CUDA quand MFEM l'est) :

- **gradients** (`batch::computeUnknownsDerivatives`,
  `batch::computeGradients`) : restriction aux éléments (ordre natif),
  dérivées physiques aux points par l'interpolateur de MFEM
  (`QuadratureInterpolator::PhysDerivatives`, rangement `byVDIM`), une
  fois par règle d'intégration pour tous les matériaux, puis un noyau par
  matériau qui écrit les gradients de la loi (déformation avec √2, ou
  gradient de la transformation) dans le tableau du matériau, dans son
  rangement, gradients macroscopiques compris ;
- **forces intérieures** (`batch::computeElementsGeometry`,
  `batch::addInnerForces`) : gradients des fonctions de forme sur
  l'élément de référence (`DofToQuad`), jacobiens et déterminants
  (`GeometricFactors`), un fil par nœud d'élément, `Bᵀ σ w det J` avec le
  même opérateur `B` que les intégrateurs générés, écrit dans un vecteur
  par éléments au rangement de la restriction de MFEM ;
- interface des intégrateurs : `supportsKernels`,
  `computeGradientsByKernels`, `addInnerForcesByKernels` (base CRTP
  standard : mécanique 3D, déformations planes, contraintes planes,
  petites et grandes transformations, sans FBar, règle d'intégration
  commune à tous les éléments du matériau, maillage d'une seule
  géométrie) ;
- paramètre provisoire `GradientsComputation` (`ElementByElement` par
  défaut, `Kernels`), avec `BehaviourIntegration=Batch`, variable
  d'environnement `MFEM_MGIS_GRADIENTS_COMPUTATION`.

`QuadratureInterpolator::MultTranspose` n'est pas implémenté dans MFEM
4.10 : le noyau des forces intérieures est écrit ici.

Tests : gradients et forces intérieures des noyaux égaux à ceux du calcul
élément par élément à 1e-16 près (hexaèdres et tétraèdres, ordres 1 et 2,
petites et grandes transformations, orthotrope, déformations planes, deux
matériaux avec gradients macroscopiques), sur l'hôte, sur le périphérique
`debug` de MFEM (qui a trouvé une lecture hôte sans synchronisation dans
le test lui-même) et sur le GPU ; calculs complets avec noyaux égaux au
calcul point par point à 1e-12 près. Suite de mfem-mgis 142/142 ; avec
tous les noyaux activés par les variables d'environnement, seuls les refus
attendus échouent (thermique, endommagement micromorphe).

Étape suivante : l'opérateur propre (résidu sans la boucle de MFEM sur les
éléments, tangente sans matrice ou par matrices élémentaires calculées par
noyaux), puis la thermique et FBar.

### 2026-10-09 (nuit), lot L3 : tangente par noyaux, géométrie, relecture

- **Géométrie sans nœuds** : `Mesh::GetGeometricFactors` (et
  `QuadratureInterpolator::PhysDerivatives`, qui l'appelle) donne des
  nœuds au maillage (`EnsureNodes`). Toutes les transformations d'éléments
  passent ensuite par eux, et les boucles restées élément par élément
  étaient 50 % plus lentes (poutre). Les facteurs géométriques des noyaux
  sont calculés sur une copie des sommets (module `Batch/ElementsGeometry`,
  gardés par règle d'intégration), et les gradients partent des dérivées
  de référence, converties par l'inverse des jacobiens dans le noyau.
- **Matrices élémentaires de la tangente par noyaux**
  (`TangentOperatorComputation`) : un fil par couple de nœuds d'un
  élément, `Σ Bᵀ K B w det J`, recopiées dans l'assemblage de MFEM. Égales
  au calcul élément par élément à 1e-16 près (CPU, `debug`, GPU). Plus
  lentes sur CPU (+24 % par appel) : elles ne vaudront qu'avec un
  assemblage sur le périphérique, à écrire (`ElementRestriction::
  FillJAndData` de MFEM ne traite que les espaces scalaires), ou sans
  matrice.
- **Résidu par noyaux sans boucle de MFEM** quand le formulaire n'a pas
  d'autre intégrateur.
- **Relecture indépendante** des noyaux : conventions de MFEM, opérateur
  B, rangements vérifiés justes ; deux défauts graves reproduits et
  corrigés : sur le GPU, les conditions de Dirichlet écrivaient les
  inconnues sur l'hôte sans le déclarer à MFEM (copie périmée lue par le
  noyau des gradients, Newton qui échoue) ; les intégrateurs régularisés
  (`Faltus2026`) étaient acceptés et perdaient leur terme de
  régularisation sans message. Corrigés aussi : refus au `setup`
  (maillage mixte, espace non pris en charge, attribut sans intégrateur),
  exceptions du Newton rapportées au lieu d'arrêter le programme, accès
  hôte non déclarés de la prédiction (plantage sur le périphérique
  `debug`), cache de géométrie sans entrée partielle, options de
  l'interpolateur partagé de MFEM remises en état. Nouveaux tests :
  refus de `Faltus2026`, `BatchIntegrationTest` sur `debug` et `cuda`,
  poutre et polycristal (périodique) en mode noyaux sur 2 processus.

Mesures (`results/L3-poutre-2026-10-09/`, `results/L3-polycristal-2026-10-09/`) :
gradients par noyaux, phase d'intégration de la poutre divisée par 2,5 à
3,4 ; résidu par noyaux divisé par 2,1 (poutre) et 2,8 (polycristal), avant
la suppression de la boucle de MFEM ; matrices élémentaires par noyaux plus
lentes sur CPU (+24 à +38 %). Calculs complets avec tous les noyaux sur le
GPU égaux au calcul point par point (`BatchIntegrationTest-cuda`).

### 2026-10-09 (nuit), lot L3 : assemblage direct de la tangente, premiers calculs sur GPU

- **Assemblage direct** (mfem-mgis b47f65c3cf) : les blocs `N x N` des
  matrices élémentaires calculés par noyaux sont ajoutés directement dans
  la matrice creuse locale (position cherchée dans la ligne, addition
  atomique, entrées hors du motif signalées), sans matrices élémentaires
  stockées. Le motif vient du premier assemblage de MFEM ; la fin de
  `GetGradient` est celle de MFEM (bloc diagonal et `PtAP` en parallèle,
  `RAP` en séquentiel).
- **Élimination des ddl imposés** : `HypreParMatrix::EliminateRowsCols`
  recopie la matrice sur l'hôte et calcule la partie éliminée, inutile
  ici (122 ms par appel sur GPU pour 111 843 inconnues). Remplacée par
  `EliminateBC` (`DIAG_ONE`), en place et sur le GPU, même matrice
  puisque le motif est symétrique.
- **Premiers calculs complets sur GPU** (poutre r4, 3 pas, 1 processus,
  `--device cuda`) : 7,8 s avec tous les noyaux, contre 24,5 s sur un
  cœur de CPU et 29,1 s sur GPU sans noyaux. Le noyau d'assemblage prend
  7 ms par appel (0,48 s sur un cœur). hypre est compilé pour le GPU
  (résolution 0,9 s au lieu de 13,4 s). Restent sur l'hôte : la
  prédiction (opérateurs linéarisés assemblés par MFEM, 4,9 s), la loi
  (`generic`) et le premier appel de la tangente (boucle de MFEM pour le
  motif, 0,9 s).
- **Sur CPU**, l'assemblage direct ne gagne rien : un fil par paire de
  nœuds recalcule l'inverse du jacobien et relit la tangente pour chaque
  paire (0,48 s par appel, contre 0,34 s pour la boucle de MFEM).
- **Seconde relecture** des noyaux : défauts de sévérité 5 vérifiés
  corrigés. Traités : échec d'un noyau sur un seul processus rendu
  collectif avant les opérations collectives suivantes (sinon blocage),
  options de l'interpolateur restaurées même sur exception, vecteur des
  forces gardé entre deux résidus, conditions de Dirichlet sans ddl sur un
  processus sans transfert inutile. Un constat de la relecture était
  faux : avec le périphérique `debug` de MFEM, tout accès sur le
  périphérique, même en lecture, protège la copie hôte
  (`MemoryManager::GetDevicePtr`) ; la relecture hôte des ddl imposés
  reste nécessaire (plantage reproduit sans elle). Les relectures hôte
  des inconnues, elles, étaient inutiles (vérifié sur `debug` et `cuda`,
  séquentiel et 2 processus) et sont retirées.
- **Piège des tests** : les exécutables de test de mfem-mgis sont hors de
  `all` (cible `build-tests`), et la bibliothèque est statique ; sans
  reconstruire `build-tests`, ctest relance les anciens binaires. Il faut
  aussi sourcer `scripts/env.sh` (le `LD_LIBRARY_PATH` du profil charge un
  autre TFEL).

Mesures : `results/L3-poutre-assemblage-2026-10-09/`. Suite : motif
construit à partir des ddl des éléments (premier appel sans boucle de
MFEM, colonnes triées), variante du noyau adaptée au CPU (un fil par
élément), test MPI des noyaux dans mfem-mgis.

### 2026-10-09 (nuit), lot L3 : motif de la tangente construit à partir des éléments

- mfem-mgis 799941c8fd : la matrice locale est créée à partir des ddl des
  éléments au premier appel de `GetGradient` (plus de boucle de MFEM).
  Premier essai avec des colonnes triées : faux en parallèle.
  `HypreParMatrix` (constructeur bloc diagonal) partage les tableaux de
  la matrice locale et met le terme diagonal en tête de chaque ligne, en
  place, sans toucher au drapeau `isSorted`. Le motif suit donc l'ordre
  de hypre (diagonale, puis colonnes triées), que ce réordonnancement
  laisse intact, et la recherche du noyau accepte tout ordre (première
  colonne, dichotomie, puis parcours).
- Le défaut n'a été vu que par les cas du prototype et la suite en mode
  noyaux (tests `...P`) : aucun test de mfem-mgis sur les noyaux ne tourne
  en parallèle. Test MPI à ajouter.
- Tangente de la poutre (3 pas) : 16 ms par appel sur GPU, 0,40 s sur un
  cœur de CPU contre 0,45 s pour la boucle de MFEM. Sur GPU, la
  prédiction fait maintenant 73 % du temps.

### 2026-10-09 (soir), tranche GPU minimale : loi intégrée sur le GPU

- MGIS a0b5140 : `make_unchecked_view`, vue sur un gestionnaire sans
  vérifier que les valeurs qu'il garde sur l'hôte (valeurs uniformes,
  stockage local) sont accessibles aux calculs de la loi ; l'appelant
  remplace les vues des tableaux et copie ces valeurs. `make_view`
  inchangé (refus). Test CUDA étendu.
- mfem-mgis ee1f31fa1c : une loi `generic-parallel` est intégrée sur les
  copies des tableaux du matériau sur le périphérique si ses calculs y ont
  accès (loi `cuda` sur GPU ; toute loi sur le périphérique `debug`, dont
  les copies sont en mémoire hôte). Vue de MGIS déplacée sur les copies
  rendues valides par MFEM (`Read`, `ReadWrite`, `Write`), propriétés,
  variables externes et masse volumique copiées sur le périphérique à
  chaque appel, rotations des lois orthotropes sur le périphérique, ainsi
  que le passage au repère global quand la loi n'a pas accès à l'hôte.
  Plus de synchronisation sur l'hôte avant l'intégration quand tout est
  sur le périphérique ; synchronisation en lecture après l'intégration si
  les forces ou la tangente restent élément par élément.
- Validation : lois `openmp` et `stlpar` sur le périphérique `debug`
  (accès à des copies invalides détectés), identiques au point par point
  au dernier bit sans noyaux ; lois `cuda` sur GPU, en séquentiel et sur 2
  processus, avec et sans noyaux, orthotropes compris. Cas `poutre-gpu`.
- Poutre r4, 3 pas, GPU : 6,25 s. Intégration 0,24 s (0,47 s avec la loi
  `generic` sur l'hôte). La prédiction fait 77 % du temps : prochain
  chantier (lot L4).

### 2026-10-09 (soir), prédiction par noyaux, poutre en 4 s sur GPU

- mfem-mgis 7159b6de4c : avec les forces et la tangente par noyaux,
  l'opérateur de la prédiction est assemblé par noyaux dans une matrice
  locale gardée d'une prédiction à l'autre, et les forces intérieures au
  début du pas aussi : plus de formes bilinéaires élément par élément.
  Système formé comme `FormLinearSystem` de MFEM, élimination en place
  (`B = Pᵀ b - K X`, puis `B[e] = X[e]`). Les conditions aux limites qui
  contribuent au second membre (pression) restent assemblées par MFEM ;
  repli sur l'ancien chemin si elles contribuent à l'opérateur.
- Défaut trouvé en passant : la liste des ddl imposés de la prédiction
  enveloppait la mémoire d'un `std::vector` ; lue sur le périphérique
  `debug`, elle faisait planter le test une fois sur quelques-unes (selon
  la disposition de la mémoire). Copie allouée par MFEM.
- Poutre r4, 10 pas (`results/L6-poutre-2026-10-09/`) : **4,06 s sur un
  GPU** contre 17,54 s sur 8 cœurs (point par point) et 14,68 s sur 8
  cœurs avec les mêmes noyaux. Il reste la résolution de hypre (61 %) et
  la prédiction (21 %, dont sa résolution). Deux processus sur le même
  GPU : 28,6 s (résolution 9 fois plus lente, non analysé).

### 2026-10-09 (nuit), relecture indépendante des lots L3, L4 (prédiction) et L6 (tranche GPU)

Relecture adverse de mfem-mgis b47f65c3cf à 7159b6de4c et de MGIS a0b5140
(lecture complète, MFEM 4.10 comparé, programmes de comparaison sur cpu,
`debug` et `cuda`, 1 à 3 processus dont un sans élément, éléments
d'ordre 2, rangement `byVDIM`). **Aucun défaut de sévérité 3 à 5.**
Prédiction par noyaux, élimination en place, motif dans l'ordre de hypre
et synchronisations après l'intégration sur le périphérique vérifiés.

Traité (mfem-mgis 9398a9e33c, decb6fd2c1) :
- tests (sévérité 2) : la prédiction de chaque pas est enregistrée et
  comparée ; cas avec pression imposée (les forces intérieures au début
  du pas comptent ; leur oubli dans la prédiction par noyaux, que les
  tests ne voyaient pas, est maintenant détecté) ; deux matériaux, l'un
  intégré sur le périphérique, l'autre sur l'hôte (nouveau maillage
  `cube-two-materials.mesh` : `cube_2mat_per.mesh` est périodique, ses
  faces x = 0 et x = 1 sont identifiées et une condition de Dirichlet n'y
  charge rien) ; loi à propriétés variables d'un point à l'autre en
  `generic-parallel` ; noyaux avec prédiction pour ces lois ;
- documentation devenue fausse (calcul dans le repère global sur le
  périphérique, commentaires de `solve` et de `MaterialArrays`) ;
- copie inutile vers l'hôte des forces tournées à chaque résidu ; fuite
  possible de la matrice de la prédiction ; fonction de chargement des
  conditions de Dirichlet évaluée sur tous les processus ; débordement
  d'entier refusé ; tests GPU jamais simultanés.

Restent (sévérité 1, notés) : `update` et `revert` sur l'hôte (aller et
retour des tableaux à chaque pas ; MGIS a les opérations CUDA pour les
faire sur le GPU) ; synchronisation de tous les tableaux quand seuls les
gradients sont calculés sur l'hôte ; matrice de la prédiction en double
de la tangente ; chaque processus réserve environ 5 Go sur le GPU
(réserve Umpire de hypre), et hypre serait compilé sans MPI compatible
GPU, ce qui pourrait expliquer la lenteur de deux processus sur un GPU.

### 2026-10-09 (nuit), polycristal sur GPU

Polycristal de mm-opera-hpc (loi de monocristal de l'UO2 compilée pour
`cuda`, 10 pas, `results/L6-polycristal-2026-10-09/`) : **13,18 s sur un
GPU contre 31,65 s sur 8 cœurs** (151 s sur un cœur), mêmes contraintes
macroscopiques à 4e-13 près de l'échelle des contraintes. La loi, qui
domine sur CPU, prend 5,7 s sur le GPU, autant que la résolution. Avec la
loi sur l'hôte, le GPU ne gagne rien (107 s). Deux processus sur le même
GPU : résolution 14 fois plus lente, comme pour la poutre.

### 2026-10-09 (nuit), fin de pas sur le GPU

mfem-mgis (commit « Update and revert on the device ») : les états d'un
matériau dont la loi est intégrée sur le périphérique sont mis à jour et
restaurés sur le périphérique (`update` et `revert` de MGIS sur une vue
des copies, par `MFEMFieldOperations`), et la synchronisation sur l'hôte
au début de chaque résolution n'est plus faite quand toutes les lois sont
intégrées sur le périphérique (sauf avec PETSc). Les tableaux restent sur
le GPU d'un pas à l'autre. Poutre r4, 10 pas, GPU : 3,3 s au lieu de
4,1 s (8 cœurs : 17,5 s).

### 2026-10-09 (nuit), profil Nsight Systems de la poutre sur GPU

Pas de copie cachée notable entre l'hôte et le GPU (146 Mo vers le GPU
sur 3 pas, surtout des données de mise en place ; un ou deux vecteurs par
itération). Le GPU calcule 0,71 s sur 1,8 s, en 95 000 lancements de
noyaux, dont 65 % de produits matrice-vecteur de hypre : la résolution est
limitée par la latence des lancements. Détails dans
`results/L6-poutre-2026-10-09/README.md`.

### 2026-10-10 (nuit), action de la tangente sans assemblage

- mfem-mgis cbe693a9a8 : action de la tangente `Σ Bᵀ K B x w det J` et sa
  diagonale par noyaux, sans assemblage (incréments des gradients par le
  noyau des gradients, produit par les blocs de la tangente à chaque
  point, noyau des forces intérieures ; diagonale par le noyau des
  matrices élémentaires restreint aux couples d'un nœud avec lui-même).
  Égales au produit par la matrice assemblée et à sa diagonale à 1e-16
  près (tous les cas de `GradientsKernelsTest`, cpu, `debug`, `cuda`).
  Contrôles communs des calculs par noyaux factorisés (`getKernelsData`).
- mfem-mgis 4d465460c3 : `TangentOperatorComputation=MatrixFree`, la
  tangente rendue par `GetGradient` est un opérateur sans matrice
  (`MatrixFreeTangentOperator`, ddl imposés par `ConstrainedOperator`),
  avec le préconditionneur `OperatorJacobiSmoother` (diagonale par
  `AssembleDiagonal`, ou `GetDiag` pour une `SparseMatrix` de MFEM, qui
  n'a pas `AssembleDiagonal`). Calculs complets égaux au point par point
  (petites et grandes déformations, orthotrope, deux matériaux, lois
  `generic-parallel`, `debug`, GPU, 2 processus). Cas `poutre-sans-matrice`.
- Poutre r4 sur GPU : 95 s avec CG + Jacobi contre 3,4 s avec la tangente
  assemblée et BoomerAMG. Juste, mais il faut un préconditionneur adapté
  (AMG gardé plusieurs itérations, multigrille polynomial).

### 2026-10-10 (nuit), préconditionneur sur la tangente assemblée

- mfem-mgis b129009b0c : préconditionneur `AssembledOperator` pour la
  tangente sans matrice. Il demande à l'opérateur sa version assemblée
  (tangente assemblée par noyaux, ddl imposés éliminés) et y construit le
  préconditionneur donné en option (BoomerAMG, Jacobi...), refait toutes
  les `UpdatePeriod` itérations de Newton et gardé entre deux. Solveur
  `cg-amg` du cas de la poutre.

### 2026-10-10 (nuit), corrections de la seconde relecture des lots L3 et L6

Relecture indépendante de cbe693a9a8 et 4d465460c3
(`scratch/nuit-2026-10-09/relecture-L3-L6-suite.md`). Corrigé dans
mfem-mgis 1bbcb877d8 (et MGIS 860b693) :

- F1 : `update` et `revert` d'une loi intégrée sur le périphérique ne
  reportaient plus les variables externes de s1 dans s0 (température de
  s0 figée : dilatation thermique recomptée à chaque pas, sans message).
  Nouvelle fonction de MGIS `updateInputValues` (propriétés matériau,
  variables externes et masse volumique, mêmes règles que `updateValues`,
  sans les tableaux), appelée sur l'hôte après la copie des tableaux sur
  le périphérique.
- F2 : après `update` sur le périphérique, les copies hôte des tableaux
  étaient périmées. Elles sont maintenant synchronisées en lecture à la
  fin de `update` et de `revert` (une copie du GPU vers l'hôte par pas).
- F3 : une exception de la tangente, de sa diagonale ou de son action dans
  `computeNewtonCorrection` (`noexcept`) finissait en `std::terminate`.
  Elle est enregistrée dans le contexte de la résolution, qui échoue
  proprement.
- F4 : diagonale sur maillage non conforme : `|Pᵀ| diag(A)`, comme MFEM.
- F5, F6, A2, Q2 : copie de la liste des ddl imposés de l'opérateur
  contraint (reconstruit quand elle change), liste du lisseur de Jacobi
  gardée en membre (lisseur réutilisé), garde du nombre de fils par nœud
  pour la diagonale, fonctions partagées au lieu de copiées, paramètre
  `TangentOperatorComputation` lu une fois (énumération), documentation
  des synchronisations.
- Tests (T1) dans `BatchIntegrationTest` : loi avec dilatation thermique
  et température variable dans le temps, états relus sur l'hôte sans
  synchronisation après chaque `update` et `revert`, diagonale de
  l'opérateur sans matrice comparée à celle de la tangente assemblée,
  refus propre de MatrixFree avec un intégrateur de bord. Chaque test a
  été vérifié par mutation : sans l'appel de `updateInputValues`, les cas
  « temperature » échouent sur `debug` ; sans la synchronisation, `debug`
  s'arrête à la lecture ; avec une diagonale faussée de 1e-6, tous les
  cas MatrixFree échouent ; sans le `try/catch`, `std::terminate`.
- Suites : 146/146 (dont 7 sur GPU), noyaux forcés 129/139 (10 refus
  attendus : chaleur, micromorphe), cas 19/19.

### 2026-10-10 (nuit), noyaux des dérivées et action fusionnée

- mfem-mgis a2602d117a : les dérivées des inconnues aux points sont
  calculées par un noyau propre (gradients des fonctions de forme de
  `mfem::DofToQuad`) au lieu de `QuadratureInterpolator::Derivatives` de
  MFEM, dont la version générique pour les simplexes coûtait 890 µs par
  appel sur la poutre (39 µs maintenant). L'action de la tangente sans
  matrice est calculée par un seul noyau par matériau, un fil par élément
  (gradient de l'incrément, produit par K, `Bᵀ`, sans tableau aux points) :
  128 µs par matériau au lieu de 350 µs. Temps de calcul du GPU de la
  poutre sans matrice (10 pas) : 2,88 s au lieu de 5,15 s.
- Le noyau fusionné lit les blocs de K point par point (accès non
  coalescents, environ 45 % de la bande passante) : le rangement
  composante par composante (paramètre `QuadratureDataLayout`, prévu en
  L1) reste à brancher dans les intégrateurs.

### 2026-10-10 (nuit), solveur de point fixe (lot L5)

mfem-mgis f3d8ed1e55 : solveur non linéaire `FixedPoint`, tiré du
prototype du 04/10 (`mfem-mgis-gpu`) et écrit dans mfem-mgis :

- méthodes `ElasticOperator` (K0 assemblé une fois par résolution ou une
  fois pour tout le calcul, résolu par le solveur du problème ou par un
  solveur propre dont la mise en place est gardée) et `ModifiedNewton` ;
- accélérations dans une fabrique extensible (`AccelerationAlgorithmFactory`)
  : `Anderson` (type II, matrice des produits scalaires mise à jour au fil
  des itérations, deux produits scalaires par itération et par différence
  gardée) et `Cast3M` (comme MTest) ;
- `HypreBoomerAMG` déclaré aussi comme solveur linéaire (un cycle d'AMG
  par itération).

Cube en traction (plasticité à écrouissage linéaire), 2 pas : Newton 8
itérations, point fixe 103 (taux observé 0,642, borne 3μ/(3μ+H) = 0,887),
Anderson 37, Cast3M 56, Newton modifié 59, un cycle d'AMG + Anderson 46.
Mêmes solutions que Newton à 1e-12 près, sur CPU, `debug` et `cuda`, à 1
et 2 processus. Sur le GPU, un cycle de BoomerAMG est bien plus faible (101
itérations au lieu de 46) : réglages de hypre sur GPU (lisseur l1-Jacobi,
grossissement PMIS) à exposer.

Relecture indépendante de b129009b0c et 1bbcb877d8 terminée
(`scratch/nuit-2026-10-10/relecture-1bbcb877d8.md`) : F1 à F6 confirmés
corrigés ; nouveaux constats C3 (écriture directe après `update` ignorée
sur le périphérique), C7 (exception sur un seul processus : blocage), C4
et C5 (matrice du préconditionneur `AssembledOperator`), C9 (tests), C8
(copie de K). Corrections en cours.

### 2026-10-10 (nuit), corrections de la relecture de 1bbcb877d8

mfem-mgis bae47cb272 :

- C3 : `update` et `revert` d'une loi intégrée sur le périphérique se
  terminent par la synchronisation des tableaux d'état sur l'hôte (copies
  du périphérique invalidées), comme sur l'hôte : une écriture directe
  dans l'état après `update` est prise en compte. C8 : les blocs de K,
  nuls, sont mis à zéro sur l'hôte au lieu d'être recopiés.
- C7 : un échec de l'assemblage de la tangente ou de la mise en place du
  solveur linéaire sur un processus est rendu collectif avant la
  résolution ; pendant la résolution, seuls les échecs signalés sur tous
  les processus (`CollectiveFailure`, levée par
  `raiseOnFailureOnAnyProcess`) sont rattrapés en parallèle.
- C4 et C5 : le préconditionneur `AssembledOperator` connaît le nombre
  d'assemblages de la tangente (compteur partagé par le problème et ses
  fournisseurs) : il se remet en place si sa matrice a été réassemblée
  par un autre appel, et garde une mise en place valide pour la
  résolution de la prédiction au lieu de remettre la période à zéro.
  Poutre r3, 4 processus, 10 pas, 25 itérations de Newton : 25, 13 et 7
  mises en place de BoomerAMG pour une période de 1, 2 et 4 (avant : une
  de plus par pas, quelle que soit la période).
- Tests : écriture directe après `update` (détecte C3 par mutation sur
  `debug`) ; diagonale de l'opérateur sans matrice vérifiée aussi avec le
  préconditionneur assemblé, et cas sans prédiction (détecte C4 par
  mutation : arrêt dans `MPI_Comm_size` à 2 processus).
- Suites : 141/141, GPU 7/7, noyaux forcés 131/141 (10 refus attendus),
  cas 19/19.
- Restent (sévérité 1) : C11 (en séquentiel, AssembledOperator n'enveloppe
  que Jacobi), C12 (préconditionneur gardé après un échec de Newton), C14
  (une réduction MPI par application de l'opérateur sans matrice), tests
  de F4, F5, F6.

### 2026-10-10 (nuit), accélérations de MTest

mfem-mgis : portage de toutes les accélérations de MTest dans la fabrique
du point fixe : `UAnderson` et `FAnderson` (méthode d'Anderson sur les
résidus du point fixe ou sur ceux du problème, accélération tous les deux
pas d'itération par défaut), `Secant`, `CrossedSecant`, `AlternateSecant`,
`IronsTuck`, `Steffensen`, `CrossedDelta2`, `Crossed2Delta`,
`Crossed2Deltabis`, `AlternateDelta2`, `Alternate2Delta`. Cube en traction
(2 pas, Newton 8 itérations, point fixe seul 103) : Anderson 37,
UAnderson et FAnderson 38, CrossedSecant 41, Alternate2Delta 46, IronsTuck
48, Cast3M 56, AlternateDelta2 60, AlternateSecant 61, Secant 65,
Crossed2Delta 70, CrossedDelta2 140, Crossed2Deltabis 374 ; Steffensen,
appliqué composante par composante, diverge sur ce problème couplé (même
formule que MTest). Toutes les autres donnent la solution de Newton à
1e-12 près.

### 2026-10-10 (nuit), options de BoomerAMG, prédiction par K0, relecture L5

- mfem-mgis d05475d265 : options de BoomerAMG (`StrengthThreshold`,
  `CoarseningType`, `AggressiveCoarseningLevels`, `InterpolationType`,
  `RelaxationType`, `RelaxationSweeps`, `CycleType`, `NumberOfCycles`,
  `MaximumNumberOfLevels`), réappliquées après chaque `SetOperator` (qui ne
  garde qu'une partie des réglages de hypre). Sur le GPU, le lisseur de
  Chebyshev (`RelaxationType=16`) rend un cycle de BoomerAMG bien plus
  efficace pour le point fixe (438 itérations au lieu de 1 231).
- Mesures L5 (`results/L5-poutre-2026-10-10-0250/`) : sur 8 cœurs, le
  point fixe (Anderson 8, un cycle d'AMG, K0 gardé) bat Newton (13,3 s
  contre 20,3 s sur machine chargée ; 9,4 s contre 16,6 s sur machine
  calme, ×1,7) ; sur le GPU, égalité (3,4 s) avec Chebyshev.
- mfem-mgis 846551e0cb : `getElasticPredictionLinearSolver` : le point fixe
  peut fournir le solveur de la prédiction élastique (gradient conjugué
  préconditionné par l'AMG de K0 déjà construite) ; gain modeste
  (prédiction 1,34 s au lieu de 1,52 s sur le GPU), la prédiction reste
  le poste principal du point fixe.
- Relecture indépendante de a2602d117a, f3d8ed1e55, bae47cb272
  (`scratch/nuit-2026-10-10/relecture-L5.md`) : noyaux et point fixe
  justes (42 configurations par périphérique, Anderson identique à une
  implantation directe, Cast3M fidèle à MTest) ; C7 pas corrigé (l'accord
  collectif pouvait s'apparier avec une autre réduction). Corrigé dans
  1d5adf1a32 : en parallèle, seules les `CollectiveFailure` sont traitées,
  toute autre exception arrête le programme ; options du point fixe
  validées avant d'être appliquées, actions appelées une fois, garde NaN,
  informations du solveur linéaire, bornes des nombres d'itérations dans
  le test (un Cast3M faux est détecté), variantes noyaux et `debug`.

### 2026-10-10 (nuit), NewtonK0, forces intérieures, relecture courte

- mfem-mgis 3e414bc352 : forces intérieures par un fil par élément (81 µs
  au lieu de 217 µs par appel sur le GPU).
- mfem-mgis 1ef05104d1 : sections de profilage de la prédiction
  (`PRED::integrate`, `forms`, `assemble`, `solve`) ; les mises en place
  ponctuelles (géométrie 214 ms, motif) faussaient les moyennes par appel.
- mfem-mgis 30beae4764 : méthode `NewtonK0` du point fixe (Newton,
  gradient conjugué préconditionné par BoomerAMG de K0 construit une
  fois) ; corrections de la relecture courte de 846551e0cb et 1d5adf1a32
  (`scratch/nuit-2026-10-10/relecture-L5-suite.md`) : solveur de prédiction
  du point fixe limité à la prédiction élastique en début de pas (N1),
  message puis `MPI_Abort` pour une exception non collective (N4), K0 et
  son solveur remplacés sans état intermédiaire (N2), test du solveur de
  prédiction par un solveur linéaire qui compte ses résolutions (mutation
  détectée).
- Mesures (`results/L5-poutre-2026-10-10-0250/README.md`) : sur le GPU,
  NewtonK0 2,6 s contre 3,16 s pour Newton (stable) ; sur 8 cœurs, point
  fixe accéléré 9,5 s, NewtonK0 11,6 s, Newton 13,5 s.

### 2026-10-10 (nuit), NewtonK0 sur le polycristal, relecture K0

- mfem-mgis 9249013150 : `PeriodicNonLinearEvolutionProblem` prend les
  paramètres du problème (par exemple `NonLinearSolver`). mm-opera-hpc
  6fdd03c puis 155202f : option `--non-linear-solver` du polycristal
  (`Newton` ou `NewtonK0`, K0 obtenu par une intégration de la loi, qui ne
  fournit pas d'opérateur de prédiction).
- Mesures (`results/L5-polycristal-2026-10-10-0340/README.md`, script
  `cas/mesures-L5-polycristal.sh`) : sur 8 cœurs, NewtonK0 24,9 s contre
  29,6 s (×1,19, résolution -36 %) ; sur le GPU, 11,3 s contre 11,9 s
  (×1,06), la loi domine. Mêmes itérations de Newton, contraintes égales à
  2e-13 près.
- Relecture courte de 30beae4764 (`scratch/nuit-2026-10-10/relecture-K0.md`,
  36 configurations NewtonK0 en séquentiel et sur 2 processus, à 5e-11 de
  Newton) : rien au-dessus de la sévérité 2. Corrigé dans 9355256a9f :
  `TangentLinearSolver` et `PredictionLinearSolver` nomment un solveur
  itératif de MFEM de la fabrique (gradient conjugué par défaut ; GMRES,
  FGMRES, BiCGSTAB, MINRES, SLI), le gradient conjugué demande une tangente
  symétrique définie positive (documenté) ; `None` pour `Acceleration` ;
  `FGMRESSolver` dans la fabrique ; régularisation d'Anderson bornée à
  1e-14 ; tests : chaque solveur itératif, tangente sans assemblage,
  options changées entre les pas, prédiction à gradients constants ; six
  mutations détectées (`scratch/nuit-2026-10-10/mutations-K0/`). Option
  `-tls` de la poutre (solveur itératif de NewtonK0).

### 2026-10-10 (nuit), Newton inexact, relecture K0b

- Relecture de 9355256a9f (`scratch/nuit-2026-10-10/relecture-K0b.md`) :
  rien au-dessus de la sévérité 2. Corrigé dans 82087f9f83 : prédiction
  séquentielle sans noyaux avec la diagonale à un aux ddl imposés (un
  solveur itératif préconditionné par K0 la résout en une itération :
  GMRES faisait 2 21 itérations de point fixe, 2 2 maintenant), options de
  base appliquées seulement si toutes sont valides, messages, tests
  (remplacement du solveur de K0 entre deux pas avec un solveur qui compte
  les appels après sa destruction, solveur et tolérance par défaut de
  NewtonK0, prédiction par GMRES) ; toutes les mutations détectées.
- mfem-mgis 72f2e6e29b et 517c0c0b08 : Newton inexact, option
  `AdaptiveLinearSolverTolerance` du Newton et de la méthode `NewtonK0`
  (choix 1 et 2 d'Eisenstat et Walker, réglages de MFEM, garde de Kelley
  contre la sur-résolution, tolérance du solveur linéaire restaurée après
  chaque résolution, norme du résidu linéaire dans les informations
  d'itération) ; test unitaire `AdaptiveLinearSolverToleranceTest`.
- Mesures (`results/L5-poutre-2026-10-10-0429/README.md`, script
  `cas/mesures-L5-ew.sh`, options `-alst` et `-alsti` de la poutre) : sur
  le GPU, NewtonK0 avec Eisenstat-Walker 1 et une tolérance initiale de
  0,1 prend 1,96 s contre 3,12 s pour Newton (×1,59) et 2,57 s pour
  NewtonK0 à tolérance fixe ; sur 8 cœurs, 9,46 s, au niveau du meilleur
  point fixe. Avec Newton classique, Eisenstat-Walker ne paie pas (chaque
  itération de plus reconstruit BoomerAMG).
- Relecture du Newton inexact (`scratch/nuit-2026-10-10/relecture-EW.md`) :
  rien au-dessus de la sévérité 3 (couverture). Corrigé dans 4b72095ee0 :
  garde du choix 1 d'Eisenstat et Walker comme dans l'article et PETSc,
  tolérance initiale bornée par la maximale, tests qui recalculent la
  suite des tolérances à partir des informations d'itération (dix
  mutations détectées), section du glossaire.
- mfem-mgis b4127a6c12 : lisseur de Chebyshev de MFEM
  (`OperatorChebyshevSmoother`, ordres 1 à 5) dans la fabrique, test
  `OperatorSmoothersTest` (Laplacien : 46 itérations de CG avec Jacobi, 27
  et 17 avec Chebyshev d'ordre 2 et 4). Mesures sans matrice sur le GPU
  (`results/L5-poutre-2026-10-10-0456/README.md`, option `-ls cg-chebyshev`
  de la poutre) : Chebyshev plus lent que Jacobi (chaque ordre ajoute une
  action sans matrice) ; NewtonK0 + Eisenstat-Walker 1 sans matrice à 3 %
  de la tangente assemblée (2,05 s contre 1,99 s).
- Polycristal (mm-opera-hpc, options `--adaptive-linear-solver-tolerance`
  et `--adaptive-linear-solver-initial-tolerance`), mesures
  `results/L5-polycristal-2026-10-10-0502/README.md` : NewtonK0 avec
  Eisenstat-Walker 1 est le meilleur sur le GPU (9,58 s contre 11,88 s pour
  Newton, ×1,24) et sur 8 cœurs (23,46 s contre 29,60 s) ; la résolution
  devient négligeable sur le GPU et la loi domine (7,1 s sur 9,6 s).

### 2026-10-10 (matin), thermique par noyaux, relecture C

- Relecture de b4127a6c12, 4b72095ee0 et du polycristal
  (`scratch/nuit-2026-10-10/relecture-C.md`) : rien au-dessus de la
  sévérité 2. Corrigé : mfem-mgis 6cae30af99 (invariance des lisseurs par
  mise à l'échelle diagonale, qui détecte une diagonale ignorée ou un
  lisseur figé ; choix 1 sans résidu linéaire ; glossaire), mm-opera-hpc
  976a60b (erreurs des options), README des mesures (mémoire non réduite
  sans matrice, chiffres arrondis).
- mfem-mgis 0a26a38be0 : thermique stationnaire non linéaire par noyaux
  (lot L7.3) : genres de gradient de température 2D et 3D, nombre de
  composantes des inconnues séparé de la dimension dans tous les noyaux,
  bloc dérivée du flux par rapport à la température (matrice, diagonale,
  action sans matrice), valeurs des fonctions de forme, température aux
  points recopiée dans la variable externe (nouvelle méthode des
  intégrateurs thermiques générés et du générateur). Test
  `HeatTransferKernelsTest` : 3D et déformations planes, assemblée et sans
  matrice, avec et sans prédiction, lois generic et generic-parallel,
  séquentiel, 2 processus, `debug` et GPU, écarts au point par point de
  1e-16 à 3e-14 ; trois mutations détectées. Suite à noyaux forcés : seuls
  les 6 tests micromorphes restent refusés (10 avant). La tangente
  thermique est définie négative : le gradient conjugué ne s'applique pas
  sans matrice (GMRES dans le test). Relecture indépendante en cours.


### 2026-10-10 (matin), mesures thermiques, jacobiens par noyau

- Mesures de la thermique stationnaire par noyaux (`cas/thermique`, r5,
  279 873 inconnues, script `cas/mesures-L7-thermique.sh`,
  `results/L7-thermique-2026-10-10-0542/README.md`) : sur le GPU, la
  première intégration coûtait 1,76 s, calcul unique de la géométrie sur
  l'hôte (espace d'éléments finis d'ordre un : 1,0 s ; projection des
  coordonnées : 0,83 s ; facteurs géométriques de MFEM sur le GPU :
  0,02 s).
- mfem-mgis 8b0b47acfb : jacobiens des maillages sans nœuds calculés par
  un noyau à partir des sommets des éléments et de l'élément fini des
  transformations de MFEM (`GetTransformationFEforElementType`). Test
  `ElementsGeometryTest` : segments, triangles, quadrangles, tétraèdres,
  hexaèdres, prismes, pyramides, réguliers et déformés, inconnues d'ordre
  un et deux, maillage courbe, hôte, `debug` et GPU, écarts aux
  transformations de MFEM au plus 2e-15 ; six mutations détectées.
  Thermique sur le GPU : 3,01 s vers 1,25 s, ×3,8 contre 8 processus avec
  les noyaux (`results/L7-thermique-2026-10-10-0606/README.md`). Poutre
  r4 : NewtonK0 + Eisenstat-Walker 1 de 2,14 s à 1,90 s sur le GPU
  (`results/L5-poutre-2026-10-10-0610/README.md`) ; polycristal inchangé
  (`results/L5-polycristal-2026-10-10-0612/README.md`). Relecture
  indépendante en cours.
