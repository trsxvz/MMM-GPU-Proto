# Journal : mfem-mgis-gpu-proto

## ÉTAT

- FAIT : fork de TFEL (`tfel/`, branche `generic-parallel-behaviours`,
  partie de master `fec0ade99`, 01/10/2026), compilation installée dans
  `install/tfel` (`build-tfel.sh`, réinstallation : `make -j16 install` dans
  `build/tfel`, 10 à 60 s).
- FAIT : jalon 1, code généré compilable pour le GPU.
- FAIT : jalon 2, interface `generic-parallel` des lois (petites
  déformations et lois générales), backends `cuda`, `hip` (non testé, pas de
  GPU AMD) et `stlpar`.
- FAIT : jalon 3, tests `tfel-check` installés
  (`test-suite/behaviours/generic-parallel/{cuda-nvcc,stlpar-parunseq-gcc}`) :
  J2 (3D, déformations planes), Norton implicite (3D, contraintes planes),
  Hosford (valeurs propres) et StrainHardeningCreep3 (loi de matériau sur
  GPU, grandeurs) comparés à l'interface `generic` (écart ≤ 6e-13), deux
  rangements mémoire, trois types de calcul ; les deux configurations
  passent (CUDA 1 min 31, CPU 32 s).
- Balayage des 226 lois de test de MFront (`scripts/compile-sweep.sh`,
  résultats `results/compile-sweep-cuda-2026-10-05.tsv`), passage final :
  **160 compilent pour le GPU sans aucun appel hôte** (99 au premier
  passage) ; 33 grandes déformations (non supportées) ; 16 bloquées par du
  code utilisateur (`throw`, `std::cout`, `static const` calculé) ; 11 avec
  appels hôte restants : systèmes de glissement des monocristaux (6),
  code utilisateur (`raise_if`, modèles externes : 4), conversion
  d'opérateur tangent en grandes déformations (1) ; 6 hors sujet
  (fragments, modèle, erreur de fichier).
- FAIT : suite de tests de TFEL complète, 6963/6963 (pas de régression
  CPU).
- Commits locaux (7) dans `tfel/`, série de patchs dans `tfel/proto-patches/`.
- FAIT : backend HIP testé sur le GPU NVIDIA (hipcc 5.7 + CUDA 12.9).
- FAIT : lois de matériau sur GPU, valeurs propres (8 solveurs), pile CUDA,
  doc `docs/web/generic-parallel-behaviour-interface.md`.
- RESTE : systèmes de glissement (monocristaux) ; jalon 4 (TDLS) ;
  jalon 5 (grandes déformations) ; puis MGIS et mfem-mgis.

## OBJECTIFS ET PRINCIPES

- Portabilité CPU et GPU propre, flexible et idiomatique de tout
  TFEL/MFront (propriétés matériau, lois, modèles...), puis MGIS et
  mfem-mgis, avec pour but final une infrastructure complète
  mfem-mgis / MGIS / TFEL-MFront sur GPU.
- Rester dans l'esprit du code existant, mais corriger directement les
  défauts évidents, propres et sans casse, y compris dans le code existant.
- Pas d'optimisation fine GPU pour l'instant (pile, launch bounds...) :
  tout doit compiler et les performances ne doivent pas être mauvaises.
- Couvrir tout TFEL/MFront (grandes déformations, monocristaux, fonctions
  d'initialisation, post-traitements, rotations, modèles autonomes, etc.),
  sans code spaghetti.
- Tout doit fonctionner de -O1 à -O3, sur tous les backends : stdpar (gcc
  et nvc++), CUDA, HIP, Kokkos, SYCL ; OpenMP (hôte et déport) en plus.
- Organisation : un dossier par bibliothèque adaptée (`tfel/`, `mgis/`,
  `mfem-mgis/`) ; la série de patchs de TFEL depuis master est dans
  `patches/tfel/` (lien `tfel/proto-patches`).
- Tenir `JOURNAL.md` et `DESIGN.md` à jour au fur et à mesure : ils
  permettent de reprendre le travail.

## PLAN DE TRAVAIL (cocher au fur et à mesure)

A. Robustesse
- [x] A1 Suite de tests TFEL complète (05/10, après les 5 commits) :
      **6963/6963 réussis** avec un environnement propre
      (`env -u LD_LIBRARY_PATH`, PATH sans spack) et régénération forcée
      du code généré (`touch build/tfel/mfront/src/mfront`). Les échecs
      BarBehaviour2/4 et TDLS des passes précédentes venaient de l'ancienne
      libTFELMFront de spack. Commande : voir l'entrée du 05/10
      « réorganisation ».
- [x] A2 Interface `generic-parallel` des propriétés matériau (commit
      0006) : `stlpar` vraiment parallèle (itérateur de comptage partagé,
      `MFront/GenericParallel/CountingIterator.hxx`) ; `errno` remplacé par
      la détection des résultats invalides dans les noyaux (statut -4 :
      `errno` est par thread, absent sur GPU, et modifié par TBB/CUDA à
      l'initialisation, d'où les faux échecs du premier appel, cas UO2) ;
      bornes de la sortie écrites hors du tableau des statuts et jamais
      analysées ; en-têtes des backends installés à plat (les deux
      `AbstractBackend.hxx` s'écrasaient). Tests `stlpar-parunseq-gcc` et
      `cuda-nvcc` des propriétés matériau : OK.
B. Couverture MFront
- [~] B1 Tests numériques : Runge-Kutta, lois générales non mécaniques,
      zones cohésives, énergies, vitesse du son, politiques de bornes :
      faits (entrée du 05/10 nuit), matrice à repasser.
- [~] B2 Fonctions d'initialisation et post-traitements : faits, 118/118
      comparaisons stlpar-gcc ; matrice complète des backends à repasser.
- [~] B3 Rotations des lois orthotropes (sur les champs du lot) : faites
      (module partagé avec l'interface `generic`, vecteurs compris),
      stlpar-gcc OK, matrice à repasser.
- [~] B4 Monocristaux (systèmes de glissement, matrices d'interaction) :
      faits (générés et classes TFEL), CUDA et stlpar OK, matrice à
      repasser. Restent les polycristaux à orientations lues dans un
      fichier (`PolyCrystalsSlidingSystems`, entrées/sorties hôte).
- [~] B5 Grandes déformations (Green-Lagrange, Hencky, lois en F),
      conversions d'opérateurs tangents sur GPU : faits ; CUDA, HIP,
      Kokkos, SYCL GPU OK ; nvc++ corrigé (voir entrée du 05/10 soir),
      matrice complète à repasser.
- [~] B6 Modèles autonomes (`@DSL Model`) : interface `generic-parallel`
      des modèles faite (entrée du 05/10 nuit, B6) ; stlpar et CUDA OK,
      matrice à repasser.
- [~] B7 TDLS dans les lois implicites sur GPU : fonctionne sans
      changement de code (entrée du 05/10 nuit, B7) ; tests `tdls` par
      configuration ajoutés, matrice à repasser.
C. Backends
- [x] C1 HIP sur GPU NVIDIA (commit 0007) : configuration `hip-hipcc`
      ajoutée à `build-tfel.sh` (`HIP_PLATFORM=nvidia`, CUDA 12.9 du SDK
      HPC via `CUDA_PATH`) ; test `behaviours/generic-parallel/hip-hipcc`
      36/36 comparaisons OK (87 s), propriétés matériau `hip-hipcc` OK ;
      cuda-nvcc et stlpar-parunseq-gcc (lois et propriétés) relancés : OK.
- [x] C2 `stlpar` sur GPU avec `nvc++ -stdpar=gpu` : mode
      `-gpu=mem:managed` (le mode `unified` par défaut est mal compilé par
      nvc++ 25.11), option `device_stack_size` (16 Ko dans la
      configuration `stlpar-parunseq-nvhpc-gpu`) ; tests lois (646 s, la
      compilation nvc++ est lente) et propriétés matériau OK.
- [x] C3 Kokkos 5.0.2 (partagé, sans RDC, CUDA sm_120 + OpenMP + Serial,
      `install/kokkos`) : backend `kokkos` (option `execution_space`),
      tests `kokkos` (CUDA, 151 s) et `kokkos-host` (OpenMP, 71 s) OK.
- [x] C4 SYCL : AdaptiveCpp 25.10 compilé depuis les sources
      (`deps/adaptivecpp-25.10.0`, `build/adaptivecpp`,
      `install/adaptivecpp`) contre un LLVM 20.1.8 construit par spack
      (`deps/spack-sycl/llvm.spec`) ; backend `sycl` (mises à jour
      atomiques, transferts automatiques, `device_stack_size`) ; tests
      `sycl-default-acpp` (GPU) et `sycl-default-acpp-omp` (CPU) OK.
- [~] C5 OpenMP (backend `openmp`, hôte et déport `target`) : fait (entrées
      du 05/10 soir, C5) ; `openmp-host-gcc` activée (284/284) ;
      `openmp-offload-nvhpc-gpu` non activée (deux défauts de nvc++ 25.11,
      lois et propriétés matériau justes une à une).
- Matrice complète (`scripts/run-generic-parallel-tests.sh`,
  `results/generic-parallel-tests-2026-10-05.txt`) : 8 tests de lois et
  6 de propriétés matériau actifs, tous OK après les corrections ci-dessous
  (dernière passe en cours au moment de l'écriture).
D. Intégration
- [ ] D1 MGIS (proposition de conception détaillée :
      `scratch/d1-mgis-design.md`, 10 questions ouvertes) :
      chargement d'une bibliothèque `generic-parallel`,
      intégration d'un lot décrit par des vues. Point de conception : MGIS
      range chaque propriété matériau et variable externe séparément
      (uniforme ou champ, mémoire locale ou externe), alors que
      `mfront_gpb_BatchData` attend un champ unique à pas fixes par groupe.
      Choix à faire : champs par variable dans `BatchData` (tableau de
      `ConstField`, le plus fidèle à MGIS et sans copie) ou assemblage d'un
      bloc par MGIS avant l'appel.
- [ ] D2 mfem-mgis : assemblage GPU, variantes (Newton, points fixes,
      Ginkgo, hypre GPU). Constat (lecture de `BehaviourIntegratorBase`) :
      l'intégration se fait point par point dans la boucle d'assemblage
      (`performsLocalBehaviourIntegration`, propriétés matériau évaluées
      par point dans un espace de travail) ; sur GPU il faudra intégrer tout
      le lot d'un matériau avant l'assemblage (deux phases).

## Organisation

- `tfel/` : TFEL complet adapté (clone du dépôt officiel de TFEL, branche
  `generic-parallel-behaviours`, base master `fec0ade99`, commits locaux).
  Le dépôt du prototype suit aussi les sources de `tfel/` (copie du dernier
  commit, sans l'historique de TFEL) et la série `git format-patch` depuis
  la base dans `patches/tfel/` (lien `tfel/proto-patches`). Après chaque
  commit dans `tfel/` : `scripts/sync-tfel.sh` met les deux à jour, puis
  commit dans le dépôt du prototype.
- `mgis/` : MGIS (dépôt officiel, base master `818f0f9`, branche
  `gpu-proto`), non modifié pour l'instant.
- `mfem-mgis/` : mfem-mgis (dépôt officiel, base master `4c8a91d815`,
  branche `gpu-proto`), non modifié pour l'instant.
- `build/tfel`, `install/tfel` : compilation et installation de TFEL
  (`build-tfel.sh`, puis `make -j16 install` dans `build/tfel`).
- `scripts/` : `compile-sweep.sh` (balayage des lois de test),
  `host-calls.sh` (appels hôte depuis le code GPU),
  `annotate-function-templates.py`.
- `results/` : résultats des balayages. `scratch/` : sondes jetables.
- `DESIGN.md` : décisions de conception (anglais).
- Outils disponibles sans sudo : nvcc 13 (view spack tfelgpu_local), SDK
  NVIDIA HPC 25.11 (`/opt/nvidia/hpc_sdk/Linux_x86_64/25.11`, `nvc++`,
  CUDA 12.9 et 13.0), `hipcc` 5.7 (fonctionne sur GPU NVIDIA avec le CUDA
  12.9 du SDK : `HIP_PLATFORM=nvidia HIP_COMPILER=nvcc CUDA_PATH=...`),
  Kokkos (view spack, CUDA sm_120). SYCL : AdaptiveCpp à compiler après
  `sudo apt install libclang-20-dev libclang-cpp20-dev libomp-20-dev lld-20
  libpolly-20-dev libzstd-dev libedit-dev libboost-context-dev
  libboost-fiber-dev` (paquets système installés le 05/10).

## Entrées

### 2026-10-05 — état des lieux de TFEL master

- Interface `generic-parallel` des propriétés matériau (#966) :
  backends CUDA, HIP, C++ standard, SYCL ; deux prototypes (avec et sans pas
  de mémoire) ; paramètres dans une structure POD passée par valeur au kernel ;
  erreurs dans un tableau ; métadonnées exportées ; tests dans le répertoire
  d'installation (`test-suite/material-properties/generic-parallel`).
- Lois de comportement : le code généré est déjà presque compatible GPU
  (méthodes `TFEL_HOST_DEVICE`, aucune exception). Branche
  `411-introduce-a-derived-generic-interface-for-gpus` : ébauche de 2023,
  abandonnée.

### 2026-10-05 — sonde J2 (DSL par défaut) dans un kernel CUDA

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

### 2026-10-05 — jalon 1 terminé : code généré compilable pour le GPU

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

### 2026-10-05 — jalon 2 : interface `generic-parallel` des lois

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

### 2026-10-05 — jalon 3 et balayage des lois de test

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

### 2026-10-05 — lois de matériau, valeurs propres et pile sur GPU

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

### 2026-10-05 — réorganisation, commits locaux, piège d'environnement

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

### 2026-10-05 — HIP sur GPU NVIDIA (C1, commit 0007)

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

### 2026-10-05 — `stlpar` sur GPU avec nvc++ 25.11 (C2, en cours)

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

### 2026-10-05 — backends Kokkos et SYCL, factorisation, corrections CMake

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

### 2026-10-05 — matrice complète des backends, corrections

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

### 2026-10-05 (soir) — B5 grandes déformations, B2 initialisations et post-traitements

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

### 2026-10-05 (soir, suite) — B3 rotations des lois orthotropes

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

### 2026-10-05 (nuit) — B4 monocristaux, défauts trouvés par nvcc

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

### 2026-10-05 (nuit, suite) — B1 couverture numérique, nouveaux défauts nvcc et clang

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

### 2026-10-05 (nuit, suite) — B6 modèles autonomes

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

### 2026-10-05 (nuit, suite) — B7 solveur TDLS sur GPU

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

### 2026-10-05 (nuit, suite) — variables de comportement et balayage complet

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

### 2026-10-05 (nuit, suite) — brique StandardElastoViscoPlasticity et en-têtes de TFEL

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

### 2026-10-05 (nuit, suite) — commits 0016 à 0018

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

### 2026-10-05 (soir) — C5 backend OpenMP, points d'entrée factorisés

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

### 2026-10-05 (soir, suite) — propriétés matériau : backends Kokkos et OpenMP

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

### 2026-10-05 (soir, suite) — matrice c, OpenMP, modèles compilés par nvc++

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
