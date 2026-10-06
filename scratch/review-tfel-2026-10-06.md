# Relecture de conception de la branche generic-parallel-behaviours (TFEL)

Date : 2026-10-06. Diff relu : fec0ade99..generic-parallel-behaviours (40 commits, 364 fichiers, +18421/-4875).

## FICHIERS LUS / RESTANTS

LUS :
- DESIGN.md, docs/web/generic-parallel-behaviour-interface.md
- GenericParallelBehaviour/BatchData.h, Integrate.hxx, GPU.hxx, ParallelSTL.hxx, Kokkos.hxx, OpenMP.hxx, SYCL.hxx, FiniteStrainIntegrators.hxx, ModelIntegrator.hxx
- diff de GenericBehaviour/Integrate.hxx
- GenericParallelBehaviourInterface.{hxx,cxx}, GenericParallelModelInterface.cxx
- GenericParallel/Behaviour/* (backends générateurs), DeviceStackSize.hxx
- diff MaterialProperty/* (ParallelSTL, BackendBase, CxxProgrammingModelBackendBase)
- diff AbstractBehaviourInterface, BehaviourCodeGeneratorBase (.hxx, .cxx en grande partie), GenericModelInterface.hxx
- MechanicalBehaviour.hxx, TFELConfig.hxx
- config-*.json, test.check, compare-all.sh, compare-functions.sh, compare-tdls.sh, compare.cxx (structure)

- BehaviourRotationFunctions.hxx, CodeGeneratorUtilities.cxx (diff), MFrontMaterialPropertyInterface.cxx (diff, mfront::device)
- GenericBehaviour/StressMeasures.hxx, plan de GreenLagrangeStrainIntegrate.hxx
- include/TFEL : TFELConfig, ConstantObject, TrivialGenTypeBase, GenTypeBase, BoundsCheck, MechanicalBehaviour, GenericFixedSizeArray, TensorViewFromStensor, StensorComputeEigenVectors (commit 9cd913eb2), FiniteStrainBehaviourTangentOperator, tmatrix (change_basis), LinearHomogenizationSchemes, DDIF2Base (échantillons)
- cmake/modules/tfel.cmake, ImplicitCodeGeneratorBase (diff), PerformanceProfiling, RungeKuttaDSLBase, BehaviourDSLCommon (diffs)
- messages des commits de contournement (9cd913eb2, a22a3cff6, 58f4ced56, f008859df, b0d52da8f, 494929641, 1b621828d)

RESTANTS : aucun (synthèse classée écrite en fin de fichier)

## CONSTATS (au fil de la lecture)


### Lot 1 : en-têtes d'exécution (GenericParallelBehaviour/*)

- [MODÉRÉ, duplication] Chaque backend réécrit les mêmes quatre fonctions
  d'entrée (integrate, executeInitializeFunction, executePostProcessing,
  rotate) qui ne font que lier `makeApply` : GPU.hxx:280-364,
  ParallelSTL.hxx:103-203, OpenMP.hxx (même structure), Kokkos.hxx (idem).
  Environ 100 lignes par backend, documentation comprise. Le code généré
  pourrait appeler directement `mfront::gpb::integrateBatch<...>(s, d, p,
  policy, gpu::makeApply(...))` : chaque backend ne fournirait que son
  `apply`/`makeApply`. La conception annoncée (DESIGN.md, "each backend
  provides a single primitive") est juste, mais la mise en oeuvre ne va pas
  au bout.
- [MODÉRÉ, duplication] SYCL.hxx:369-581 garde ses propres points d'entrée
  (4 fonctions de 40 à 50 lignes) qui répètent chacun : check, test de la
  queue, création de DeviceBatchData, transfert, dispatch, copyOutputs. La
  seule raison est le transfert des champs. On pourrait transférer d'abord,
  puis appeler `integrateBatch` avec le BatchData transféré. Incohérence
  mineure : `std::unique_ptr<DeviceBatchData>` pour trois fonctions,
  `std::optional<DeviceCopies>` pour `rotate` (SYCL.hxx:563).
- [MODÉRÉ, duplication] Integrate.hxx : `dispatch`, `dispatchInitializeFunction`,
  `dispatchPostProcessing` (Integrate.hxx:1081-1092, 1241-1258, 1286-1303)
  répètent le même aiguillage direct/copie. `PointIntegrator`,
  `InitializeFunctionExecutor`, `PostProcessingExecutor` répètent les mêmes
  membres (BatchData, InputFields, parameters, policy). Un aiguillage
  générique (`std::bool_constant<direct>`) suffirait.
- [MINEUR, duplication] Integrate.hxx:110-210 : `gather`/`gatherValues` et
  `scatter`/`scatterValues` font la même chose, une version sur
  `LocalValues<N>`, l'autre sur un pointeur nu.
- [MODÉRÉ, sur-ingénierie / boilerplate] FiniteStrainIntegrators.hxx : trois
  structures (FiniteStrainIntegrator, GreenLagrangeStrainIntegrator,
  LogarithmicStrainIntegrator), plus StandardIntegrator dans Integrate.hxx,
  qui ne diffèrent que par l'espace de noms appelé. 300 lignes dont 80 %
  de documentation répétée. Le style "stratégie = struct de fonctions
  statiques templates" est raisonnable, mais ces structures auraient leur
  place dans `mfront::gb::*` (DESIGN.md parle d'ailleurs de
  `mfront::gb::StandardIntegrator`, qui n'existe pas : le type est dans
  `mfront::gpb`).
- [IMPORTANT, intégration MGIS/mfem-mgis] GPU.hxx:178-245 : chaque appel fait
  cudaGetDevice, cudaDeviceGetAttribute, cudaMallocAsync,
  cudaFuncGetAttributes, cudaDeviceGetLimit (sous mutex), lancement, copie,
  cudaFreeAsync puis `cudaStreamSynchronize`. L'appel est donc toujours
  bloquant et le bilan (OutputStatus) est toujours rapatrié sur l'hôte. Pour
  mfem-mgis avec assemblage sur GPU, cela impose une synchronisation par
  appel et par partie de matériau, sans moyen de l'éviter. Même chose pour
  SYCL (`wait_and_throw` 3 fois, SYCL.hxx:338-343) et Kokkos (réduction
  dans un scalaire). Aucune variante asynchrone n'est prévue dans le
  contrat C.
- [IMPORTANT, intégration] GPU.hxx:49-83 et DeviceStackSize.hxx : la
  bibliothèque générée modifie une limite globale du contexte CUDA
  (`cudaDeviceSetLimit(cudaLimitStackSize)`) à l'insu de l'appelant. Cela
  change la mémoire locale réservée pour tous les noyaux du processus
  (MFEM, hypre). C'est justifié techniquement mais doit être documenté
  comme effet de bord et pouvoir être désactivé par l'appelant (seul le
  backend stlpar/openmp/sycl a une option, CUDA/HIP l'applique toujours).
- [MODÉRÉ, contrat C] BatchData.h : pas de numéro de version ni de taille de
  structure. DESIGN.md:55-57 dit que la structure "can evolve without
  breaking the prototype", mais un MGIS compilé contre une ancienne
  version passera silencieusement une structure de mauvaise taille. Un
  champ `version` (ou un symbole exporté donnant la version du contrat)
  coûte peu et protège MGIS.
- [MINEUR, contrat C] `statuses` est un `int*` alors que tout le reste est un
  champ avec strides (BatchData.h:152). `rdts` est un champ, `statuses` non.
  Conventions de statut inversées par rapport au `generic` pour les
  fonctions d'initialisation et de post-traitement (1 = succès contre 0),
  documentée mais source d'erreurs pour MGIS.
- [MINEUR] Les foncteurs passés aux noyaux copient tout le `BatchData`
  (environ 25 champs) plus `InputFields` (4 x N champs de 24 octets plus
  tailles) plus les paramètres. Pour une loi avec beaucoup de propriétés
  matériau, la taille des arguments de noyau grossit vite (limite 4 Ko
  avant CUDA 12.1). Rien de bloquant, mais à surveiller.
- [BIEN] Integrate.hxx : séparation nette "vue d'un point" (`PointView`),
  "fonction appliquée en un point" (foncteurs trivialement copiables),
  "primitive du backend" (`apply`). Pas de lambda étendue requise. Le
  choix direct/copie est fait une fois sur l'hôte. Les commentaires
  expliquent les choix et les contournements de compilateurs.
- [BIEN] Réduction par warp puis atomiques seulement si nécessaire
  (GPU.hxx:121-156), clé ordonnée pour le min de rdt (Integrate.hxx:765-778).
  Simple et correct.

### Lot 2 : pilotes du `generic` (GenericBehaviour/Integrate.hxx)

- [BIEN] Le découpage "coeur TFEL_HOST_DEVICE + fonction historique qui décode
  K" est propre et évite toute duplication de la physique entre `generic` et
  `generic-parallel`. `ComputationOptions` remplace des décodages ad hoc de
  `K[0]` dispersés. `MethodCaller<m>` garde l'API par pointeur de membre sur
  l'hôte.
- [MODÉRÉ, processus] Le commit 494929641 mélange le refactoring et six
  corrections de comportement du `generic` sur CPU (statut `if (r)` des
  pilotes grandes déformations, DTAU_DDF en Hencky, variables internes de
  fin de pas dans les post-traitements, etc.). Ce sont de vraies
  corrections, bien décrites dans le message, mais elles contredisent le
  principe 4 de DESIGN.md ("No change of the CPU behaviour of existing
  interfaces") et devraient être des commits (voire des PR) séparés pour
  l'amont : elles concernent tous les utilisateurs du `generic`.
- [MINEUR] `callAndReportExceptions(d, f, on_exception)` : le troisième
  argument est `[] {}` dans 3 appels sur 4. Une surcharge à deux
  arguments serait plus lisible.

### Lot 3 : générateurs (interfaces et backends)

- [MODÉRÉ, cohérence d'architecture] Deux philosophies coexistent dans la
  même famille `generic-parallel`. Les propriétés matériau (amont) écrivent
  les boucles et noyaux sous forme de chaînes dans le générateur
  (MaterialProperty/ParallelSTLBackend.cxx:68-120, OpenMPBackend.cxx et
  KokkosBackend.cxx ajoutés par la branche dans ce style). Les lois de
  comportement mettent toute la logique dans des en-têtes templates et le
  générateur n'écrit qu'un appel (Behaviour/*Backend.cxx, 40 à 70 lignes
  chacun). La seconde approche est nettement meilleure, mais la branche a
  étendu la première (Kokkos, OpenMP pour les propriétés matériau) au lieu
  de converger. Deux hiérarchies de backends, deux fabriques
  (`generic_parallel::material_property::BackendFactory` et
  `generic_parallel::behaviour::BackendFactory`), même clé d'options
  `backend` dans les fichiers de configuration.
- [MODÉRÉ, duplication] Analyse des options dupliquée entre les deux
  hiérarchies : `execution_policy` (MaterialProperty/ParallelSTLBackend.cxx:
  29-46 et Behaviour/ParallelSTLBackend.cxx:27-49, mêmes chaînes, messages
  légèrement différents), `device_stack_size` (5 copies), `data_transfer_policy`
  (SYCL). Une petite fonction commune `parseExecutionPolicy` suffirait.
- [MINEUR, duplication] GenericParallelBehaviourInterface.cxx et
  GenericParallelModelInterface.cxx répètent `setOptions` (identiques au nom
  de classe près), `getLibraryName`, `getHeaderFileName`,
  `getSourceFileName`, `getTargetsDescription`, et l'écriture des points
  d'entrée (vérification du pointeur nul, `handleExceptions`, paramètres,
  `writeCall`).
- [MODÉRÉ, conception fragile] BehaviourCodeGeneratorBase.cxx,
  `writeInterfacesCode` : la compatibilité `generic` + `generic-parallel`
  repose sur la déduplication de blocs de code par égalité de texte. Un
  espace ou un commentaire de différence suffit à produire des définitions
  en double (erreur de compilation peu lisible pour l'utilisateur). Une
  relation explicite ("cette interface réutilise les membres de
  `generic`", par exemple une méthode `getBaseInterfaceForGeneratedMembers`)
  serait plus robuste. C'est pragmatique et court, donc acceptable pour un
  prototype, mais c'est le genre de mécanisme implicite que l'amont refusera
  probablement.
- [IMPORTANT, conception] `areDynamicallyAllocatedVectorsAllowed()` et
  `getParametersSourceInBehaviourConstructor()` sont des propriétés d'une
  interface (AbstractBehaviourInterface.hxx:93-111) alors que la classe de
  la loi est commune à toutes les interfaces. Ajouter `generic-parallel`
  change donc la classe générée pour `generic`, `castem`, etc. dans le même
  appel de mfront (pas de vecteurs dynamiques, constructeur supplémentaire).
  C'est voulu, mais cela devrait être documenté comme tel et nommé comme une
  contrainte sur la classe (par exemple "requiresDeviceCompatibleClass").
- [IMPORTANT, sémantique] Paramètres traités comme variables statiques
  (BehaviourCodeGeneratorBase.cxx, `writeBehaviourParameters`) : en code
  device, `static constexpr` devient un membre `const` non statique sous
  `#ifdef TFEL_DEVICE_CODE`. La classe n'a donc pas la même disposition
  mémoire ni la même définition dans les passes hôte et device de nvcc,
  clang et hipcc (violation de l'ODR). Cela marche parce que l'objet loi
  n'est jamais transmis de l'hôte au device, mais c'est fragile : le jour
  où MGIS ou mfem-mgis voudra construire ou copier une loi côté hôte pour
  le device, cela cassera silencieusement. Plus simple et sûr : la même
  déclaration dans les deux passes (membre `const` partout, au prix de
  quelques octets), ou une constante device de portée namespace sur le
  modèle de `getConstantObject`. Lignes : BehaviourCodeGeneratorBase.cxx:3809-3814.
- [IMPORTANT, sémantique et tests] TFELConfig.hxx:383-390 :
  `TFEL_DEVICE_CODE` est défini dans toute l'unité de traduction avec
  `nvc++ -stdpar=gpu`, `nvc++ -mp=gpu` et AdaptiveCpp SSCP. Avec les fichiers
  de configuration fournis, la bibliothèque `generic` compilée dans le même
  appel (`--interface=generic --interface=generic-parallel`) l'est aussi avec
  ces options : elle perd les exceptions, les avertissements de bornes, les
  paramètres lus depuis les initialiseurs dans certains chemins, etc. Dans
  les tests de ces configurations (stlpar-parunseq-nvhpc-gpu/test.check,
  openmp-offload-nvhpc-gpu, sycl-default-acpp), la "référence" n'est donc
  pas indépendante : les chemins device sont comparés aux chemins device,
  avec le même compilateur. Le nom `TFEL_DEVICE_CODE` est aussi trompeur
  dans ce cas ("peut être du code device"). Il faudrait compiler la
  référence avec le compilateur hôte habituel, ou au moins le dire dans les
  tests.
- [MINEUR, bug de contrat] Points d'entrée générés
  (GenericParallelBehaviourInterface.cxx, `write_function`, et
  GenericParallelModelInterface.cxx:318-323, GenericParallelBehaviourInterface.cxx:725-730) : si `mfront_batch_data` est
  nul, `report(-2)` est appelé sans `initialize` : `number_of_failures`,
  `number_of_unreliable_results` et `rdt` du statut restent non
  initialisés. Les fonctions de rotation, elles, appellent `initialize`
  avant.
- [MINEUR, désordre] GenericParallelBehaviourInterface.cxx:102-111 : bloc de
  documentation orphelin (celui de `getPointDataSizes`) suivi d'un autre
  bloc ; la fonction documentée est ligne 149. Préfixes de noms
  générés incohérents : `mgb_d`, `mgpb_parameters`, `mfront_parameters`,
  `mfront_data`.
- [MINEUR] GenericParallelModelInterface.cxx:344-346 : la politique hors
  bornes passée aux modèles est `None` et ignorée (ModelIntegrator.hxx:59,
  paramètre anonyme). Le commentaire le dit, mais l'API transmet une valeur
  inutile.
- [MINEUR, API] KokkosBackend : la chaîne `execution_space` est recopiée
  telle quelle comme argument template dans le code généré
  (Behaviour/KokkosBackend.cxx:39-43), sans validation.
- [BIEN] `GPUBackend` factorise CUDA et HIP (une seule classe, deux noms,
  deux extensions). La fabrique et la validation des options suivent
  exactement le modèle des propriétés matériau. Le générateur de
  l'interface fait 840 lignes lisibles, sans physique.
- [BIEN] Les modèles réutilisent `writeGenericModelClass` avec une option
  `device` plutôt qu'une copie du générateur `generic`.

### Lot 4 : tests

- [MODÉRÉ, maintenabilité] La liste des 29 fichiers .mfront est recopiée
  à l'identique dans les 15 test.check (vérifié : même empreinte md5 pour
  les 15 configurations). Ajouter un cas de test oblige à modifier 16
  fichiers. compare-all.sh centralise les comparaisons mais pas la
  génération.
- [MODÉRÉ, maintenabilité] compare-all.sh passe à la main les tailles de
  chaque loi (`6 6 2 7 1 36 ...`) alors que la bibliothèque `generic`
  exporte ces informations (symboles lus par MGIS). Toute modification d'une
  loi de test oblige à recalculer ces nombres. Pour les rotations, 8
  arguments factices (`0 0 0 0 0 0 0 - 0 --rotation=6`) montrent que
  l'interface en ligne de commande de compare.cxx est à bout de souffle.
- [MINEUR] compare.cxx (1235 lignes) mêle cinq environnements d'exécution
  par macros (Kokkos, SYCL, SYCL mémoire hôte, CUDA/HIP, CPU). C'est
  raisonnable pour un banc de test, mais un futur test via MGIS
  remplacerait avantageusement une grande partie de ce code.
- [BIEN] La référence est l'interface `generic` point par point, sur deux
  pas de temps, trois dispositions mémoire, trois types de calcul, avec
  échecs partiels, énergies, vitesse du son, fonctions d'initialisation,
  post-traitements et rotations. La couverture fonctionnelle est très
  bonne.

### Lot 5 : en-têtes TFEL (include/TFEL)

Mesure : 163 fichiers, +2362/-1836 lignes ; 757 lignes ajoutées sur 2362
contiennent une annotation TFEL_HOST_DEVICE. Le reste : reformatage induit,
remplacement d'exceptions, et une poignée de vrais changements.

- [BIEN] Les contournements de compilateurs dans les en-têtes TFEL sont peu
  nombreux et tous commentés avec version et symptôme :
  ConstantObject.hxx:29,44 (nvcc, nvc++), BoundsCheck.hxx:49 (nvcc),
  TensorViewFromStensor.hxx:57 (nvc++ 25.11, switch), GenericFixedSizeArray.ixx:36
  (nvc++ 25.11, initialisation par valeur). Les messages de commit
  (58f4ced56, a22a3cff6, f008859df) donnent les reproducteurs. Pas
  d'accumulation sauvage.
- [BIEN] BoundsCheck.hxx : passage de `std::string` à `const char*` et des
  spécialisations par dimension à une boucle : -365/+199 lignes, plus simple
  qu'avant. Retour booléen sur device, exception sur l'hôte : compromis
  clair et documenté.
- [BIEN] TrivialGenTypeBase (union étiquetée, 249 lignes, test unitaire
  TrivialGenTypeBaseTest.cxx) : justifié (std::variant non supporté par nvc++
  sur device), isolé, même interface que GenTypeBase.
- [MODÉRÉ, invasif] Ruptures d'API/ABI dans le coeur de TFEL, toutes
  justifiées mais à négocier avec l'amont :
  * MechanicalBehaviour.hxx : suppression des fonctions virtuelles (commit
    1b621828d), remplacées par un concept. Tout code tiers qui écrit
    `override` sur ces méthodes ne compile plus ;
  * BoundsCheck : signature `const char*` au lieu de `std::string` ;
  * src/Material/LogarithmicStrainHandler.cxx supprimé : deux symboles
    exportés disparaissent de libTFELMaterial ;
  * GenericFixedSizeArray : constructeurs de copie et opérateurs
    d'affectation supprimés (pour être trivialement copiable).
- [MODÉRÉ, changement de résultats CPU] StensorComputeEigenVectors.hxx
  (commit 9cd913eb2) : pour une valeur propre multiple en 3D, les vecteurs
  propres sont désormais calculés par Jacobi (`fses::syevj3`) au lieu de
  l'algorithme analytique, sur l'hôte aussi. C'est un changement
  numérique (cas dégénéré) pour tous les utilisateurs de TFEL, motivé par
  le device. Il aurait pu être limité au code device, ou au moins être
  signalé comme changement de résultats.
- [MINEUR, churn] GenTypeBase.hxx : toutes les méthodes reçoivent
  TFEL_HOST_DEVICE et `std::get` devient `*std::get_if`, alors que
  GenTypeBase (std::variant) n'est plus utilisé en code device depuis
  l'introduction de TrivialGenTypeBase dans le même commit (329a740e9).
  Modification non nécessaire d'un en-tête central.
- [MINEUR] ConstantObject.hxx:57-70 : le type de retour de
  `getConstantObject` dépend du compilateur (`T` avec nvc++, `const T&`
  ailleurs). La règle "lier à une référence constante" est documentée, mais
  c'est un piège d'API (un `auto x = ...` copie partout, un `auto& x = ...`
  ne compile qu'ailleurs que chez nvc++).
- [MINEUR, cosmétique] Une partie du diff est du reformatage clang-format
  (clauses `requires`, ordre `[[nodiscard]] TFEL_HOST_DEVICE`). L'ordre des
  attributs est une vraie correction (grammaire C++), mais la
  réindentation des `requires` (LinearHomogenizationSchemes.ixx, etc.)
  alourdit la relecture amont.

### Lot 6 : points qui gêneront MGIS et mfem-mgis

- [IMPORTANT] Appels toujours synchrones et bilan toujours rapatrié sur
  l'hôte (voir lot 1). mfem-mgis sur GPU voudra enchaîner intégration,
  assemblage et solveur sur le même flux sans synchronisation.
- [IMPORTANT, sémantique silencieuse] Sur device, plusieurs valeurs
  diffèrent silencieusement de l'hôte : paramètres des lois matériau
  appelées (`mfront::device`, valeurs par défaut,
  MFrontMaterialPropertyInterface.cxx, `writeDeviceMaterialPropertyBody`),
  paramètres des variables de comportement (valeurs par défaut), bornes
  standard des lois matériau non vérifiées, avertissements de bornes
  supprimés. Documenté dans generic-parallel-behaviour-interface.md:453-474,
  mais aucune détection à l'exécution : un utilisateur de MGIS qui modifie
  un paramètre d'une loi matériau aura des résultats GPU différents des
  résultats CPU sans aucun message.
- [MODÉRÉ] Deux bibliothèques (`libBehaviour.so` et
  `libGenericParallelBehaviour-<backend>.so`) exportent les mêmes noms C
  (fonction d'intégration, `setParameter`, `setOutOfBoundsPolicy`) et
  contiennent chacune leur copie de l'objet compilé de la loi (Makefile
  généré : `X.o` lié dans les deux). L'état (paramètres, politique hors
  bornes) est donc dupliqué : si MGIS se replie sur `generic` pour une
  partie des points, il devra répercuter chaque `setParameter` dans les deux
  bibliothèques. Avec un chargement RTLD_GLOBAL (ou une édition de liens
  directe), les noms identiques peuvent même s'interposer.
- [MODÉRÉ] Pas de version du contrat C (voir lot 1).
- [MINEUR] Les messages d'erreur par point sont perdus
  (`bd.error_message = nullptr`, Integrate.hxx:531, `PointView`) : en cas
  d'échec sur CPU (stlpar, openmp hôte), le message de l'exception de la
  loi n'est jamais remonté, seulement "integration failed at some
  integration points". Garder le premier message (au moins sur l'hôte)
  aiderait beaucoup au débogage dans MGIS.
- [MINEUR] GPU.hxx:54 : le mutex de `ensureStackSize` est statique dans un
  modèle de fonction, donc un mutex par noyau. La séquence lire-comparer-
  écrire de la limite de pile n'est pas protégée entre noyaux différents ni
  entre bibliothèques : avec des appels concurrents depuis plusieurs fils
  hôte (matériaux différents dans mfem-mgis), une limite plus faible peut
  écraser une limite plus forte juste avant un lancement.
- [BIEN] Le contrat C par champs (pointeur + 2 strides), avec propriétés
  matériau et variables externes "par variable", colle exactement à
  `MaterialDataManager` de MGIS (valeurs uniformes ou par point) : MGIS
  pourra passer ses données sans copie. Le symbole
  `<behaviour>_generic_parallel_backend` et la réutilisation des symboles
  du `generic` simplifient le chargement. Le rôle de chaque bibliothèque
  (DESIGN.md:234-265) est clair et sain : TFEL n'alloue pas, MGIS ne
  dépend d'aucun outil GPU.

### Divers vérifiés (cosmétique)

- cmake/modules/tfel.cmake:554 : `_tfel-generic-parallel_behaviours_only_configurations`
  est une liste vide, utilisée dans un test : échafaudage mort.
- RungeKuttaDSLBase.cxx:753-760 et RungeKuttaCodeGeneratorBase.cxx:2254-2265
  écrivent encore leurs `#ifdef TFEL_DEVICE_CODE` et leurs branches
  `useQt` à la main au lieu d'utiliser `makeHostAndDeviceCode` et
  `getInvalidStandardTangentOperatorFlagCondition` introduits par la
  branche : harmonisation incomplète.
- DESIGN.md:94-97 cite `mfront::gb::StandardIntegrator` et
  `finite_strain::Integrator`, qui n'existent pas sous ces noms
  (`mfront::gpb::StandardIntegrator`, `mfront::gpb::FiniteStrainIntegrator`).

### Compléments (vérifiés)

- [MODÉRÉ, cohérence] Conventions de statut différentes dans la même
  famille `generic-parallel` : propriétés matériau (`mfront_gmp_OutputStatus`,
  GenericMaterialProperty/OutputStatus.h:26-75) : 0 succès, 1 prudence,
  -2 exception, -4 résultat invalide, -6 pile (ajouté par la branche) ;
  lois (`mfront_gpb_OutputStatus`, BatchData.h:187-209) : 1 succès,
  0 peu fiable, -1 échec, -2 arguments invalides, -3 erreur d'exécution
  ou exception. Le choix "comme `generic`" se défend pour les lois, mais
  -2 n'a pas le même sens d'une structure à l'autre et MGIS manipulera les
  deux.
- [MODÉRÉ, périmètre] Six backends et quinze configurations testées. La
  plupart des contournements et des chemins particuliers viennent de
  nvc++ (21 mentions dans les lignes ajoutées hors tests : switch, init
  par valeur, référence à un drapeau local, pointeurs de membres,
  assertions, mémoire `managed` imposée, pile via dlopen) et d'AdaptiveCpp
  (pas de réductions, pas de récursion, pas de fonctions virtuelles).
  Pour l'amont et pour MGIS, un premier lot CUDA/HIP + stlpar/OpenMP hôte
  serait beaucoup plus facile à faire accepter et à maintenir ; les autres
  pourraient suivre.
- [MODÉRÉ, coût de compilation] Chaque point d'entrée est instancié deux
  fois (`direct` vrai ou faux, Integrate.hxx:1081-1092), donc deux noyaux
  complets par fonction, initialisations et post-traitements compris. Le
  gain du mode direct n'est mesuré nulle part (JOURNAL.md ne donne pas de
  comparaison). Pour des lois lourdes (monocristaux), cela double un temps
  de compilation déjà long. À mesurer, ou à rendre optionnel.

## SYNTHÈSE CLASSÉE

Verdict général : la conception de fond est bonne (une seule
implémentation par point partagée avec `generic`, logique dans des
en-têtes, contrat C simple et adapté à MGIS). Les défauts sérieux sont des
effets de bord sur l'existant et des divergences silencieuses hôte/device,
pas la structure du nouveau code. La duplication est réelle mais locale et
facile à résorber.

### A. Problèmes sérieux

A1. `TFEL_DEVICE_CODE` défini dans toute l'unité de traduction pour nvc++
-stdpar/-mp et AdaptiveCpp SSCP (TFELConfig.hxx:383-390). Avec les
fichiers de configuration fournis, la bibliothèque `generic` du même appel
est compilée avec ces options (Makefile généré : `X.o` lié dans
libBehaviour.so et dans libGenericParallelBehaviour-stlpar.so). Elle perd
alors exceptions et avertissements, et le constructeur par défaut des
variables de comportement prend les valeurs par défaut des paramètres même
sur l'hôte (BehaviourCodeGeneratorBase.cxx:2446-2455). Les tests de ces
configurations (stlpar-parunseq-nvhpc-gpu, openmp-offload-nvhpc-gpu,
sycl-default-acpp) comparent donc chemin device contre chemin device, avec
le même compilateur : la référence n'est pas indépendante.

A2. Divergences silencieuses hôte/device, documentées mais jamais
détectées : paramètres des lois matériau appelées (`mfront::device`,
valeurs par défaut, MFrontMaterialPropertyInterface.cxx,
`writeDeviceMaterialPropertyBody`), paramètres des variables de
comportement, bornes standard non vérifiées, avertissements supprimés. Un
`setParameter` côté MGIS donnera des résultats GPU différents des
résultats CPU sans aucun message. Il faudrait au minimum un refus (ou un
avertissement) à l'exécution quand un paramètre modifié n'est pas
transmis au device.

A3. Appels toujours bloquants, effets de bord globaux. GPU.hxx:178-245 :
cudaMallocAsync, copie, lancement, copie, cudaFreeAsync puis
cudaStreamSynchronize à chaque appel ; SYCL.hxx:338-343 : trois
`wait_and_throw` ; Kokkos : réduction dans un scalaire. Le bilan est
toujours rapatrié sur l'hôte. La limite de pile du contexte CUDA est
modifiée sans possibilité de refus pour CUDA/HIP (GPU.hxx:49-83), et via
dlopen de libcudart pour les autres (DeviceStackSize.hxx). Pour mfem-mgis
avec assemblage sur GPU, c'est la première chose à revoir : une variante
asynchrone (bilan laissé en mémoire device, synchronisation à la charge
de l'appelant) et un réglage de pile explicite.

A4. La classe de loi n'a pas la même définition dans les passes hôte et
device : paramètres traités comme variables statiques déclarés `static
constexpr` sur l'hôte et membres `const` non statiques sur le device
(BehaviourCodeGeneratorBase.cxx:3809-3814). Violation de l'ODR pour nvcc,
clang et hipcc. Cela marche tant que l'objet loi n'est jamais construit
sur l'hôte puis copié vers le device ; ce sera faux le jour où MGIS ou
mfem-mgis le feront.

A5. Changements de comportement de l'existant mêlés au prototype,
contraires au principe 4 de DESIGN.md ("No change of the CPU behaviour of
existing interfaces") : six corrections des pilotes `generic` dans le
commit 494929641, corrections des rotations (79eacfee1), vecteurs propres
d'une valeur propre multiple calculés par Jacobi sur l'hôte aussi
(9cd913eb2), ruptures d'API (MechanicalBehaviour sans fonctions
virtuelles, 1b621828d ; BoundsCheck en `const char*`), rupture d'ABI
(LogarithmicStrainHandler.cxx supprimé de libTFELMaterial). Chaque
changement est justifié et bien expliqué dans son message de commit, mais
l'ensemble doit être séparé en PR distinctes pour l'amont et le principe
de DESIGN.md corrigé.

### B. Problèmes modérés

B1. Deux architectures dans la même famille `generic-parallel` :
propriétés matériau générées sous forme de chaînes (la branche y ajoute
Kokkos et OpenMP dans ce style), lois par en-têtes templates. Deux
fabriques, options analysées deux fois (`execution_policy` :
MaterialProperty/ParallelSTLBackend.cxx:29-46 et
Behaviour/ParallelSTLBackend.cxx:22-49 ; `device_stack_size` : six
copies), deux conventions de statut (`mfront_gmp_OutputStatus` : 0 succès,
-2 exception ; `mfront_gpb_OutputStatus` : 1 succès, -2 arguments
invalides).

B2. Duplication dans les en-têtes d'exécution : quatre points d'entrée
recopiés par backend alors qu'ils ne font que lier `makeApply`
(GPU.hxx:280-364, ParallelSTL.hxx:103-203, OpenMP.hxx, Kokkos.hxx) ;
quatre points d'entrée SYCL de 40 à 50 lignes répétant la même séquence
(SYCL.hxx:369-581) ; trois `dispatch*` identiques (Integrate.hxx:1081-1092,
1241-1258, 1286-1303) ; quatre structures d'intégration qui ne diffèrent
que par l'espace de noms appelé (FiniteStrainIntegrators.hxx, 300 lignes).
Le code généré pourrait appeler `mfront::gpb::integrateBatch` avec le
`makeApply` du backend : 300 à 400 lignes de moins.

B3. Mécanismes implicites dans le générateur : compatibilité
`generic` + `generic-parallel` obtenue en dédupliquant des blocs de code
par égalité de texte (`writeInterfacesCode`, BehaviourCodeGeneratorBase.cxx) ;
propriétés d'interface (`areDynamicallyAllocatedVectorsAllowed`,
AbstractBehaviourInterface.hxx:93-111) qui modifient la classe partagée par
toutes les interfaces du même appel. Pragmatique, mais fragile et peu
lisible.

B4. Contrat C et intégration MGIS : pas de version du contrat (alors que
DESIGN.md:53-57 compte sur l'évolution de la structure) ; état dupliqué
entre libBehaviour.so et libGenericParallelBehaviour-*.so (paramètres,
politique hors bornes, mêmes noms C), à répercuter par MGIS en cas de
repli sur `generic` ; message d'erreur par point perdu
(Integrate.hxx:531).

B5. Tests difficiles à maintenir : la liste des 29 fichiers .mfront
recopiée à l'identique dans 15 test.check ; tailles de chaque loi passées à
la main dans compare-all.sh alors que la bibliothèque `generic` les
exporte ; arguments positionnels factices pour les rotations
(`0 0 0 0 0 0 0 - 0 --rotation=6`).

B6. Périmètre : six backends et quinze configurations. La plupart des
contournements viennent de nvc++ (21 mentions dans les lignes ajoutées) et
d'AdaptiveCpp. Un premier lot CUDA/HIP + stlpar/OpenMP hôte serait plus
facile à faire accepter.

B7. Double instanciation direct/copie de chaque point d'entrée, sans mesure
du gain : temps de compilation doublé pour les lois lourdes.

### C. Mineur et cosmétique

- Integrate.hxx:110-213 : `gather`/`gatherValues`, `scatter`/`scatterValues`
  en double.
- Pointeur de batch nul : `report(-2)` sans `initialize`, champs du statut
  non initialisés (GenericParallelBehaviourInterface.cxx:725-730,
  GenericParallelModelInterface.cxx:318-323).
- Bloc de documentation orphelin (GenericParallelBehaviourInterface.cxx:102-111) ;
  préfixes générés incohérents (`mgb_d`, `mgpb_parameters`, `mfront_parameters`).
- GenTypeBase.hxx annoté TFEL_HOST_DEVICE et modifié (`*std::get_if`) alors
  qu'il n'est plus utilisé sur device.
- `getConstantObject` : type de retour qui dépend du compilateur
  (ConstantObject.hxx:57-69).
- GPU.hxx:54 : mutex statique par instanciation de modèle, la séquence
  lire-comparer-écrire de la pile n'est pas protégée entre noyaux.
- Modèles : politique hors bornes `None` passée puis ignorée
  (GenericParallelModelInterface.cxx:344-346, ModelIntegrator.hxx:59).
- KokkosBackend : `execution_space` recopié sans validation dans le code
  généré.
- tfel.cmake:554 : liste vide `_tfel-generic-parallel_behaviours_only_configurations`.
- RungeKutta : `#ifdef TFEL_DEVICE_CODE` et branches `useQt` écrits à la
  main au lieu des nouveaux utilitaires (RungeKuttaDSLBase.cxx:753-760,
  RungeKuttaCodeGeneratorBase.cxx:2254-2265).
- DESIGN.md:94-97 : noms de types qui n'existent pas
  (`mfront::gb::StandardIntegrator`).
- Taille des arguments de noyau (BatchData + InputFields + paramètres
  copiés par valeur) à surveiller pour les lois à nombreuses propriétés.
- Reformatage clang-format de clauses `requires` qui alourdit le diff.

### D. Ce qui est bien fait

- Une seule implémentation par point : les pilotes `generic` sont coupés en
  un coeur TFEL_HOST_DEVICE (options explicites, arguments du constructeur
  transmis) et les fonctions historiques. Aucune physique dupliquée.
- Séparation nette vue d'un point / foncteur trivialement copiable /
  primitive `apply` du backend ; pas de lambdas étendues ; choix
  direct/copie fait une fois sur l'hôte.
- Contrat C par champs (pointeur + deux pas) avec propriétés matériau et
  variables externes "par variable" : correspond exactement à
  `MaterialDataManager`, MGIS pourra passer ses données sans copie.
  Rôles de TFEL, MGIS et mfem-mgis clairs (DESIGN.md:234-265).
- Réduction par warp et atomiques seulement si nécessaire ; clé ordonnée
  pour le minimum de rdt.
- Contournements de compilateurs peu nombreux dans les en-têtes TFEL
  (quatre), tous commentés avec version et symptôme ; reproducteurs dans
  results/nvhpc-bugs.
- BoundsCheck simplifié (-365/+199 lignes) ; TrivialGenTypeBase justifié et
  testé ; `change_basis` ajouté avec test ; rotations partagées entre
  `generic` et `generic-parallel` (BehaviourRotationFunctions) avec
  corrections du `generic`.
- Modèles réalisés en réutilisant `writeGenericModelClass` (option
  `device`) ; backends CUDA et HIP factorisés dans `GPUBackend`.
- Messages de commit exemplaires (motif, défaut de compilateur,
  corrections) ; documentation utilisateur complète ; couverture
  fonctionnelle des tests très bonne (deux pas, trois dispositions, trois
  types de calcul, échecs partiels, énergies, vitesse du son, fonctions
  d'initialisation, post-traitements, rotations, modèles).

### E. Ordre de travail suggéré avant MGIS

1. Variante asynchrone et réglage de pile explicite (A3).
2. Détection des paramètres non transmis au device (A2).
3. Même définition de classe hôte/device (A4).
4. Référence `generic` compilée par le compilateur hôte dans les tests
   (A1).
5. Factoriser les points d'entrée et les `dispatch` (B2), version du
   contrat C (B4).
6. Découper la branche en PR amont : corrections du `generic`, TFEL
   device, interface (A5).
