# D1, MGIS : intégration d'un lot avec une bibliothèque `generic-parallel` (proposition de conception)

Statut : proposition, aucune modification de code. Bases : MGIS master `818f0f9` (`mgis/`), TFEL proto branche `generic-parallel-behaviours` (`tfel/`). Rédigée le 2026-10-05, mise à jour le 2026-10-06 avec l'état de TFEL (section 0).

## 0. Ce que TFEL fournit (mise à jour du 2026-10-06)

- Propriétés matériau et variables externes **variable par variable**
  (`material_properties_by_variable`, `external_state_variables_by_variable`) :
  chaque `FieldHolder` de MGIS devient directement un champ (uniforme : pas
  nul ; en place : `{p+b*s, s, 1}`). Le planificateur du 2.2 et l'espace de
  travail ne sont plus nécessaires.
- Symboles : `<b>_mfront_interface` = "GenericParallel",
  `<b>_generic_parallel_backend` (nom du backend),
  `<b>_generic_parallel_api_version` (= `MFRONT_GPB_API_VERSION`). Le
  symbole d'accessibilité de la mémoire hôte proposé au 1.6 n'est pas
  nécessaire : le backend suffit (table du 3.1).
- Options d'exécution (`mfront_gpb_ExecutionOptions`, membre `execution`
  du lot et des rotations) : file, statistiques asynchrones
  (`mfront_gpb_BatchStatistics`, fonctions C d'initialisation et de
  décodage), `preserve_stack_size`. `BatchExecutionPolicy` (4.1) les
  reprend.
- Messages d'erreur par point non transmis (décision documentée) : la
  bibliothèque `generic` compagnon les donne (question 5).
- Données partagées (`@SharedData`) gérées par la bibliothèque générée :
  rien à faire dans MGIS.
- Chargement depuis un programme compilé par g++ : les bibliothèques nvc++
  exigent que `libm` soit chargée (toujours le cas d'un programme C++
  utilisant `<cmath>`, défaut d'empaquetage de NVHPC 25.11).

### 0.1 Constats des essais d'exécution de toutes les lois (R3, 2026-10-07)

- Les tailles d'un point tiennent compte des tableaux de gradients et de
  forces thermodynamiques : MGIS doit lire les tailles dans la description
  de la loi (`getArraySize` des variables), comme pour l'interface
  `generic`.
- Les paramètres changés par `<b>_setParameter` (ou
  `<b>_<hypothèse>_setParameter` quand la loi ne supporte pas toutes les
  hypothèses) sont lus sur l'hôte à chaque appel : MGIS les règle comme
  aujourd'hui, dans la bibliothèque `generic-parallel` elle-même (les
  bibliothèques `generic` et `generic-parallel` d'une même loi ont chacune
  leurs paramètres, modèles compris depuis R3).
- Les bornes physiques des arguments des propriétés matériau sont
  exportées même sans bornes standard.
- Une loi peut demander une pile de plus de 100 Ko par fil sur GPU
  (polycristaux) : la réserve vaut alors presque toute la mémoire du GPU,
  quel que soit le nombre de points. Un autre processus sur le même GPU
  suffit à la faire échouer (statut -3, « out of memory »). Plusieurs
  rangs MPI par GPU fonctionnent avec les lois usuelles ; MGIS devrait
  exposer la pile demandée (`localSizeBytes` des noyaux, à exporter par
  TFEL) pour que l'appelant puisse le prévoir.
- Les résultats sur GPU diffèrent de ceux du CPU par les arrondis
  (contraction FMA de nvcc). Aux limites de convergence des solveurs
  locaux, un point peut échouer d'un côté et pas de l'autre : les tests de
  MGIS doivent comparer avec une tolérance et sur des chargements loin de
  ces limites.

## Résumé

1. **Chargement.** MGIS détecte la bibliothèque par le symbole `<b>_mfront_interface` == "GenericParallel" (`LibrariesManager::getInterface`, `src/LibrariesManager.cxx:655-662`, jamais appelé aujourd'hui).
   - C'est indispensable : aujourd'hui MGIS charge sans erreur une bibliothèque `generic-parallel` isotrope. Les points d'entrée portent les mêmes noms mais pas la même signature, donc l'intégration plante (comportement indéfini). Une bibliothèque orthotrope échoue au chargement (pas de `rotateArrayOf*`).
   - Pour une telle bibliothèque, `Behaviour` reçoit des pointeurs de lot (`batch_b`, `batch_f`, trois rotations) et `b` reste nul.
2. **Données.** Mappage direct, sans copie, de tous les blocs, propriétés matériau et variables externes comprises (champs par variable de TFEL, section 0).
3. **Mémoire.** Phase 1 : `std::pmr::memory_resource*` dans les initialiseurs, stockage privé en `std::pmr::vector`, aucune rupture d'API, mémoire accessible depuis l'hôte (managée, USM partagée, Kokkos `SharedSpace`). Phase 2 : mémoire uniquement device, avec copie et remplissage fournis par la ressource.
4. **API.** `integrate(const BatchExecutionPolicy&, MaterialDataManager&, const BatchIntegrationOptions&, real dt)`, variantes avec intervalle et `Context`. Statuts -2/-3 levés comme erreurs (pas des échecs de convergence). rdt initial 10 comme le chemin par point. Basculement automatique des surcharges existantes.

## 1. Chargement

### 1.1 Constat

- Une bibliothèque `generic-parallel` exporte les symboles du `generic` (`_mfront_mkt`, `_api_version`, noms et types des variables, paramètres, `_setParameter`, `_setOutOfBoundsPolicy`). `<b>_mfront_interface` y vaut "GenericParallel" (`GenericParallelBehaviourInterface.cxx:356-358`, écrit par `SymbolsGenerator.cxx:745`) contre "Generic" pour le `generic` (`GenericBehaviourInterface.cxx:119-121`).
- Points d'entrée de **mêmes noms** mais de **signatures différentes** :
  - `<b>_<H>` (`GenericParallelBehaviourInterface.cxx:695-704`) ;
  - `<b>_<H>_InitializeFunction_<f>` (706-724) ;
  - `<b>_<H>_PostProcessing_<p>` (725-741) ;
  - `<b>_<H>_rotateGradients`, `_rotateThermodynamicForces[_<mesure>]`, `_rotateTangentOperatorBlocks[_<opérateur>]` (751-781) ; **pas** de `rotateArrayOf*`.
- Nom de la bibliothèque : "GenericParallel" + base + "-" + backend (431-443), par exemple `libGenericParallelBehaviour-cuda.so`. Aucun symbole ne donne le backend.
- **Risque actuel** :
  - `load` (`src/Behaviour.cxx:68-95`) appelle `lm.getBehaviour` (`LibrariesManager.cxx:271-282`) et range la fonction de lot dans `Behaviour::b` (`BehaviourFctPtr`) ; `integrate(BehaviourDataView&, const Behaviour&)` (`Integrate.ixx:22-25`) l'appelle avec un mauvais argument : comportement indéfini.
  - Comportement orthotrope : échec de `load` sur `getRotateArrayOf*` (`Behaviour.cxx:108-121`, `LibrariesManager.cxx:352-369`).
  - `loadBehaviourDescription` (`BehaviourDescription.cxx:124-364`) ne lit jamais l'interface ; la vérification d'`_api_version` (151-158) passe.

### 1.2 Détection

- `BehaviourDescription` (`BehaviourDescription.hxx:35-59`) : champ `std::string interface;` rempli par `lm.getInterface(l, b)` (`LibrariesManager.hxx:93`) ; absent = "" = `generic`.
- `std::string backend;` lu dans un nouveau symbole TFEL `<b>_mfront_gpb_backend` (1.6), sinon déduit du suffixe du nom de bibliothèque (informatif).
- Prédicat public : `bool isGenericParallelBehaviour(const BehaviourDescription&) noexcept;`

### 1.3 Ajouts à `Behaviour` (rétrocompatibles au niveau source, ABI cassée)

Nouvel en-tête C `include/MGIS/Behaviour/BatchBehaviourFctPtr.hxx` (structure de `BehaviourFctPtr.hxx:26-81`), avec une **copie des structures C** de `BatchData.h` (comme `mgis_bv_StateView` duplique `mfront_gb_State`), pour éviter une dépendance de compilation de MGIS à TFEL : `mgis_bv_ConstField`, `mgis_bv_Field`, `mgis_bv_BatchInitialStateView`, `mgis_bv_BatchStateView`, `mgis_bv_BatchDataView`, `mgis_bv_BatchRotationDataView`, `mgis_bv_BatchOutputStatus`. Types compatibles (`mgis_size_type` et `mfront_gb_size_type` = `size_t` ; `mgis_real` = `double` sauf `MGIS_REAL_TYPE`, à vérifier par `static_assert`).

```c
typedef int (*mgis_bv_BatchBehaviourFctPtr)(mgis_bv_BatchOutputStatus* const,
                                            const mgis_bv_BatchDataView* const);
typedef int (*mgis_bv_BatchBehaviourInitializeFctPtr)(
    mgis_bv_BatchOutputStatus* const, const mgis_bv_BatchDataView* const,
    const mgis_bv_ConstField* const /* inputs, may be null */);
typedef int (*mgis_bv_BatchBehaviourPostProcessingFctPtr)(
    mgis_bv_BatchOutputStatus* const, const mgis_bv_Field* const,
    const mgis_bv_BatchDataView* const);
typedef int (*mgis_bv_BatchRotationFctPtr)(
    mgis_bv_BatchOutputStatus* const, const mgis_bv_BatchRotationDataView* const);
```

`Behaviour.hxx` :
- `BehaviourInitializeFunction` (41-58) : `+ BatchBehaviourInitializeFctPtr batch_f = nullptr;`
- `BehaviourPostProcessing` (61-78) : `+ BatchBehaviourPostProcessingFctPtr batch_f = nullptr;`
- `Behaviour` (81-161) : `+ BatchBehaviourFctPtr batch_b = nullptr;`, `+ BatchRotationFctPtr batch_rotate_gradients_ptr, batch_rotate_thermodynamic_forces_ptr, batch_rotate_tangent_operator_blocks_ptr = nullptr;`

Les métadonnées (`inputs`, `outputs`, `Behaviour.cxx:84,91`) restent chargées comme aujourd'hui. Pour une bibliothèque `GenericParallel`, `b`, `f` et `rotate_*_ptr` restent nuls.

### 1.4 `LibrariesManager` (`LibrariesManager.hxx` vers 152-253)

Nouvelles fonctions : `getBatchBehaviour(l, b, h)`, `getBatchBehaviourInitializeFunction(l, b, f, h)`, `getBatchBehaviourPostProcessing(l, b, f, h)`, `getBatchRotationFunction(l, b, h, std::string_view name)`, `getGenericParallelBackend(l, b)`. Chacune vérifie `getInterface(l, b) == "GenericParallel"`. Symétriquement, `getBehaviour`, `getBehaviourInitializeFunction`, `getBehaviourPostProcessing`, `getRotate*` (`LibrariesManager.cxx:227-575`) lèvent une erreur pour une interface "GenericParallel". Logique de suffixe (mesure de contrainte, opérateur tangent, 413-424 et 513-526) factorisée.

### 1.5 Flux de chargement (`src/Behaviour.cxx`)

- `load_behaviour` (68-95) : si `isGenericParallelBehaviour(d)`, remplir `batch_b` et `batch_f` ; sinon code actuel.
- Rotations (108-121 et 171-188) : branche lot (`rotateGradients`, `rotateThermodynamicForces[_<CauchyStress|PK2Stress|PK1Stress>]`, `rotateTangentOperatorBlocks[_<dsig_dF|dPK2_dEGL|dPK1_dF|dtau_ddF>]`).
- `Model` (= `Behaviour`) : `libGenericParallelModel-<backend>.so` couvertes sans travail supplémentaire (`Model.cxx:39-45`).
- Gardes du chemin par point : `Integrate.ixx:22-25` (`if (b.b == nullptr) raise("...use the batch API")`), `executeInitializeFunction`/`executePostProcessing` par point (`Integrate.cxx:555-582, 1006-1018`), rotations par point (`Behaviour.cxx:387-394, 492-499, 607-614`).
- `setParameter` (`LibrariesManager.cxx:931-1002`) fonctionne tel quel ; chaque bibliothèque a sa copie des paramètres (à fixer dans les deux dans les tests).

### 1.6 Symboles TFEL

Faits (section 0) : `_generic_parallel_backend` et
`_generic_parallel_api_version`.

## 2. Correspondance des données

### 2.1 Table (zc = zéro copie)

| Champ `BatchData` | Source MGIS | Mappage |
|---|---|---|
| `number_of_points` | `m.n` (`MaterialDataManager.hxx:238`) ou `e-b` |  |
| `dt` | argument |  |
| `options[0]` | `encodeBehaviourIntegrationOptions(opts)` (`Integrate.cxx:40-46`, +100 si vitesse du son) |  |
| `options[1..2]` | `b.options` (`Behaviour.hxx:160`) | 0 si vide |
| `rdt` | `opts.rdt` (défaut 10, 4.4) |  |
| `tangent_operators` | `m.K`, `K_stride` | zc `{K+b*Ks, Ks, 1}` ; nul si pas d'opérateur |
| `rdts` | nouveau `m.time_step_scaling_factors` (optionnel) | zc |
| `speed_of_sound` | `m.speed_of_sound` | zc si demandée |
| `statuses` | nouveau `m.integration_statuses` | zc |
| `s0/s1.gradients`, `thermodynamic_forces`, `internal_state_variables` | `MaterialStateManager.hxx:176-216` | zc |
| `stored_energy`, `dissipated_energy` | 186-190 si `b.computes*` | zc sinon nul |
| `mass_density` | `optional<MutableFieldHolder>` (208) | uniforme : espace de travail (`ps=0`) ; `span` : zc ; `vector` : zc si mémoire hôte acceptée sinon copie |
| `material_properties_by_variable` | map par variable (199), scalaires | un champ par variable (2.2) |
| `external_state_variables_by_variable` | map par variable (226), uniformes ou champs | un champ par variable (2.2) |
| `queue` | `policy.queue` |  |

### 2.2 Propriétés matériau et variables externes

Champs par variable de TFEL (section 0) : pas de planificateur ni d'espace
de travail. La masse volumique uniforme reste un champ de pas nul.

### 2.3 Options A et B (`JOURNAL.md`, D1)

| | B : MGIS assemble un champ par groupe | A : champ par variable dans `BatchData` |
|---|---|---|
| Changement TFEL | aucun | ABI de `BatchData.h`, `PointView`, `check()`, `areAllFieldsStoredPointByPoint`, `PointTraits`, copies SYCL, modèles, `compare.cxx`, matrice complète |
| Changement MGIS | planificateur et espace de travail (~250 lignes) | descripteurs (~120 lignes) |
| Mémoire en plus | ≤ 2·n·(nmp+nesv) réels (cas 5) | ≈ 0 |
| Temps par appel | cas 5 : copie O(n·nvals) et transfert sur GPU | rassemblement dans le noyau |
| Champs par variable sur device | non (sauf cas 4) | oui |
| Paramètres de noyau | inchangés | +48 o par variable (SYCL 1 Ko garanti, HIP 4 Ko) ou tableau en mémoire device |

Ordres de grandeur (RTX 5080, 100 003 points) : cas 5 avec 4 propriétés matériau et 1 variable externe en champs ≈ 8 Mo par appel, ~1-2 ms séquentiel, ~0,3-0,6 ms PCIe : ×5-10 sur J2 GPU, +20-40 % sur Norton GPU. Cas 2-4 ≈ 0. Sur CPU, pas de régression par rapport au chemin par point (qui copie déjà ces valeurs, `Integrate.cxx:154-168`).

**Décision : les deux formes cohabitent dans TFEL.** Le champ par groupe
reste valable. Les champs par variable sont optionnels
(`material_properties_by_variable`, `external_state_variables_by_variable`).
Les deux formes sont converties côté hôte dans la même description
(`VariableFields`), passée par valeur aux noyaux : un seul chemin de code
par point. Pour MGIS, chaque `FieldHolder` devient directement un champ
(uniforme : pas nul ; en place : `{p+b*s, s, 1}`), sans espace de travail
ni copie. Le planificateur du 2.2 n'est plus nécessaire pour les
propriétés matériau et les variables externes.

## 3. Résidence mémoire

### 3.1 Exigences par backend

| Backend | Mémoire acceptée | Ressource conseillée | `queue` |
|---|---|---|---|
| stlpar gcc/clang/icpx, openmp-host-gcc, kokkos-host | hôte | défaut | nul |
| sycl (CPU ou GPU) | USM ; hôte accepté avec transferts automatiques | `sycl::malloc_shared(q)` | `sycl::queue*` obligatoire |
| cuda, hip | device ou managée | `cudaMallocManaged` / `hipMallocManaged` | stream ou nul |
| kokkos (CUDA) | accessible depuis l'espace d'exécution | `kokkos_malloc<SharedSpace>` | `ExecSpace*` ou nul |
| stlpar nvc++ GPU, openmp nvc++ | managée (MGIS compilé par g++ : son `malloc` n'est pas managé) | `cudaMallocManaged` | nul |

### 3.2 Phase 1 (sans dépendance GPU dans MGIS)

```cpp
// MaterialStateManagerInitializer
std::pmr::memory_resource* memory_resource = nullptr;  // nul : new_delete_resource()
// MaterialDataManagerInitializer : propagée à s0/s1
std::pmr::memory_resource* memory_resource = nullptr;
MaterialDataManager(const Behaviour&, const size_type, std::pmr::memory_resource&);
std::pmr::memory_resource& getMemoryResource() const noexcept;
```

- Stockage privé en `std::pmr::vector<real>` ; vues publiques (`std::span`) inchangées : aucune rupture d'API.
- `MutableFieldHolder::value` inchangé (mfem-mgis fait `std::get<std::vector<real>>`) : vecteurs LOCAL copiés dans l'espace de travail sauf si la bibliothèque accepte la mémoire hôte ; mémoire EXTERNAL utilisée en place.
- Mémoire **accessible depuis l'hôte** exigée en phase 1 (MGIS lit et écrit sur l'hôte : F = 1, `update`/`revert`, `updateValues`, extraction, HDF5, chemin par point).
- Adaptateurs optionnels « header-only » : `MGIS/Memory/{CUDAManagedMemoryResource, HIPManagedMemoryResource, KokkosMemoryResource<Space>, SYCLSharedMemoryResource(queue&)}.hxx`.

### 3.3 Phase 2

Ressource étendue (`isHostAccessible`, `copy`, `fill`, `prefetch`) ; `update`/`revert`/`updateValues` sur le device ; fonctions purement hôte en erreur si la mémoire n'est pas accessible ; tableau propriétaire au lieu de `std::pmr::vector`. Le cas 5 devient impossible avec des sources uniquement device, d'où l'option A.

## 4. API

### 4.1 Types (`MGIS/Behaviour/BatchIntegrate.hxx`)

```cpp
struct BatchExecutionPolicy { void* queue = nullptr; };
struct BatchIntegrationOptions {
  IntegrationType integration_type = IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR;
  bool compute_speed_of_sound = false;
  real rdt = 10;
  bool store_statuses = true;
  bool store_time_step_scaling_factors = false;
};
struct MGIS_EXPORT BatchIntegrationResult : BehaviourIntegrationResult {
  size_type number_of_failures = 0;
  size_type number_of_unreliable_results = 0;
};
```

### 4.2 Fonctions

```cpp
MGIS_EXPORT BatchIntegrationResult integrate(const BatchExecutionPolicy&, MaterialDataManager&,
                                             const BatchIntegrationOptions&, const real dt);
MGIS_EXPORT BatchIntegrationResult integrate(const BatchExecutionPolicy&, MaterialDataManager&,
                                             const BatchIntegrationOptions&, const real dt,
                                             const size_type b, const size_type e);
MGIS_EXPORT int integrate(const BatchExecutionPolicy&, MaterialDataManager&, const IntegrationType, const real dt);
MGIS_EXPORT std::optional<BatchIntegrationResult> integrate(Context&, const BatchExecutionPolicy&,
    MaterialDataManager&, const BatchIntegrationOptions&, const real) noexcept;
MGIS_EXPORT BatchIntegrationResult executeInitializeFunction(const BatchExecutionPolicy&,
    MaterialDataManager&, std::string_view, std::span<const real> inputs = {});
MGIS_EXPORT BatchIntegrationResult executePostProcessing(const BatchExecutionPolicy&,
    std::span<real> outputs, MaterialDataManager&, std::string_view);
// rotations : surcharges de Behaviour.hxx:290-507 avec BatchExecutionPolicy en tête
```

Plage `[b, e)` (champs décalés, rassemblement limité), `m.n == 0` : succès immédiat ; fonctions d'initialisation à `dt = 0` ; matrices de rotation construites dans l'espace de travail (9 ou 9n valeurs).

### 4.3 Statuts

- `status` ∈ {1, 0, -1} → `exit_status` ; `time_step_increase_factor = s.rdt` (minimum sur le lot) ; nombres d'échecs recopiés ; `n` : premier indice à -1 (sinon dernier à 0) via `integration_statuses` ; `error_message = s.msg` (messages par point perdus).
- **-2 et -3 levés comme erreurs** (variante `Context` : enregistrés), pour ne pas réduire le pas en boucle sur un bogue.
- Initialisation et post-traitements : 1 = succès en `generic-parallel`, 0 en `generic` ; traduits en 1/-1.
- Tous les points sont intégrés (le chemin par point s'arrête au premier échec).

### 4.4 Facteur d'échelle du pas de temps

Le chemin par point part de 10 (`Integrate.cxx:384`) ; `MaterialDataManager::rdt` est documenté mais **jamais lu**. Défaut 10 pour des résultats identiques ; question ouverte.

### 4.5 Basculement automatique

Les surcharges existantes (`integrate`, `executeInitializeFunction`, `executePostProcessing`, rotations) passent par le lot si `batch_*` est non nul (politique par défaut, erreur -2 explicite pour SYCL sans file) ; surcharges `ThreadPool&` : pool ignoré ; `integrate_debug` : erreur ou bibliothèque `generic` compagnon (question ouverte) ; la nouvelle API sur un comportement `generic` se rabat sur une boucle par point.

## 5. Tests (MGIS, `tests/`)

Compilation des `.mfront` de test avec `--interface=generic-parallel` par configuration (`stlpar-parunseq-gcc` toujours actif, option `MGIS_GENERIC_PARALLEL_CONFIGURATIONS`), bibliothèques `generic` comme référence.

1. `GenericParallelLoadTest` : interface, backend, pointeurs, métadonnées, rotations, API par point en erreur.
2. `BatchIntegrateTest` : deux gestionnaires (`generic` et `generic-parallel`), états aléatoires, deux pas, trois types d'intégration ; s1, K, `exit_status`, `rdt` à 1e-12 ; lois Plasticity, Norton, Gurson, TensorialExternalStateVariableTest, ComputeSpeedOfSound, ImplicitNortonWithEnergies, FiniteStrainSingleCrystal, `ode_rk54`.
3. `BatchDataMappingTest` : champs uniformes, en place et mixtes pour les propriétés matériau et les variables externes.
4. Plages `[b, e)`.
5. Échecs (`NortonFailure.mfront`).
6. Erreurs (variable non définie, SYCL sans file, pas de confusion avec -1).
7. Fonctions d'initialisation et post-traitements.
8. Rotations.
9. Ressource mémoire (compteur d'allocations, non-régression de toute la suite).
10. Variantes GPU (adaptateurs managés).
11. Disposition mémoire (`static_assert` entre copie C et `mfront_gpb_BatchData`).
12. Basculement (`IntegrateTest*` rejoués sur `generic-parallel`).

## 6. Plan, risques, questions

### 6.1 Étapes

1. (TFEL, optionnel) symboles `_mfront_gpb_backend`, `_mfront_gpb_host_memory_accessible`, `_mfront_gpb_api_version`.
2. MGIS : en-tête C `BatchBehaviourFctPtr.hxx` ; test 11.
3. MGIS : interface et backend dans `BehaviourDescription`, membres lot de `Behaviour`, fonctions lot de `LibrariesManager`, branches de `load`, gardes du chemin par point ; test 1.
4. MGIS : `std::pmr` ; test 9 ; suite complète.
5. MGIS : `BatchIntegrate` (politique, options d'exécution, résultat, champs par variable) ; tests 2-6 avec stlpar-gcc.
6. MGIS : basculement des surcharges existantes ; test 12.
7. MGIS : initialisation et post-traitements en lot ; test 7.
8. MGIS : rotations en lot ; test 8.
9. MGIS : adaptateurs mémoire ; configurations GPU des tests (test 10).
10. Documentation, `DESIGN.md`, `JOURNAL.md`.

### 6.2 Risques

Noms identiques et signatures différentes (gardes) ; dérive de la copie C (test 11, version) ; coût du cas 5 sur GPU ; migrations de pages en mémoire managée ; spans EXTERNAL non accessibles par le device (faute non détectable) ; sémantique (tous les points intégrés, messages perdus, pas d'`integrate_debug`) ; SYCL (file obligatoire) ; Kokkos partagé et initialisé par l'appelant ; nvc++ (mémoire `cudaMallocManaged`) ; espace de travail non sûr entre threads ; tailles de K en grandes déformations ; paramètres en deux copies dans les tests.

### 6.3 Questions ouvertes

1. B maintenant et A (ABI de `BatchData`) plus tard sur mesure, ou A tout de suite ?
2. Ajouter les symboles TFEL (backend, accessibilité de la mémoire hôte, version) ?
3. rdt initial : garder 10 ou lire `MaterialDataManager::rdt` dans les deux chemins ?
4. Basculement automatique des surcharges existantes vers le lot ? Repli de la nouvelle API sur la boucle par point pour un comportement `generic` ?
5. Bibliothèque `generic` « compagnon » (pour `integrate_debug` et les messages par point) ?
6. Mémoire : phase 1 `std::pmr` puis phase 2 device ? Adaptateurs GPU dans MGIS (header-only) ou dans mfem-mgis ?
7. `queue` seulement dans `BatchExecutionPolicy`, ou aussi membre par défaut du `MaterialDataManager` ?
8. Liaisons C, Python et Fortran de l'API lot : en D1 ou plus tard ?
9. Nom de l'API : surcharge de `integrate` avec `BatchExecutionPolicy`, ou `integrateBatch` ?
10. Ajouter `setOutOfBoundsPolicy(const Behaviour&, ...)` à MGIS ?

### Fichiers clés

- `mgis/src/Integrate.cxx`
- `mgis/src/Behaviour.cxx`
- `mgis/src/LibrariesManager.cxx`
- `mgis/include/MGIS/Behaviour/MaterialDataManager.hxx`
- `tfel/mfront/include/MFront/GenericParallelBehaviour/BatchData.h`
