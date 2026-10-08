# D1, MGIS : intégration d'un lot avec une bibliothèque `generic-parallel`

Statut : plan révisé le 2026-10-08, avant tout code. Bases : MGIS master
`818f0f9` (`mgis/`, à jour avec l'amont), TFEL du proto (`tfel/`, R1 à R6).
La version du 2026-10-06 (analyse détaillée de MGIS, numéros de lignes)
reste dans l'historique git de ce fichier.

## 1. Ce que TFEL fournit maintenant

- Contrat C versionné : `MFront/GenericParallelBehaviour/BatchData.h`
  (`MFRONT_GPB_API_VERSION`) pour les lois et les modèles,
  `MFront/GenericParallelMaterialProperty/MaterialProperty.h`
  (`MFRONT_GPMP_API_VERSION`) pour les propriétés matériau.
- Symboles de chaque bibliothèque : `<b>_mfront_interface` =
  `GenericParallel`, `<b>_generic_parallel_backend`,
  `<b>_generic_parallel_api_version` (lois, modèles et propriétés
  matériau).
- Propriétés matériau et variables externes données groupe par groupe ou
  variable par variable : chaque `FieldHolder` de MGIS devient un champ,
  sans copie (uniforme : pas nul ; tableau : `{p + b * s, s, 1}`).
- Lot vide accepté même avec des pointeurs nuls ; options de calcul
  invalides refusées pour le lot (`-2`) ; statut de sortie toujours
  complet.
- Exécution : file ou flux (`queue`), statistiques asynchrones
  (`mfront_gpb_BatchStatistics`, fonctions C de décodage), file SYCL
  ordonnée exigée en asynchrone, pile des fils GPU réglée par la
  bibliothèque.
- Données partagées (`@SharedData`) copiées par la bibliothèque, une fois
  par GPU : rien à faire dans MGIS, plusieurs GPU par processus possibles.
- Paramètres : réglés dans la bibliothèque `generic-parallel` elle-même
  (`<b>_setParameter`), lus sur l'hôte à chaque appel.

## 2. Ce que MGIS fournit déjà et qui simplifie tout

- **Stockage fourni par l'appelant.** `MaterialDataManagerInitializer` et
  `MaterialStateManagerInitializer` acceptent des tableaux externes pour
  les gradients, les forces, les variables internes, les énergies, `K` et
  la vitesse du son ; les `FieldHolder` des propriétés matériau et des
  variables externes acceptent des `std::span`. Un code qui alloue lui-même
  en mémoire GPU (mfem-mgis avec la mémoire de MFEM) peut donc déjà tout
  placer : **aucune nouvelle API mémoire n'est nécessaire pour mfem-mgis**.
- Stockage interne privé (`std::vector`) : il peut changer d'allocateur
  sans rupture d'API.
- `LibrariesManager::getInterface` existe (jamais appelé).
- `MGIS_HOST_DEVICE` et un module `mgis::function` avec des algorithmes
  CUDA existent déjà : D1 n'en dépend pas, mais n'empêche pas de s'y
  raccorder plus tard.

Ce qu'attend mfem-mgis (D2) : aujourd'hui l'intégration est faite point
par point au milieu de la boucle d'assemblage
(`BehaviourIntegratorBase::performsLocalBehaviourIntegration`). Sur GPU,
l'assemblage devient trois phases : calcul des gradients en tous les points
(sur le GPU), intégration du lot, assemblage. MGIS doit donc :

1. intégrer un lot (ou une plage `[b, e)`) sur des tableaux qui restent sur
   le GPU ;
2. ne jamais lire ni écrire ces tableaux sur l'hôte pendant l'intégration ;
3. permettre de soumettre plusieurs matériaux avant une seule
   synchronisation ;
4. mettre à jour l'état (`update`, `revert`) sans rapatrier les données.

## 3. Principes

1. **MGIS reste une bibliothèque hôte sans dépendance GPU** (compilée par
   g++). Elle charge, vérifie, décrit les données et appelle ; le calcul
   est dans la bibliothèque générée.
2. **Zéro copie** entre `MaterialDataManager` et `BatchData`, propriétés
   matériau et variables externes comprises.
3. **MGIS n'alloue jamais dans la mémoire d'un GPU.** La mémoire vient de
   l'appelant (tableaux externes), ou d'une ressource mémoire qu'il donne
   (section 5, confort pour les codes simples et les tests).
4. **Le chemin de lot ne touche pas au contenu des tableaux sur l'hôte**
   (pas de copie, pas de lecture des statuts par point) : la mémoire peut
   être purement GPU.
5. **Erreurs explicites au lieu de comportements indéfinis.** MGIS connaît
   le backend (symbole) et l'emplacement des tableaux, déclaré à la
   construction (section 4.2) : un backend GPU avec de la mémoire hôte
   ordinaire est refusé avec un message clair. Les points d'entrée par
   point refusent une bibliothèque `generic-parallel`, et inversement.
6. **Basculement automatique.** Les fonctions existantes (`integrate`,
   `executeInitializeFunction`, `executePostProcessing`, rotations)
   reconnaissent une bibliothèque `generic-parallel` : un code existant
   fonctionne sans changement dès que sa mémoire convient au backend.
7. **Statuts.** `1`, `0`, `-1` : statut d'intégration habituel ; `-2`
   (arguments) et `-3` (erreur du backend) : exceptions, ou erreurs du
   `Context` dans les variantes `noexcept`. Un bogue n'est jamais pris pour
   un échec de convergence (qui ferait réduire le pas sans fin).

## 4. API

### 4.1 Chargement

- `BehaviourDescription` : `interface` (`""`/`Generic` ou
  `GenericParallel`), `backend` ; version d'API vérifiée au chargement.
- `Behaviour` : pointeurs de lot (`batch_b`, fonctions d'initialisation et
  post-traitements de lot, trois rotations). Pour une bibliothèque
  `generic-parallel`, les pointeurs par point restent nuls.
- `LibrariesManager` : `getBatchBehaviour`, `getBatchInitializeFunction`,
  `getBatchPostProcessing`, `getBatchRotation` ; les fonctions par point
  vérifient l'interface.
- Modèles (`libGenericParallel<...>-<backend>.so`) : même chemin, sans
  travail supplémentaire.
- MGIS garde sa propre copie des structures C (comme pour l'interface
  `generic`, sans dépendance de compilation à TFEL) ; un test vérifie la
  disposition mémoire contre `BatchData.h`.

### 4.2 Emplacement des données et politique d'exécution

L'emplacement des données est une propriété des données, déclarée à la
construction ; la politique d'exécution ne décrit que l'exécution.

```cpp
namespace mgis::behaviour {
  //! \brief location of the arrays of a material data manager
  enum struct MemoryLocation {
    HOST,    //!< ordinary memory of the host (default)
    DEVICE,  //!< memory of the device, not accessible by the host
    SHARED   //!< accessible by the host and the device (managed, unified)
  };
  // MaterialDataManagerInitializer: MemoryLocation location = HOST;
  //! \brief execution of a batch (generic-parallel libraries)
  struct BatchExecutionPolicy {
    //! \brief queue, stream or execution space (backend dependent)
    void* queue = nullptr;
    //! \brief see mfront_gpb_ExecutionOptions
    bool preserve_stack_size = false;
  };
}
```

MGIS croise l'emplacement et le backend : mémoire `HOST` refusée par un
backend GPU (sauf SYCL avec transferts automatiques) ; fonctions hôte
(extraction, sauvegarde, valeurs locales) refusées pour `DEVICE`.

`MaterialDataManager` porte une politique par défaut
(`m.batch_execution_policy`), utilisée par le basculement automatique ;
les surcharges qui prennent une politique la remplacent pour un appel.

### 4.3 Fonctions

```cpp
// synchronous integration of all the points, or of [b, e)
BehaviourIntegrationResult integrate(const BatchExecutionPolicy&,
                                     MaterialDataManager&,
                                     const BehaviourIntegrationOptions&,
                                     real dt);
BehaviourIntegrationResult integrate(const BatchExecutionPolicy&,
                                     MaterialDataManager&,
                                     const BehaviourIntegrationOptions&,
                                     real dt, size_type b, size_type e);
// Context variants (noexcept), initialize functions, post-processings
// and rotations with the same pattern
```

- Résultat : `BehaviourIntegrationResult` existant, enrichi du nombre
  d'échecs et de résultats peu fiables. Le premier point en échec (`n`)
  n'est connu que si l'appelant fournit un tableau de statuts.
- Statuts et facteurs de pas par point : tableaux optionnels du
  `MaterialDataManager`, fournis comme `K`
  (`useExternalArrayOfIntegrationStatuses`, `...OfTimeStepScalingFactors`).
- Facteur de pas initial : le même que le chemin par point (résultats
  identiques).
- `ThreadPool&` : ignoré pour une bibliothèque `generic-parallel` (le
  parallélisme est celui du backend).
- `integrate_debug` : refusé pour une bibliothèque `generic-parallel`
  (utiliser la bibliothèque `generic` de la même loi).

### 4.4 Asynchrone (pour plusieurs matériaux)

```cpp
// submits the integration and returns at once (CUDA, HIP, Kokkos, SYCL
// with an in-order queue); the statistics must be accessible by the device
void submit(const BatchExecutionPolicy&, MaterialDataManager&,
            const BehaviourIntegrationOptions&, real dt,
            mfront_gpb_BatchStatistics&);
// once the queue is synchronized by the caller
BehaviourIntegrationResult getIntegrationResult(
    const mfront_gpb_BatchStatistics&);
```

### 4.5 Mise à jour de l'état

`update` et `revert` copient `s1` dans `s0` (ou l'inverse) sur l'hôte.
Le propriétaire de la mémoire fournit la primitive de copie, MGIS garde
l'algorithme (quels blocs copier) :

```cpp
// copy: callable copying n reals (the copy of MFEM for mfem-mgis)
template <typename CopyFunction>
void update(MaterialDataManager&, const CopyFunction&);
template <typename CopyFunction>
void revert(MaterialDataManager&, const CopyFunction&);
```

Sans primitive, `update`/`revert` copient sur l'hôte et refusent la
mémoire `DEVICE`.

## 5. Ressource mémoire (confort)

Pour les codes qui laissent MGIS allouer : `std::pmr::memory_resource*`
optionnel dans les initialiseurs, utilisé pour tout le stockage interne
(états, `K`, vitesse du son, statuts). Ressources fournies en en-têtes
seuls, compilées par le code appelant (MGIS n'en dépend pas) : mémoire
managée CUDA et HIP, USM partagée SYCL, espace partagé Kokkos. Les
`FieldHolder` gardent `std::vector` (API publique utilisée par mfem-mgis) :
avec un backend GPU, les propriétés matériau et variables externes doivent
être uniformes ou des `std::span` en mémoire accessible, sinon erreur.

## 6. Tests

Lois de test de MGIS compilées avec `generic` (référence) et
`generic-parallel` pour les configurations choisies à la configuration de
MGIS (`stlpar-parunseq-gcc` par défaut ; `cuda-nvcc`, `kokkos`,
`sycl-default-acpp` en option). Les tests GPU allouent en mémoire managée
par l'API C du runtime (`cudaMallocManaged` appelée depuis du code g++, sans
nvcc) et tournent **sans HMM**, pour attraper tout accès à la mémoire de
l'hôte.

1. Chargement : interface, backend, version, pointeurs, refus croisés des
   API par point et par lot.
2. Disposition mémoire des structures copiées.
3. Intégration : lot contre boucle par point sur des états identiques
   (tolérance d'arrondi, chargements loin des limites de convergence) ;
   trois types d'intégration ; lois de MGIS (Plasticity, Norton, Gurson,
   variables externes tensorielles, vitesse du son, énergies, grandes
   déformations, Runge-Kutta, monocristal).
4. Champs : propriétés matériau et variables externes uniformes, en place
   et mixtes ; plages `[b, e)` ; lot vide.
5. Échecs et erreurs : échec en quelques points, file SYCL absente,
   mémoire hôte refusée sur GPU, `-2`/`-3` jamais pris pour `-1`.
6. Fonctions d'initialisation, post-traitements, rotations.
7. Basculement : les tests d'intégration existants rejoués sur
   `generic-parallel`.
8. Asynchrone et mise à jour sur le GPU.
9. Ressource mémoire : suite complète avec la ressource managée.

## 7. Étapes

Chaque étape : code, tests, suite complète de MGIS, commit.

1. Chargement, détection, gardes ; structures C et test de disposition.
2. Intégration synchrone (politique, correspondance des données,
   vérifications de mémoire, statuts, plages) et basculement automatique ;
   tests CPU (`stlpar`, `openmp`) puis GPU (`cuda`, mémoire managée).
3. Fonctions d'initialisation, post-traitements, rotations.
4. Ressource mémoire et ressources en en-têtes seuls.
5. Asynchrone, `update`/`revert` sur le GPU.
6. Propriétés matériau `generic-parallel` dans `mgis::material_property`
   (évaluation d'un tableau de points), si mfem-mgis en a besoin.
7. Documentation. Liaisons C, Python et Fortran : plus tard.

Les étapes 1 à 3 donnent un MGIS utilisable par mfem-mgis sur CPU et sur
GPU (mémoire fournie par l'appelant) ; 4 et 5 servent la performance et le
confort ; 6 dépend de D2.

## 8. Risques

- Noms de points d'entrée identiques, signatures différentes : gardes et
  test 1.
- Dérive des structures copiées : version vérifiée au chargement, test 2.
- Mémoire managée touchée par l'hôte entre deux pas : migrations de pages
  (performances) ; d'où la copie sur le GPU de la section 4.5 et la règle
  « le chemin de lot ne touche pas aux tableaux ».
- Sémantique différente du chemin par point : tous les points sont
  intégrés (le chemin par point s'arrête au premier échec), messages par
  point non transmis.
- Pile GPU des polycristaux (près de 15 Go réservés) : un seul processus
  par GPU pour ces lois.

## 9. Décisions à prendre

1. Basculement automatique des fonctions existantes (recommandé : oui).
2. Politique par défaut portée par `MaterialDataManager` (recommandé : oui).
3. Asynchrone et copie sur le GPU dès D1 (étape 5) ou après les premiers
   essais de mfem-mgis (recommandé : après l'étape 3, selon les mesures).
4. Propriétés matériau `generic-parallel` dans MGIS (étape 6) : selon D2.
