# D2, MGIS : API sur des vues (proposition, version 4)

> Document historique. Proposition de conception écrite avant le travail
> sur les vues de MGIS, après trois relectures (annexe). L'état actuel est
> décrit dans [mgis.md](../mgis.md#vues-et-rangements), le cheminement
> dans le [journal](../../JOURNAL.md).

Écarts avec la proposition, constatés en codant :

- les ressources mémoire CUDA, HIP et Kokkos ne pouvaient pas servir de
  ressource amont à un pool (alignement refusé au-delà de celui de leur
  mémoire) : corrigé par sur-allocation, sans quoi le remède de 3.5 ne
  marchait pas ;
- l'extraction d'une variable interne sur une vue d'état passe par des
  opérations (`extractInternalStateVariable(sortie, état, nom, ops)`) et
  marche donc aussi sur GPU. Des accesseurs par variable sur une vue
  d'état (`getGradientView`, etc.) s'ajoutent, comme sur le gestionnaire ;
- les surcharges des rotations et de `evaluate` sur tableaux construisent
  des vues et appellent les fonctions sur vues : un seul chemin de code ;
- une copie d'une vue dans elle-même ne fait rien (les copies de la
  bibliothèque standard et des runtimes GPU interdisent les recouvrements).

Statut : proposition, aucun code écrit. Première étape du
[plan pour mfem-mgis](../docs/plan-mfem-mgis.md). Cette version intègre
trois relectures indépendantes (version 1 : constats A1 à C5 ; version 2 :
constats N1 à N10 ; version 3 : constats F1 à F6), tous vérifiés dans le
code. Leur traitement est en annexe. La dernière relecture conclut que la
conception tient et qu'aucun point restant n'impose de changer l'API.

## 1. Besoin

mfem-mgis veut ranger les données aux points de Gauss de deux façons, au
choix, sur CPU comme sur GPU :

- **point par point** : les composantes d'un point sont contiguës
  (`valeur(i, j) = v[i * nc + j]`) ;
- **composante par composante** : les valeurs d'une composante sont
  contiguës sur tous les points du matériau (`valeur(i, j) = v[j * n + i]`),
  accès coalescés sur GPU. C'est le rangement *strided coalesced* de
  `mgis::function` et de TFEL (`tfel::math::StridedCoalescedView`).

MGIS doit calculer sur ces données sans copie et sans surcoût, avec toutes
ses fonctions, et finir un pas de temps sur GPU.

**TFEL n'a rien à changer** (vérifié par les deux relectures). Son contrat
décrit chaque groupe par un pointeur et deux pas. Ses noyaux lisent et
écrivent ainsi tous les groupes, les entrées des fonctions
d'initialisation, les sorties des post-traitements et les rotations. Le
backend SYCL copie les champs à pas. Les tableaux de champs par variable
sont lus sur l'hôte pendant l'appel, en asynchrone aussi, pour tous les
backends (foncteurs passés par valeur ; OpenMP et `stlpar` synchrones).

## 2. Ce que montrent MFEM 4.9 et mfem-mgis

- L'ordre des degrés de liberté (`FiniteElementSpaceOrdering`, `byNODES`
  ou `byVDIM`) est indépendant du rangement aux points de Gauss.
  `ElementRestriction` rend toujours des vecteurs élémentaires
  `(ND, VDIM, NE)`.
- Le rangement `byNODES` de `QuadratureInterpolator` est par blocs
  d'élément (`q + NQ * (c + K * e)`). Il ne sert pas tel quel : les
  gradients de la loi demandent un noyau de conversion, et les points d'un
  matériau n'y sont pas contigus. Ce noyau écrit dans le rangement choisi.
- mfem-mgis : `Material` dérive de `MaterialDataManager`. Le chemin actuel
  suppose le point par point à une quarantaine d'endroits et garde ce
  rangement. Le rangement composante par composante sert aux nouveaux
  chemins par lot. mfem-mgis donne ses propriétés matériau et variables
  externes avec la politique `NOUPDATE`.

## 3. Proposition

### 3.1 `FieldView` : un champ à pas, avec ses étendues

Le pendant à deux dimensions de `std::span`, avec le sens de `std::mdspan`
et `layout_stride` (absent de C++20 et de gcc 14) :

```cpp
namespace mgis {
  //! rangement des tableaux d'un gestionnaire ; COMPONENT_BY_COMPONENT est
  //! le rangement « strided coalesced » de mgis::function
  enum struct FieldLayout { POINT_BY_POINT, COMPONENT_BY_COMPONENT };

  //! valeur(i, j) = data()[i * getPointStride() + j * getComponentStride()]
  template <typename ValueType>  // real ou const real
  struct FieldView {
    //! vue vide : jamais un champ uniforme par défaut
    constexpr FieldView() noexcept = default;
    //! préconditions : étendue suffisante, pas cohérents (voir plus bas)
    static constexpr bool checkPreconditions(AbstractErrorHandler&,
                                             std::span<ValueType>, size_type n,
                                             size_type nc, size_type ps,
                                             size_type cs);
    constexpr FieldView(std::span<ValueType>, size_type n, size_type nc,
                        size_type ps, size_type cs);
    template <bool doPreconditionsCheck>
    constexpr FieldView(const PreconditionsCheck<doPreconditionsCheck>&,
                        std::span<ValueType>, size_type n, size_type nc,
                        size_type ps, size_type cs);
    //! conversion gratuite d'une vue modifiable en vue constante
    template <typename OtherValueType>
    requires(std::is_same_v<const OtherValueType, ValueType>)
    constexpr FieldView(const FieldView<OtherValueType>&) noexcept;
    // data(), getNumberOfPoints(), getNumberOfComponents(),
    // getPointStride(), getComponentStride(), empty(), operator()(i, j)
   private:
    std::span<ValueType> values;
    size_type n = 0, nc = 0, point_stride = 0, component_stride = 0;
  };

  FieldView<T> makeFieldView(std::span<T>, FieldLayout, size_type n,
                             size_type nc);                        // dense
  FieldView<T> makePointByPointFieldView(std::span<T>, size_type n,
                                         size_type nc, size_type ps); // ps >= nc
  FieldView<T> makeComponentByComponentFieldView(std::span<T>, size_type n,
                                                 size_type nc, size_type ld); // ld >= n
  FieldView<T> makeUniformFieldView(std::span<T>, size_type n); // pas de point nul
  FieldView<T> getComponents(FieldView<T>, size_type offset, size_type count);
}
```

- Gestion des erreurs comme dans `mgis::function` : `checkPreconditions`
  statique, et constructeur à `PreconditionsCheck` qui appelle
  `check_preconditions<FieldView>` (vérification contournable quand
  l'appelant l'a déjà faite). Sans dériver de `PreconditionsChecker`, qui
  n'a ni constructeur par défaut ni affectation : `FieldView` reste un type
  régulier, comme `std::span`.
- Préconditions : un ensemble vide (`n == 0` ou `nc == 0`) est accepté
  d'abord, sans calcul d'étendue ; sinon l'étendue couvre
  `(n - 1) * ps + (nc - 1) * cs + 1` valeurs.
- **Sorties** (contrôle à l'appel) : aucune valeur ne doit servir à deux
  couples `(i, j)`. Contrôle exact et en temps constant d'une condition
  suffisante : en rangeant les deux dimensions d'étendue supérieure à un
  par pas croissant, le petit pas est non nul et le grand pas est au moins
  le petit pas fois son étendue. Les entrelacements exotiques sont donc
  refusés pour les sorties. Les entrées peuvent être uniformes.
- **Groupes vides** (modèles sans gradients, lois sans variables internes,
  énergies non calculées) : la vue doit être vide, sans autre contrôle.
- Interopérabilité avec `mgis::function`, dans un seul sens : une
  `FunctionView` et une vue *strided coalesced* (à laquelle on ajoute un
  accès public à ses données) se convertissent en `FieldView`. Les vues
  `CoalescedMemoryAccess` (un pointeur par composante) n'ont pas
  d'équivalent : un seul pointeur par champ, comme le contrat de TFEL.

### 3.2 `MaterialDataView` : la description d'un matériau sur des tableaux d'autrui

`MaterialDataManager` possède ses tableaux. `MaterialDataView` décrit des
tableaux de l'appelant, avec n'importe quels pas. Ce sont des types
réguliers (copiables et affectables) : la loi est gardée par pointeur,
comme `BehaviourDataView` garde ses tableaux.

```cpp
namespace mgis::behaviour {
  struct MaterialStateView {
    MaterialStateView(const Behaviour&, size_type n,
                      std::pmr::memory_resource* = nullptr);
    FieldView<real> gradients, thermodynamic_forces,
        internal_state_variables, stored_energies, dissipated_energies;
    // une vue par variable, dans l'ordre de la loi, vides au départ
    std::pmr::vector<FieldView<const real>> material_properties;
    FieldView<const real> mass_density;
    std::pmr::vector<FieldView<const real>> external_state_variables;
    const Behaviour* b;  // jamais nul
    size_type n;
  };
  // par nom, avec variantes Context, comme pour MaterialStateManager
  void setMaterialProperty(MaterialStateView&, std::string_view,
                           FieldView<const real>);
  void setExternalStateVariable(MaterialStateView&, std::string_view,
                                FieldView<const real>);
  void setMassDensity(MaterialStateView&, FieldView<const real>);

  struct MaterialDataView {
    MaterialDataView(const Behaviour&, size_type n, MemoryLocation,
                     std::pmr::memory_resource* = nullptr);
    MaterialStateView s0, s1;
    FieldView<real> K;               // si un opérateur est demandé
    FieldView<real> speed_of_sound;  // si demandée
    //! politique des fonctions qui n'en prennent pas, comme le gestionnaire
    BatchExecutionPolicy batch_execution_policy;
    const Behaviour* b;  // jamais nul
    size_type n;
    MemoryLocation memory_location;
  };
  //! vue sur les tableaux d'un gestionnaire, dans son rangement, avec sa
  //! politique ; les vecteurs par variable sont alloués par la ressource
  MaterialDataView make_view(MaterialDataManager&,
                             std::pmr::memory_resource* = nullptr);
}
```

- Vecteurs `std::pmr` : la copie d'une vue alloue avec la ressource par
  défaut, le déplacement garde la ressource. Une vue construite sur un
  tampon de pile ne doit pas sortir de sa portée, et deux vues de
  ressources différentes ne s'échangent pas (`swap`). Seules les fonctions
  du gestionnaire utilisent un tampon de pile, pour une vue qui ne sort
  pas de l'appel. `make_view` sans ressource utilise la ressource par
  défaut.
- `make_view` vérifie que les valeurs gardées par le gestionnaire sur
  l'hôte (scalaires, stockage local) sont lisibles par les calculs de la
  loi, comme aujourd'hui : une vue obtenue par `make_view` ne contourne
  aucun contrôle.
- Membres publics, comme le gestionnaire. Chaque fonction vérifie la
  cohérence de la vue (même loi et même nombre de points dans la vue et
  ses états, nombre de vues par variable, étendues de chaque champ,
  emplacement). Un contrôle par variable, sur l'hôte, par appel.
- Trois niveaux : `MaterialDataManager` (possède), `MaterialDataView`
  (décrit, dure tout le calcul), et à chaque appel la structure C
  `BatchDataView` (contrat de TFEL). Copier une vue copie des
  descripteurs, jamais des valeurs.

**Fonctions sur une vue.** Les mêmes que sur le gestionnaire. La forme
« sans politique » du gestionnaire, qui bascule entre lot et calcul point
par point, prend sur une vue un **espace de travail explicite** :

| Fonction | Formes sur une vue |
|---|---|
| `integrate` | `(politique, v, opts, dt[, b, e])` et `(politique asynchrone, ...)` : bibliothèque de lots ; `(espace de travail, v, opts, dt[, b, e])` : toute bibliothèque (lot avec la politique de la vue, ou point par point sur l'hôte) ; `integrate_debug(espace de travail, v, ...)` |
| `executeInitializeFunction` | mêmes formes, entrées en `FieldView` ou absentes |
| `executePostProcessing` | mêmes formes, sorties en `FieldView` |
| rotations (gradients, forces, blocs de K) | matrices en `FieldView` ou à partir des axes (`RotationMatrix2D/3D`), en place ou non ; synchrone, asynchrone, sans politique |
| `evaluate` (propriétés matériau) | sortie et entrées en `FieldView`, avec ou sans `OutOfBoundsPolicy` ; synchrone, asynchrone, sans politique |

- Le nombre de points vient des étendues : plus de `n` à part.
- **Un seul chemin de code** : les fonctions du gestionnaire appellent
  celles de la vue. Sans politique : `make_view(m)` et l'espace de travail
  du fil appelant. Avec `ThreadPool` : un espace de travail par fil, sur
  des plages. Ce qui est propre au gestionnaire reste chez lui (allocation
  de K à la demande, contrôle des valeurs gardées sur l'hôte).
- **Espace de travail** (`BehaviourIntegrationWorkSpace`) : un par fil,
  fourni par l'appelant. Il garde un pointeur vers la loi pour laquelle il
  a été construit (vérifié à chaque appel). Il reçoit, alloués une fois à
  sa construction, les tampons d'un point pour `s0` et `s1` : gradients,
  forces, variables internes, énergies, K, plus grande entrée de fonction
  d'initialisation, plus grande sortie de post-traitement.
- Les fonctions du gestionnaire sans politique aiguillent une loi
  `generic-parallel` vers le lot **avant** de prendre un espace de travail,
  comme aujourd'hui : ni verrou ni allocation sur ce chemin.
- **Bibliothèque `generic`** : calcul point par point sur l'hôte. Deux
  versions de la boucle, choisies une fois par appel : pointeurs directs si
  toutes les composantes sont contiguës (le cas actuel, inchangé), sinon
  valeurs rassemblées dans l'espace de travail puis recopiées. Le
  rassemblement reproduit le calcul actuel (options d'intégration écrites
  dans `K[0]`, propriétés et variables externes lues avec leurs pas). Les
  rotations d'une bibliothèque `generic` sur des champs non contigus font
  de même avec un tampon d'un point (sur la pile, sur le tas au-delà).
- Pas de `ThreadPool` sur une vue : le parallélisme CPU passe par les
  backends `openmp` et `stlpar`, ou par l'appelant avec ses espaces de
  travail.
- Les anciennes formes à `IntegrationType`, gardées pour la compatibilité
  (et utilisées par les liaisons C), restent sur le gestionnaire seulement :
  les vues sont une nouvelle API.
- Fonctions hôte sur une vue, comme sur le gestionnaire :
  `extractInternalStateVariable(sortie, état, nom)`,
  `convertFiniteStrainStress` et `convertFiniteStrainTangentOperator`
  (sortie en `FieldView`), lisibles par l'hôte seulement.

### 3.3 Rangement d'un `MaterialDataManager`

```cpp
auto i = MaterialDataManagerInitializer{};
i.layout = FieldLayout::COMPONENT_BY_COMPONENT;  // POINT_BY_POINT par défaut
```

Le rangement est rangé dans `MaterialDataManager::layout` et
`MaterialStateManager::layout` (nouveau paramètre de son constructeur),
comme l'emplacement mémoire. **Tout ce qui est donné par point au
gestionnaire suit ce rangement** :

- tableaux d'état, K, vitesse du son ;
- variables externes à plusieurs composantes données par un `std::span`
  (les propriétés matériau sont scalaires) ;
- entrées des fonctions d'initialisation et sorties des post-traitements
  données par un `std::span` ;
- sorties de `extractInternalStateVariable` et des conversions des grandes
  transformations, qui reçoivent aussi des surcharges en `FieldView`.

Ce qui change dans MGIS :

- **F mis à l'identité** (lois `F-Cauchy`) avec les pas du rangement ;
- `updateValues(o, i)` vérifie que les deux états ont le même rangement ;
- `save` et `restore` écrivent et lisent toujours point par point : une
  reprise ne dépend pas du rangement ;
- accès dans le rangement du gestionnaire : `getGradientsView(s)`,
  `getGradientView(s, nom)`, de même pour les forces, les variables
  internes et les blocs de K ;
- les surcharges de rotation à `std::span` restent point par point (elles
  ne voient pas le gestionnaire) : la doc le dit et renvoie aux surcharges
  en `FieldView` ;
- `MaterialFunctionManager` : constructeur avec initialiseur ; rangé
  composante par composante, il rend des vues *strided coalesced* par de
  nouveaux accesseurs de taille fixe (`getStridedCoalescedGradient<N>`,
  etc.). Ses accesseurs actuels restent point par point (erreur claire
  sinon) : `mgis::function` n'a pas de vue *strided coalesced* de taille
  dynamique ;
- liaisons : le rangement est exposé dans les initialiseurs et par un
  accesseur en C, Fortran et Python. Python rend des tableaux numpy avec
  les pas du rangement (pybind11). FEniCS construit toujours son
  gestionnaire point par point : rien à faire. Julia n'expose aucun
  tableau. Les accesseurs `*_stride` gardent leur sens (nombre de valeurs
  par point).

`update` et `revert` copient des tableaux entiers : ils ne dépendent pas du
rangement.

### 3.4 Fin de pas sur GPU : la même procédure, des copies faites par l'appelant

`update` et `revert` font plus que copier (`updateValues` crée ou efface
des entrées, change un scalaire en tableau, remet la masse volumique à
zéro, respecte `shall_be_updated`). Pour garder exactement ce sens, MGIS
garde sa procédure et confie seulement les copies de valeurs à des
opérations fournies par l'appelant, sur le modèle de
`std::pmr::memory_resource` :

```cpp
namespace mgis {
  struct FieldOperations {
    //! copie (mêmes étendues), éventuellement asynchrone sur la file de
    //! l'appelant ; les deux opérandes sont dans l'emplacement du
    //! gestionnaire ou de la vue
    void copy(FieldView<real>, FieldView<const real>);
    void fill(FieldView<real>, const real);
    virtual ~FieldOperations();
   private:
    // interface publique non virtuelle, comme std::pmr::memory_resource
    virtual void do_copy(FieldView<real>, FieldView<const real>) = 0;
    virtual void do_fill(FieldView<real>, const real) = 0;
  };
  struct HostFieldOperations final : FieldOperations { /* std::copy */ };
}
namespace mgis::behaviour {
  void update(MaterialDataManager&, FieldOperations&);
  void revert(MaterialDataManager&, FieldOperations&);
  void updateValues(MaterialStateManager&, const MaterialStateManager&,
                    FieldOperations&);  // avec variantes Context
  void update(const MaterialDataView&, FieldOperations&);  // états et K
  void revert(const MaterialDataView&, FieldOperations&);
  void updateValues(const MaterialStateView&, const MaterialStateView&,
                    FieldOperations&);  // tableaux d'état
}
```

- **Une seule procédure** : `update(m)` garde son contrôle d'accès par
  l'hôte, puis appelle `update(m, HostFieldOperations)`. De même pour
  `revert(m)` et `updateValues(o, i)`. Le sens actuel est conservé cas par
  cas (redimensionnement, masse volumique, insertions, effacements,
  `shall_be_updated` ; vérifié par la relecture).
- Règle des valeurs : un tableau de l'appelant vers un tableau de
  l'appelant passe par les opérations. Tout le reste se fait sur l'hôte
  par MGIS : valeurs qu'il garde (scalaires, stockage local), et toute
  copie entre une valeur qu'il garde et un tableau de l'appelant, **dans
  les deux sens**. Ces copies mixtes sont refusées, avec un message, si
  l'emplacement du gestionnaire n'est pas lisible par l'hôte. Les
  opérations ne reçoivent donc que des tableaux de l'emplacement du
  gestionnaire : `FieldView` n'a pas besoin de porter d'emplacement.
- Ordre et asynchrone : MGIS fait sa part sur l'hôte pendant l'appel, les
  copies confiées aux opérations peuvent partir plus tard sur la file.
  Préconditions, écrites dans la doc : aucun calcul en attente sur le
  gestionnaire au moment de l'appel ; synchroniser la file avant de lire
  sur l'hôte un tableau copié, ou avant un nouvel `update` ou `revert` si
  des copies mixtes sont possibles (emplacement partagé).
- Vue : tableaux d'état (`updateValues` entre deux états d'une vue) et
  remise à zéro de K. Les propriétés matériau et les variables externes
  d'une vue appartiennent à l'appelant.
- Opérations fournies en en-têtes, compilées par l'appelant avec son outil
  GPU, à côté des ressources mémoire (`MGIS/FieldOperations/`) : CUDA et
  HIP (un flux), SYCL (une file), Kokkos (une instance), OpenMP déporté
  (synchrone, numéro de périphérique). Copie contiguë directe quand les
  deux rangements sont denses, un noyau sinon.
- Pas d'échange des états : il n'équivaut pas à `update` (une prédiction
  n'écrit pas `s1`, les points en échec non plus, mfem-mgis lit `s1` pour
  le résidu de prédiction), pour un gain d'une copie par pas de temps.

### 3.5 Corrections du code existant (étape D1)

- Chaque appel asynchrone alloue son statut par la ressource de la
  politique : avec une ressource managée CUDA, un `cudaMallocManaged` et
  un `cudaFree` par appel. Remède standard, sans nouvelle API : une
  `std::pmr::synchronized_pool_resource` au-dessus de la ressource managée,
  recommandée par la doc, utilisée par les tests et les exemples. La doc
  rappelle la précondition déjà écrite (l'objet en attente vit jusqu'à la
  fin du calcul) : avec un pool, la violer corrompt l'appel suivant au lieu
  d'être lente. Les grands blocs (matrices construites à partir des axes)
  dépassent la taille des blocs du pool : il faut les construire une fois
  avec `buildRotationMatrices`, comme la doc le conseille déjà.
- `evaluate` alloue deux `std::vector` par appel : remplacés par un tampon
  sur la pile.

### 3.6 Coût

- **Noyaux de TFEL inchangés** : seul le calcul des adresses change avec le
  rangement.
- **Appels du gestionnaire : aucune allocation dans le cas courant.**
  `make_view` et les tableaux de champs de la structure C sont construits
  dans un tampon sur la pile (`std::pmr::monotonic_buffer_resource`).
  Aujourd'hui, le gestionnaire alloue quatre `std::vector` par appel.
- **Vue construite par l'appelant** : allocation une fois, à sa
  construction. Rien par appel.
- **Calcul point par point** : rien de plus pour le cas actuel ; tampons
  de l'espace de travail alloués une fois.
- **Contrôles** : un par variable et par champ, sur l'hôte, par appel.
- **Fin de pas** : une copie par tableau d'état et par pas de temps (pas
  par itération), et un appel virtuel par copie.
- **À mesurer, ce n'est pas un surcoût** : point par point contre
  composante par composante, selon la loi et la machine.

## 4. Exemples

mfem-mgis, chemin tout GPU, gestionnaire rangé composante par composante :

```cpp
auto i = MaterialDataManagerInitializer{};
i.layout = FieldLayout::COMPONENT_BY_COMPONENT;
i.location = MemoryLocation::DEVICE;
i.s0.gradients = std::span(g0, n * gsize);  // mémoire de MFEM
// ... tous les tableaux donnés
auto m = MaterialDataManager(b, n, i);
setExternalStateVariable(m.s0, "Temperature", std::span(T0, n),
                         MaterialStateManager::EXTERNAL_STORAGE);
setExternalStateVariable(m.s1, "Temperature", std::span(T1, n),
                         MaterialStateManager::EXTERNAL_STORAGE);
auto p = BatchExecutionPolicy{};
p.queue = stream;
const auto r = integrate(p, m, opts, dt);  // inchangé
auto ops = mgis::field_operations::CUDAFieldOperations{stream};
update(m, ops);  // fin de pas sur le GPU, dans l'ordre du flux
```

Code qui gère lui-même ses données :

```cpp
auto v = MaterialDataView(b, n, MemoryLocation::DEVICE);
v.s0.gradients = makeComponentByComponentFieldView(std::span(g0, gsize * ld), n, gsize, ld);
v.s1.gradients = makeComponentByComponentFieldView(std::span(g1, gsize * ld), n, gsize, ld);
// ... forces, variables internes, K
const auto E = makeUniformFieldView(std::span<const real>(E_d, 1), n);
setMaterialProperty(v.s0, "YoungModulus", E);
setMaterialProperty(v.s1, "YoungModulus", E);
auto p = BatchExecutionPolicy{};
p.queue = stream;
const auto r = integrate(p, v, opts, dt);
```

## 5. Symétrie à vérifier avant de dire « fini »

| Axe | Valeurs |
|---|---|
| Objet | loi, modèle, propriété matériau |
| Fonction | intégration, initialisation, post-traitement, rotations (3 sortes, 3 façons de donner les matrices), évaluation, `update`, `revert`, `updateValues`, extraction, conversions des grandes transformations |
| Support | vue, gestionnaire |
| Rangement | point par point, composante par composante, pas quelconques et `ld > n` (vues) |
| Exécution | sans politique (espace de travail), synchrone, asynchrone, `ThreadPool` (gestionnaire) |
| Bibliothèque | `generic-parallel` (hôte et GPU), `generic` |
| Points | tous, plage `[b, e)`, ensemble vide, groupes vides |
| Emplacement | hôte, GPU, partagé, refus attendus |
| Liaisons et reprises | C, Fortran, Python ; `save` et `restore` |

## 6. Écarté, avec la raison

- **Rangement par blocs d'élément** (`byNODES` de `QuadratureInterpolator`) :
  troisième pas et contrat TFEL v3, sans usage (section 2).
- **Statuts par point** : le contrat de TFEL en a pour les lois, pas pour
  les propriétés matériau, avec des conventions différentes, et le calcul
  point par point s'arrête au premier échec. La limite actuelle reste pour
  les deux. Ajout possible plus tard, ensemble.
- **Échange des états** : section 3.4.
- **`ThreadPool` sur une vue** : section 3.2.
- **Variantes `Context` des fonctions de calcul** : les fonctions de calcul
  du gestionnaire n'en ont pas. Celles qui en ont (affectations,
  `updateValues`) en ont aussi sur les vues.
- **Recouvrement entre deux sorties différentes** : non vérifié (coûteux en
  général), contrat documenté. Le recouvrement d'une sortie avec elle-même
  est vérifié (3.1).
- **Ressource mémoire OpenMP** : OpenMP n'a pas de mémoire partagée
  portable ; les opérations de fin de pas OpenMP existent, pas la
  ressource (situation actuelle inchangée).

## 7. Ce qui revient à mfem-mgis

- un pas de composante dans `PartialQuadratureFunction`, ou les vues
  *strided coalesced* de `mgis::function`, pour les sorties d'un matériau
  rangé composante par composante ;
- des tableaux lisibles par le GPU pour `PartialQuadratureSpace` ;
- en assemblage partiel, la gestion des éléments de chaque matériau par
  l'intégrateur (MFEM ignore les marqueurs).

## 8. Étapes

1. `FieldView`, `MaterialDataView`, fonctions sur les vues, fonctions du
   gestionnaire réécrites dessus, espace de travail explicite, tampons sur
   la pile. Les 45 tests actuels inchangés et verts.
2. Rangement du gestionnaire (3.3), `MaterialFunctionManager`, liaisons.
3. `FieldOperations`, `update`, `revert` et `updateValues` avec opérations,
   en-têtes CUDA, HIP, SYCL, Kokkos et OpenMP déporté. Corrections de 3.5.
4. Tests de chaque case de la matrice de symétrie avec quelques lois
   (élasticité, Norton, loi orthotrope, grandes transformations
   `F-Cauchy`, modèle, loi avec fonction d'initialisation et
   post-traitement, propriété matériau), sur CPU (`stlpar`, `openmp`,
   `generic`) et CUDA, dont `ld > n`, prédiction, fin de pas sur GPU, et
   GPU sans HMM. Résultats identiques bit à bit entre vue et gestionnaire,
   et entre les deux rangements. Une seconde construction de MGIS avec les
   liaisons et HDF5 (aujourd'hui désactivés) teste ce qui y est touché. Doc
   web et Doxygen de MGIS, doc du prototype, journal.
5. Relecture indépendante, corrections, seconde relecture du diff.

## Annexe : traitement des relectures

### Version 1

| Constat | Traitement |
|---|---|
| A1 F initialisé point par point | 3.3 |
| A2 espace de travail d'une vue | argument explicite (3.2) |
| A3 échange non équivalent à `update` | échange retiré (3.4) |
| A4 `update` des propriétés et variables externes | même procédure, copies par opérations, règle écrite (3.4) |
| B1 TFEL et statuts par point | statuts retirés pour les deux (section 6) |
| B2 surcharges manquantes | rotations par axes et en place, `evaluate` sans politique de bornes, politique de la vue, setters et `updateValues` avec `Context` |
| B3 statuts par point | retirés |
| B4 rangement des sorties | tout ce qui est donné par point suit le rangement (3.3) |
| B5 pièges de `FieldView` | vue vide par défaut, étendues, conversion vers `const` |
| B6 allocations | tampons sur la pile, pool pour l'asynchrone (3.5, 3.6) |
| B7 doublons et analogies | trois niveaux (3.2), types réguliers |
| B8 interopérabilité `mgis::function` | un seul sens, accès public aux données, accesseurs nommés |
| B9 code qui suppose le point par point | rangement exposé ou refusé dans chaque liaison (3.3) |
| B10 surcharges en trop | pas de `ThreadPool` sur une vue |
| C1 à C5 | ensemble vide d'abord, recouvrement d'une sortie, rangement dans le gestionnaire, fabriques `make…`, `getComponents(vue, début, nombre)`, `n` tiré des vues |

### Version 2

| Constat | Traitement |
|---|---|
| N1 pointeur d'espace de travail dans la vue | retiré : argument explicite, un par fil, loi vérifiée, tampons dimensionnés à la construction |
| N2 allocation de `make_view` | vecteurs `std::pmr`, tampon sur la pile dans les fonctions du gestionnaire |
| N3 fin de pas contre `updateValues` | même procédure avec opérations de copie, règle des valeurs gardées par MGIS, `updateValues` avec opérations, contrôle du rangement entre états, K pour une vue |
| N4 recouvrement | contrôle exact d'une condition suffisante |
| N5 rotations | surcharges `std::span` documentées point par point ; tampon d'un point pour `generic` |
| N6 groupes vides | règle écrite (3.1) |
| N7 pool | `synchronized_pool_resource`, précondition rappelée, grands blocs |
| N8 affirmations | confirmées vraies |
| N9 types non réguliers | pointeur vers la loi, membres non constants, cohérence vérifiée à chaque appel |
| N10 mineurs | `memory_location`, propriétés matériau scalaires, `PreconditionsChecker`, paramètre du constructeur de `MaterialStateManager`, numpy avec pas, opérations hôte publiques (`HostFieldOperations`), plus de fonction `get…` qui modifie, exemples corrigés |

### Version 3

| Constat | Traitement |
|---|---|
| F1 fin de pas : copie d'une valeur gardée par MGIS vers un tableau de l'appelant, ordre en asynchrone, contrôle d'accès de `update(m)` | copies mixtes sur l'hôte dans les deux sens, refusées hors mémoire lisible par l'hôte ; préconditions d'ordre écrites ; contrôle d'accès gardé (3.4) |
| F2 `FieldView` dérivée de `PreconditionsChecker` non régulière | plus d'héritage, même mécanisme par `check_preconditions` (3.1) |
| F3 `make_view` public contourne `checkHostValues` | contrôle dans `make_view` (3.2) |
| F4 sémantique des vecteurs `std::pmr` | écrite (3.2) |
| F5 asymétries : `updateValues` sur une vue, extraction et conversions sur une vue, formes `IntegrationType` | ajoutées, ou écrites comme réservées au gestionnaire (3.2, 3.4) |
| F6 détails de mise en œuvre : aiguillage avant l'espace de travail, `K[0]`, tampons de `s0` et `s1`, pas des évaluateurs, tampon d'une taille qui dépend de la loi, pointeur vers la loi dans l'espace de travail, interface non virtuelle de `FieldOperations` | écrits (3.2, 3.4) |
