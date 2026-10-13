# MGIS : détails techniques

Ce document décrit ce que le prototype change dans MGIS, et pourquoi. MGIS
charge les bibliothèques générées par MFront et gère les données matériau
(`MaterialDataManager`). Le prototype lui apprend à utiliser les
bibliothèques de l'interface `generic-parallel` : lois, modèles et
propriétés matériau calculés par lots, sur CPU ou sur GPU.

- Base : MGIS master `e8ef4b4` (13/10/2026 ; branche du prototype rebasée sur cette version, partie de `818f0f9` du 02/10/2026).
- Sources modifiées : dépôt `mgis/`, branche `gpu-proto`. Série de patchs :
  `patches/mgis/` (`git am`).
- La documentation de l'API, en anglais, est dans l'arbre de MGIS (voir
  [la dernière section](#documentation-dans-larbre-de-mgis)).

## Sommaire

1. [Ce qui change en bref](#ce-qui-change-en-bref)
2. [Chargement](#chargement)
3. [Intégration d'un lot](#intégration-dun-lot)
4. [Fonctions d'initialisation, post-traitements, rotations](#fonctions-dinitialisation-post-traitements-rotations)
5. [Modèles](#modèles)
6. [Propriétés matériau](#propriétés-matériau)
7. [Mémoire](#mémoire)
8. [Politique d'exécution](#politique-dexécution)
9. [Vues et rangements](#vues-et-rangements)
10. [Erreurs](#erreurs)
11. [Corrections de défauts existants](#corrections-de-défauts-existants)
12. [Limites actuelles](#limites-actuelles)
13. [Exemple](#exemple)

## Ce qui change en bref

- MGIS reconnaît une bibliothèque `generic-parallel` au chargement et
  vérifie sa compatibilité.
- Les fonctions existantes (`integrate`, `executeInitializeFunction`,
  `executePostProcessing`, rotations) basculent seules sur le calcul par
  lot. Le code appelant ne change pas.
- Chaque fonction de lot existe aussi avec une politique d'exécution, qui
  choisit la file (flux CUDA, file SYCL, etc.). Le calcul est synchrone par
  défaut, asynchrone avec une politique asynchrone.
- Les données ne sont jamais copiées : le lot est décrit directement sur
  les tableaux du `MaterialDataManager`.
- Les tableaux peuvent être rangés composante par composante, ou
  appartenir à l'appelant et être décrits par des vues à pas quelconques,
  sur GPU compris. La fin de pas se fait alors sur le GPU (voir
  [Vues et rangements](#vues-et-rangements)).
- L'appelant déclare où sont ses tableaux. MGIS vérifie qu'ils sont
  accessibles par les calculs de la bibliothèque.
- MGIS ne dépend d'aucune bibliothèque GPU. Les tableaux viennent de
  l'appelant, ou d'une ressource mémoire qu'il fournit.
- Profilage : `setProfilingSynchronization` donne une fonction appelée au
  début et à la fin de chaque section chronométrée active, avant la
  lecture du temps. Les calculs lancés de façon asynchrone, par exemple
  sur un GPU, sont alors comptés dans la section qui les lance. Aucune
  fonction par défaut.

## Chargement

Au chargement d'une loi (`mgis::behaviour::load`), MGIS lit les symboles de
description de la bibliothèque. Pour une bibliothèque `generic-parallel`,
il vérifie la version du contrat de données et lit le backend et l'espace
mémoire. Ils sont rangés dans la description de la loi (`interface`,
`backend`, `memory_space`). `isGenericParallelBehaviour` dit si une loi
calcule par lots.

MGIS garde sa propre copie des structures C de TFEL
(`MGIS/Behaviour/BatchDataView.hxx`), comme pour l'interface `generic` :
il ne dépend pas de TFEL à la compilation. Un test vérifie que la
disposition mémoire de ces copies est identique à celle de TFEL.

Avant ce travail, une bibliothèque `generic-parallel` se chargeait sans
erreur. Ses fonctions de lot étaient alors appelées avec les arguments d'un
point : comportement indéfini. Désormais, chaque accès du `LibrariesManager`
refuse l'autre type de bibliothèque avec un message clair.

## Intégration d'un lot

```cpp
auto r = integrate(policy, m, opts, dt);        // tous les points
auto r = integrate(policy, m, opts, dt, b, e);  // points [b, e)
```

La version sans politique (`integrate(m, opts, dt)`) utilise la politique
du gestionnaire (`m.batch_execution_policy`). Les versions avec un
`ThreadPool` font un seul appel de lot : le parallélisme est celui du
backend.

Correspondance des données :

- les gradients, les forces thermodynamiques, les variables internes et les
  énergies sont des champs rangés point par point ;
- les propriétés matériau et les variables externes sont données variable
  par variable : chacune est uniforme ou donnée par point, comme dans le
  `MaterialStateManager` ;
- la masse volumique n'est transmise que si la vitesse du son est
  demandée : c'est la seule à l'utiliser.

Différence avec le calcul point par point : tous les points sont intégrés,
même si certains échouent. Le résultat donne l'indice absolu du premier
point en échec (`n`) et le nombre d'échecs dans le message.

## Fonctions d'initialisation, post-traitements, rotations

Elles suivent le même modèle que l'intégration :

```cpp
executeInitializeFunction(policy, m, name[, inputs][, b, e]);
executePostProcessing(outputs, policy, m, name[, b, e]);
rotateGradients(policy, dest, behaviour, src, r[, location]);  // et les autres rotations
```

Les versions sans politique basculent seules. Pour les rotations, les
matrices peuvent être données par point ou uniformes. Les formes
`RotationMatrix2D` et `RotationMatrix3D` construisent les matrices sur
l'hôte à chaque appel : pour des rotations répétées sur GPU, mieux vaut les
construire une fois avec `buildRotationMatrices`, dans une mémoire
partagée donnée par une ressource, et les donner en tableau. Une rotation
se fait en place si la source et la destination sont le même tableau.

Les entrées des fonctions d'initialisation et les sorties des
post-traitements doivent être placées comme les tableaux du gestionnaire.
`allocatePostProcessingVariables` les alloue sur l'hôte (et refuse une loi
qui ne lit pas la mémoire de l'hôte), ou avec une ressource donnée. Les
rotations, avec ou sans politique, prennent l'emplacement de leurs
tableaux en dernier argument (`HOST` par défaut), comme `evaluate`. Une
ressource nulle désigne partout la ressource par défaut de la
bibliothèque standard.

## Modèles

L'interface `generic-parallel` des modèles exporte les mêmes symboles que
celle des lois. `mgis::model::load` passe par le chargement des lois : un
modèle `generic-parallel` s'utilise donc exactement comme une loi, avec
les mêmes fonctions de lot et les mêmes paramètres (`setParameter`).

## Propriétés matériau

```cpp
const auto mp = mgis::material_property::load(library, "YoungModulus");
const auto inputs = std::array{std::span<const real>(temperatures)};
auto r = evaluate(values, mp, inputs);                    // file par défaut
auto r = evaluate(policy, values, mp, inputs, location);  // file choisie
```

Au chargement, MGIS lit l'interface, la version des fonctions, le backend
et l'espace mémoire. Il charge la fonction d'un ensemble de points
(`batch_fct`) au lieu de celle d'un point (`fct`). `evaluate` évalue une
propriété sur un ensemble de points. Chaque entrée est uniforme (une
valeur) ou donnée par point. Comme `integrate`, `evaluate` fait un seul
appel pour une propriété `generic-parallel` et calcule point par point sur
l'hôte sinon.

Le résultat (`MaterialPropertyEvaluationResult`) est le même pour les deux
interfaces, comme `BehaviourIntegrationResult` pour les lois : statut et
rang de la variable hors bornes (sens de l'interface `generic`), indice du
premier point en échec, message. Une propriété `generic-parallel` ne rend
que les statuts `1`, `0`, `-1` et `-4`. Les autres (arguments invalides,
erreur du backend) lèvent une exception, dès la soumission pour une
évaluation asynchrone. Le calcul point par point s'arrête au
premier échec, le lot calcule tous les points : l'indice rapporté est le
même. Au chargement, MGIS lit aussi la fonction de compte rendu des
statistiques (`batch_report_fct`).

## Mémoire

### Emplacement des tableaux et espace mémoire

L'appelant déclare l'emplacement de ses tableaux (`mgis::MemoryLocation`) :

- `HOST` : mémoire ordinaire de l'hôte ;
- `DEVICE` : mémoire d'un GPU, que l'hôte ne peut pas lire ;
- `SHARED` : mémoire lisible par les deux (mémoire managée CUDA ou HIP,
  mémoire partagée SYCL, `Kokkos::SharedSpace`).

Pour un `MaterialDataManager`, l'emplacement est donné par son
initialiseur (`MaterialDataManagerInitializer::location`). Les valeurs
données par un `std::span` (stockage externe) et les tableaux passés aux
fonctions de lot du gestionnaire doivent avoir le même emplacement. Pour
une propriété matériau et pour une rotation, l'emplacement est un
argument.

MGIS compare cet emplacement à l'espace mémoire de la bibliothèque,
exporté par TFEL (`isMemoryLocationAccessible`) :

| Espace mémoire de la bibliothèque | `HOST` | `DEVICE` | `SHARED` |
|---|---|---|---|
| `host` (et bibliothèques `generic`) | oui | non | oui |
| `device` | non | oui | oui |
| `any` | oui | oui | oui |

MGIS ne connaît donc aucune liste de backends. Un nouveau backend de TFEL
n'a besoin d'aucun changement dans MGIS.

### Fonctions qui lisent les tableaux sur l'hôte

Certaines fonctions de MGIS lisent les tableaux sur l'hôte : `update`,
`revert` et `updateValues` sans opérations de copie (`FieldOperations`,
voir [Vues et rangements](#vues-et-rangements)), l'extraction des variables internes, `save`,
`restore`, les conversions des grandes transformations et le calcul point
par point. Elles refusent les tableaux placés sur un GPU au lieu de
planter.

De même, les valeurs gardées par MGIS lui-même (valeur uniforme, stockage
local `LOCAL_STORAGE`) sont en mémoire hôte. Elles sont refusées quand les
calculs ne peuvent pas lire cette mémoire. Il faut alors les donner par un
`std::span` sur une mémoire accessible, en stockage externe. Une valeur
uniforme peut être donnée par un `std::span` d'une seule valeur.

### Ressources mémoire

Un `MaterialDataManager` alloue lui-même les tableaux que son initialiseur
ne donne pas. Par défaut, il les alloue sur l'hôte. L'initialiseur peut
désormais donner une ressource mémoire standard
(`std::pmr::memory_resource*`) : MGIS alloue alors avec elle. MGIS
initialise ces tableaux sur l'hôte : la ressource doit donc allouer de la
mémoire lisible par l'hôte, c'est-à-dire partagée. Si les tableaux sont
sur un GPU seul, l'appelant les donne tous.

MGIS fournit quatre ressources, en en-têtes que le code appelant compile
avec son propre outil GPU (`MGIS/MemoryResources/`) :

| En-tête | Mémoire |
|---|---|
| `CUDAManagedMemoryResource.hxx` | mémoire managée CUDA |
| `HIPManagedMemoryResource.hxx` | mémoire managée HIP |
| `SYCLSharedMemoryResource.hxx` | mémoire partagée SYCL, d'une file donnée |
| `KokkosSharedSpaceMemoryResource.hxx` | `Kokkos::SharedSpace` |

Coût : un pointeur par tableau et un appel virtuel par allocation. Rien par
calcul.

## Politique d'exécution

`mgis::BatchExecutionPolicy` est commune aux lois, aux modèles et aux
propriétés matériau :

- `queue` : file d'exécution, de type propre au backend (`cudaStream_t`,
  `hipStream_t`, pointeur vers une `sycl::queue` ou vers une instance de
  l'espace d'exécution Kokkos). Nul désigne la file par défaut.
- `preserve_stack_size` : interdit de modifier la taille de pile des fils
  du GPU.

Avec cette politique, le calcul est synchrone : la fonction rend la main
une fois le calcul fini. C'est le défaut.

Les fonctions sans politique utilisent celle du gestionnaire
(`batch_execution_policy`), ou une politique par défaut pour les rotations
et les propriétés matériau. La politique par défaut donne une file nulle.
SYCL n'a pas de file par défaut : il faut alors passer une politique, ou
donner une file au gestionnaire.

## Calcul asynchrone

Chaque fonction de lot existe aussi avec une politique asynchrone
(`mgis::AsynchronousBatchExecutionPolicy`) : intégration, fonctions
d'initialisation, post-traitements, rotations, évaluation des propriétés
matériau. Rien ne change pour le code qui ne s'en sert pas.

```cpp
auto p = AsynchronousBatchExecutionPolicy{};
p.queue = stream;
p.resource = &resource;  // mémoire managée
p.location = MemoryLocation::SHARED;
const auto pending = integrate(p, m, opts, dt);
// ... autres calculs soumis au flux, par exemple d'autres matériaux
cudaStreamSynchronize(stream);
const auto r = pending.getResult();
```

- La fonction rend la main dès le calcul soumis. Elle rend un calcul en
  attente : `PendingBatchComputation` pour les lois et les modèles,
  `PendingBatchEvaluation` pour les propriétés matériau.
- Une fois la file synchronisée par l'appelant, `getResult` rend le
  résultat de la fonction synchrone, indice du premier point en échec
  compris. MGIS ne peut pas synchroniser la file lui-même : il ne dépend
  d'aucune bibliothèque GPU.
- Le calcul en attente garde la mémoire écrite ou lue par le GPU : le
  statut du calcul (statistiques des lois, statuts des propriétés
  matériau) et, pour les rotations, les matrices construites à partir des
  axes du matériau. Il ne doit pas être détruit avant la fin du calcul.
  Il peut être déplacé (par exemple dans un `std::optional` ou un
  conteneur), mais pas affecté : l'affectation recopierait le statut si
  les ressources mémoire différaient, pendant que le GPU l'écrit.
- Ce statut est initialisé et lu par l'hôte, écrit par le GPU. La
  politique donne la ressource mémoire qui l'alloue et son emplacement.
  MGIS vérifie qu'il est accessible par les deux : mémoire de l'hôte pour
  les calculs sur l'hôte, mémoire partagée sinon.

Les deux modes partagent tout le code, sauf l'obtention du statut : en
synchrone, la fonction générée le rend ; en asynchrone, il est lu dans les
statistiques après la synchronisation. Le résultat est donc le même.

Les rotations à partir des axes du matériau (`RotationMatrix2D`,
`RotationMatrix3D`) construisent les matrices dans la mémoire allouée par
la ressource de la politique. En asynchrone, elles marchent donc aussi
pour une bibliothèque qui calcule sur le GPU, ce que la version synchrone
refuse (matrices construites dans la mémoire de l'hôte).

Chaque appel asynchrone alloue son statut avec la ressource de la
politique. Avec une ressource de mémoire managée, c'est un
`cudaMallocManaged` et un `cudaFree` par appel. Le remède est standard :
un `std::pmr::synchronized_pool_resource` au-dessus de la ressource
managée, que la doc de MGIS recommande et que les tests utilisent. Les
ressources de MGIS acceptent pour cela tout alignement (voir
[Corrections de défauts existants](#corrections-de-défauts-existants)).
Avec un pool, détruire un calcul en attente avant la fin du calcul
corrompt l'appel suivant au lieu d'être seulement lent.

## Vues et rangements

Un `MaterialDataManager` range ses tableaux point par point par défaut.
MFEM et les codes GPU possèdent souvent leurs tableaux, rangés composante
par composante pour la coalescence. Le prototype ajoute à MGIS une API sur
des vues, qui décrit ces tableaux sans les copier. TFEL n'a pas eu à
changer : son contrat décrit déjà chaque champ par deux pas.

- **`FieldView`** : vue sur les valeurs d'un groupe de variables, avec un
  pas entre deux points et un pas entre deux composantes. Un pas de point
  nul décrit des valeurs uniformes. Fabriques : `makeFieldView`,
  `makePointByPointFieldView`, `makeComponentByComponentFieldView`,
  `makeUniformFieldView`, `getComponents`. Une sortie doit avoir des
  valeurs distinctes : MGIS le vérifie.
- **Rangement d'un gestionnaire** : `MaterialDataManagerInitializer::layout`
  (`POINT_BY_POINT` par défaut, ou `COMPONENT_BY_COMPONENT`). Tout ce qui
  est donné par point au gestionnaire suit ce rangement : états, blocs de
  la tangente, vitesse du son, variables externes à plusieurs composantes,
  entrées des fonctions d'initialisation, sorties des post-traitements,
  de l'extraction et des conversions des grandes transformations. `save`
  et `restore` écrivent toujours point par point.
- **`MaterialDataView`** : les données d'un matériau dans les tableaux de
  l'appelant. Toutes les fonctions du gestionnaire existent sur une vue.
  Celles du gestionnaire appellent celles de la vue : un seul chemin de
  code. Sans politique, elles prennent un espace de travail
  (`BehaviourIntegrationWorkSpace`), un par fil. Pas de `ThreadPool` sur
  une vue : le parallélisme CPU passe par les backends `openmp` et
  `stlpar`.
- **Rotations et propriétés matériau** : surcharges sur vues. Les
  surcharges sur tableaux passent par elles.
- **Fin de pas sur GPU** : `update`, `revert` et `updateValues` prennent
  des opérations de copie (`FieldOperations`), comme les allocations
  passent par une `std::pmr::memory_resource`. MGIS sait quels tableaux
  copier, l'appelant fournit la copie. En-têtes pour CUDA, HIP, SYCL,
  Kokkos et OpenMP, compilés par l'appelant avec sa chaîne GPU.
- **Vue sans contrôle et mise à jour des seules entrées** :
  `make_unchecked_view` rend la vue d'un gestionnaire comme `make_view`,
  sans vérifier que les valeurs gardées sur l'hôte (valeurs uniformes,
  stockage local, masse volumique) sont accessibles par les calculs.
  L'appelant remplace ensuite les vues des tableaux par celles de leurs
  copies sur le GPU, et copie ces valeurs dans une mémoire accessible.
  `updateInputValues` met à jour les propriétés matériau, les variables
  externes et la masse volumique d'un état à partir d'un autre, avec la
  même procédure et les mêmes politiques que `updateValues`, sans copier
  les tableaux de l'état (gradients, forces thermodynamiques, variables
  internes, énergies) : l'appelant les copie lui-même, par exemple sur
  leurs copies sur le GPU. `updateValues` copie maintenant les tableaux,
  puis appelle `updateInputValues`. mfem-mgis se sert des deux pour les
  données des matériaux qu'il garde sur le GPU.
- **`MaterialFunctionManager`** : constructeur avec initialiseur. Rangé
  composante par composante, il rend des vues *strided coalesced* de
  `mgis::function`.
- **Liaisons** : le rangement est exposé en C, Fortran et Python (tableaux
  numpy avec les pas du rangement).
- **Coût** : les noyaux de TFEL sont inchangés, seul le calcul des adresses
  change. Une fois les tableaux alloués, les fonctions synchrones qui
  calculent sur des points (intégration, initialisation, post-traitement,
  rotations, `evaluate`), `update` et `revert` n'allouent plus rien sur le
  tas dans le thread appelant, sur un gestionnaire comme sur une vue. Un
  test le vérifie en comptant les appels à `operator new`. Les threads du
  backend `stlpar` (TBB) allouent les leurs.
- **Points en échec** : les blocs de la tangente d'un point où
  l'intégration échoue sont non spécifiés, quel que soit le chemin.
- **Masse volumique** : exigée pour la vitesse du son d'une loi
  `generic-parallel` (TFEL l'exige), nulle si absente point par point,
  comme MGIS l'a toujours fait.

Le rangement par blocs d'élément (`byNODES` des interpolateurs de MFEM)
n'est pas couvert : il demanderait un troisième pas et un contrat TFEL
étendu. Les noyaux de MFEM lisent le rangement `byVDIM`, que les vues
couvrent.

La conception, ses variantes écartées et les relectures sont dans
[historique/conception-vues-mgis-2026-10-08.md](historique/conception-vues-mgis-2026-10-08.md).

## Erreurs

MGIS distingue deux cas, comme TFEL :

- **un échec du calcul en un point** (non convergence, borne physique) est
  un résultat : statut `-1` dans le résultat, avec l'indice du premier
  point en échec ;
- **un argument invalide ou une erreur du backend** lève une exception.
  Ce n'est jamais pris pour un échec du calcul.

Les fonctions d'un point (`integrate` sur une vue d'un point,
`integrate_debug`, fonctions d'initialisation et post-traitements d'un
point, rotations d'un point, liaisons C) refusent une bibliothèque de lots
avec un message. Elles appelaient un pointeur nul. Le diagnostic détaillé
d'un point se fait avec la bibliothèque `generic` de la même loi.

## Corrections de défauts existants

Ces défauts ont été trouvés en chemin et corrigés dans des commits séparés :

- les axes des matrices de rotation étaient refusés quand leur nombre de
  valeurs était correct (test inversé), et le nombre de points était
  comparé à un reste au lieu d'un quotient : aucune `RotationMatrix2D` ni
  `RotationMatrix3D` ne pouvait être construite ;
- une propriété matériau ne pouvait pas être donnée par un `std::span`
  d'une seule valeur, alors que l'évaluation la traitait déjà comme
  uniforme ;
- l'échec d'une fonction d'initialisation ou d'un post-traitement point
  par point ne mettait pas le statut à `-1` ;
- une valeur uniforme de propriété matériau en stockage local est rangée
  comme un scalaire, comme pour les variables externes ;
- la masse volumique peut être donnée par un `std::span` d'une seule
  valeur ;
- la restauration HDF5 d'une variable externe à plusieurs composantes en
  stockage local écrivait hors de son tableau (taille d'un point par
  composante au lieu de toutes) ;
- les ressources mémoire CUDA, HIP et Kokkos refusaient un alignement plus
  strict que celui de leur mémoire. Un pool de la bibliothèque standard
  demande des blocs alignés sur leur taille : il échouait au-delà de
  256 octets. L'alignement est maintenant obtenu par sur-allocation ;
- l'initialiseur `MaterialDataManagerInitializer` n'avait pas de
  constructeur en Python ;
- le message d'erreur des rotations des blocs de la tangente parlait des
  forces thermodynamiques.

- `MaterialAxisStorage` n'était pas exporté : copier un `RotationMatrix2D`
  ou un `RotationMatrix3D` hors de MGIS ne se liait pas. La copie d'un axe
  stocké localement pointait de plus vers le stockage de l'original ;
- les liaisons Fortran de l'opérateur tangent rendaient un tableau de
  forme fausse dès que l'opérateur n'était pas carré (et, pour un
  gestionnaire, de première étendue jamais fixée). Elles rendent la
  transposée de l'opérateur, de forme `[gs, ths]` en chaque point, et
  refusent un opérateur qui n'est pas un seul bloc ;
- la liaison Python lisait un tableau numpy à pas (`a[::2]`) comme s'il
  était contigu : il est maintenant refusé ;
- `update` et `revert` d'un gestionnaire allouaient un contexte d'exécution
  à chaque appel.

## Limites actuelles

- **Rangement par blocs d'élément** non couvert (voir
  [Vues et rangements](#vues-et-rangements)).
- **Pas de variante avec `Context`** : les fonctions de lot lèvent des
  exceptions, comme les autres fonctions sur un gestionnaire. mfem-mgis
  les enrobe dans des blocs `try`/`catch`.
- **`MaterialFunctionManager`** construit son gestionnaire sur l'hôte. Il
  n'est pas utilisé par mfem-mgis.
- **Intégration continue** : les tests `generic-parallel` de MGIS ne sont
  construits que si le TFEL trouvé fournit l'interface. La CI officielle
  de MGIS utilise TFEL master : ils n'y tourneront qu'une fois l'interface
  intégrée à TFEL. Les tests CUDA sont désactivés par défaut
  (`enable-generic-parallel-cuda-tests`). Les autres backends sont testés
  dans TFEL : MGIS ne fait que passer des pointeurs, une file et un espace
  mémoire.
- **Pas de statut par point** : le tableau optionnel de TFEL n'est pas
  exposé.
- **Liaisons** : les liaisons C refusent les bibliothèques de lots ; les
  liaisons Python et Fortran ne connaissent pas les fonctions de lot. Les
  trois exposent le rangement des gestionnaires, pas les vues.

## Exemple

Une loi compilée pour CUDA, des données en mémoire managée allouées par
MGIS, un flux CUDA choisi par l'appelant :

```cpp
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MGIS/Behaviour/MaterialDataManager.hxx"
#include "MGIS/Behaviour/Integrate.hxx"
#include "MGIS/MemoryResources/CUDAManagedMemoryResource.hxx"

using namespace mgis::behaviour;
auto resource = mgis::memory_resources::CUDAManagedMemoryResource{};
const auto b = load("libBehaviour-cuda.so", "Norton",
                    Hypothesis::TRIDIMENSIONAL);
auto i = MaterialDataManagerInitializer{};
i.location = MemoryLocation::SHARED;
i.resource = &resource;
auto m = MaterialDataManager(b, n, i);
// température uniforme, dans une mémoire lisible par le GPU
setExternalStateVariable(m.s0, "Temperature", temperature,
                         MaterialStateManager::EXTERNAL_STORAGE);
setExternalStateVariable(m.s1, "Temperature", temperature,
                         MaterialStateManager::EXTERNAL_STORAGE);
auto policy = BatchExecutionPolicy{};
policy.queue = stream;
const auto r = integrate(policy, m, BehaviourIntegrationOptions{}, dt);
```

## Documentation dans l'arbre de MGIS

La documentation de MGIS, en anglais, est à jour avec ces changements :

- `docs/web/generic-parallel.md` : page dédiée aux lots, reliée depuis
  l'accueil, le menu, l'introduction à la bibliothèque C++, la page des
  lois orthotropes et celle de l'analyse des échecs ;
- `docs/web/views.md` : page dédiée aux vues, aux rangements et aux
  opérations de fin de pas, reliée depuis l'accueil, le menu et la page
  des lots ;
- `docs/web/release-notes-3.2.md` : incompatibilités, nouvelle option de
  `cmake`, nouveautés et corrections ;
- `README.md`, `INSTALL-cmake.md`, `docs/web/install.md` : présentation,
  option `enable-generic-parallel-cuda-tests`, conditions des tests ;
- en-têtes publics : chaque fonction dit ce qu'elle fait pour une
  bibliothèque `generic-parallel` (bascule, refus, emplacement des
  tableaux). Les copies des structures de TFEL sont documentées comme
  leurs originaux.

La relecture a aussi corrigé des erreurs anciennes de la documentation :
`update` décrit à l'envers, noms des rotations des lois en grandes
déformations, paramètres mal annotés.
