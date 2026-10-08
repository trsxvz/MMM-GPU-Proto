# MGIS : détails techniques

Ce document décrit ce que le prototype change dans MGIS, et pourquoi. MGIS
charge les bibliothèques générées par MFront et gère les données matériau
(`MaterialDataManager`). Le prototype lui apprend à utiliser les
bibliothèques de l'interface `generic-parallel` : lois, modèles et
propriétés matériau calculés par lots, sur CPU ou sur GPU.

- Base : MGIS master `818f0f9` (02/10/2026).
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
9. [Erreurs](#erreurs)
10. [Corrections de défauts existants](#corrections-de-défauts-existants)
11. [Limites actuelles](#limites-actuelles)
12. [Exemple](#exemple)

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
- L'appelant déclare où sont ses tableaux. MGIS vérifie qu'ils sont
  accessibles par les calculs de la bibliothèque.
- MGIS ne dépend d'aucune bibliothèque GPU. Les tableaux viennent de
  l'appelant, ou d'une ressource mémoire qu'il fournit.

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
rotateGradients(policy, dest, behaviour, src, r);  // et les autres rotations
```

Les versions sans politique basculent seules. Pour les rotations, les
matrices peuvent être données par point ou uniformes. Les formes
`RotationMatrix2D` et `RotationMatrix3D` construisent les matrices sur
l'hôte à chaque appel : pour des rotations répétées sur GPU, mieux vaut les
construire une fois et les donner en tableau.

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
premier point en échec, message. Le calcul point par point s'arrête au
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
initialiseur (`MaterialDataManagerInitializer::location`). Pour une
propriété matériau, il est un argument d'`evaluate`.

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
`revert`, `updateValues`, l'extraction des variables internes, `save`,
`restore` et le calcul point par point. Elles refusent les tableaux placés
sur un GPU au lieu de planter.

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
  valeur.

## Limites actuelles

- **Mise à jour des états sur l'hôte.** `update` et `revert` copient les
  états sur l'hôte. Leur version sur GPU viendra avec mfem-mgis, qui
  possède la mémoire.
- **Pas de statut par point** : le tableau optionnel de TFEL n'est pas
  exposé.
- **Liaisons** : les liaisons C refusent les bibliothèques de lots ; les
  liaisons Python et Fortran ne connaissent pas encore les nouvelles
  fonctions.

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
