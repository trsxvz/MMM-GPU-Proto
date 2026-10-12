# TFEL et MFront : détails techniques

Ce document décrit ce que le prototype change dans TFEL et MFront, et
pourquoi. Il s'adresse à qui veut comprendre ou relire ces changements.
L'utilisation des interfaces, avec leurs exemples complets, est documentée
dans l'arbre de TFEL (voir [la dernière section](#documentation-dans-larbre-de-tfel)).

- Base : TFEL master `fec0ade99` (01/10/2026).
- Sources modifiées : `tfel/` (copie du dernier commit, sans l'historique de
  TFEL). Série de patchs : `patches/tfel/` (82 patchs, `git am`).
- Découpage prévu en contributions amont : [integration-amont.md](integration-amont.md).

## Sommaire

1. [Organisation des changements](#organisation-des-changements)
2. [Interface generic-parallel des lois et des modèles](#interface-generic-parallel-des-lois-et-des-modèles)
3. [Interface generic-parallel des propriétés matériau](#interface-generic-parallel-des-propriétés-matériau)
4. [Comparaison entre lois et propriétés matériau](#comparaison-entre-lois-et-propriétés-matériau)
5. [Code GPU dans TFEL et MFront](#code-gpu-dans-tfel-et-mfront)
6. [Corrections du code existant](#corrections-du-code-existant)
7. [Incompatibilités](#incompatibilités)
8. [Configurations et compilateurs](#configurations-et-compilateurs)
9. [Documentation dans l'arbre de TFEL](#documentation-dans-larbre-de-tfel)

## Organisation des changements

Les changements forment quatre familles.

| Famille | Contenu | Effet sur les résultats CPU |
|---|---|---|
| Code GPU dans TFEL | en-têtes mathématiques et matériau utilisables dans du code GPU | aucun |
| Code généré compatible GPU | classes générées par MFront, pilotes de l'interface `generic` partagés avec la nouvelle interface | aucun |
| Interfaces `generic-parallel` | nouvelle interface des lois et des modèles, extension de celle des propriétés matériau, tests | sans objet (nouvelles bibliothèques) |
| Corrections | défauts trouvés en chemin dans le code existant | corrigés, un commit chacun |

## Interface generic-parallel des lois et des modèles

MFront avait une interface `generic` qui intègre une loi en un point.
L'interface `generic-parallel` génère, pour la même loi, des fonctions qui
traitent tout un lot de points. Elle s'inspire de l'interface
`generic-parallel` des propriétés matériau, déjà présente dans TFEL.

### Fonctions générées et symboles

Pour chaque loi et chaque hypothèse de modélisation, la bibliothèque
exporte :

- la fonction d'intégration d'un lot, `<loi>_<Hypothèse>` ;
- une fonction par fonction d'initialisation,
  `<loi>_<Hypothèse>_InitializeFunction_<nom>` ;
- une fonction par post-traitement, `<loi>_<Hypothèse>_PostProcessing_<nom>` ;
- pour une loi orthotrope, les rotations des gradients, des forces
  thermodynamiques et des blocs de l'opérateur tangent (suffixées par la
  mesure de contrainte ou l'opérateur tangent en grandes déformations).

Elle exporte aussi des symboles de description, lus au chargement :

| Symbole | Valeur | Rôle |
|---|---|---|
| `<loi>_mfront_interface` | `GenericParallel` | distinguer une bibliothèque de lots d'une bibliothèque `generic` |
| `<loi>_generic_parallel_backend` | `cuda`, `hip`, `stlpar`, `kokkos`, `sycl` ou `openmp` | savoir où et comment la loi calcule |
| `<loi>_generic_parallel_api_version` | `MFRONT_GPB_API_VERSION` (2) | refuser une bibliothèque d'un contrat incompatible |
| `<loi>_generic_parallel_memory_space` | `host`, `device` ou `any` | savoir où placer les données (voir [Espace mémoire](#espace-mémoire)) |

Ces symboles permettent à l'appelant de tout vérifier au chargement, sans
compiler contre TFEL. C'est ce que fait MGIS.

### Contrat de données

Le contrat est en C, dans `MFront/GenericParallelBehaviour/BatchData.h`. Un
lot est décrit par une structure (`mfront_gpb_BatchData`) plutôt que par une
longue liste d'arguments. Une loi a une douzaine de groupes de variables, et
une structure peut évoluer sans changer le prototype des fonctions. Sa
version est exportée par les bibliothèques.

Chaque groupe de variables est un **champ** : un pointeur et deux pas.

```c
typedef struct mfront_gpb_ConstField {
  const mfront_gb_real* values;          // première valeur du premier point
  mfront_gb_size_type point_stride;      // entre deux points, 0 : uniforme
  mfront_gb_size_type component_stride;  // entre deux composantes
} mfront_gpb_ConstField;
```

Deux pas couvrent les deux rangements utiles :

- point par point, le rangement de MGIS : `point_stride` est la taille du
  groupe, `component_stride` vaut 1 ;
- composante par composante, qui donne des accès mémoire groupés sur GPU :
  `point_stride` vaut 1, `component_stride` est le nombre de points.

Les propriétés matériau et les variables externes peuvent aussi être
données **variable par variable**. Chacune est alors uniforme ou donnée par
point, indépendamment des autres. C'est la forme naturelle des données de
MGIS. Les éléments d'un tableau de variables (`a[2]`) sont des variables
distinctes, comme dans les symboles exportés par MFront.

La structure du lot contient encore le nombre de points, l'incrément de
temps, les options de calcul (les trois premières valeurs de `K` de
l'interface `generic`, même sens), le facteur de pas initial, les sorties
optionnelles (opérateurs tangents, facteurs de pas et statuts par point,
vitesse du son), les états de début et de fin de pas, et les options
d'exécution.

### Exécution en un point

Le travail en un point est écrit une fois, dans
`MFront/GenericParallelBehaviour/Integrate.hxx`, en gabarits utilisables
sur l'hôte et sur le GPU :

1. les valeurs du point sont copiées dans des tableaux locaux, de tailles
   connues à la génération du code ;
2. la vue `generic` du point est construite sur ces copies, et la stratégie
   d'intégration de la loi est appelée ;
3. les sorties sont recopiées dans les champs, et le statut du point est
   enregistré.

Un accès direct aux champs, sans copie locale, a été mesuré : aucun gain à
l'exécution, et 25 à 40 % de temps de compilation en plus avec nvcc. Il a
été retiré.

Les données d'un point sont décrites par un type généré (`PointTraits`).
Il donne les tailles maximales des données d'un point et la stratégie
d'intégration de la loi. Un type plutôt qu'une valeur, car nvcc ne sait pas
générer le lanceur hôte d'un noyau dont un paramètre de gabarit est une
classe.

La classe de loi générée est la même pour les deux interfaces.
L'interface `generic-parallel` dérive de l'interface `generic` et réutilise
les membres qu'elle ajoute à la classe. Le générateur n'écrit qu'une fois
les blocs identiques produits par plusieurs interfaces : une même loi peut
donc être générée pour les deux interfaces à la fois. L'interface ajoute
seulement un constructeur qui reçoit les paramètres.

L'interface `generic-parallel` n'a pas de code d'intégration propre. Les
pilotes de l'interface `generic` (petites déformations, grandes
déformations, stratégies Green-Lagrange et Hencky, fonctions
d'initialisation, post-traitements) ont été coupés en deux :

- un cœur utilisable sur hôte et GPU, qui prend les options de calcul
  explicitement ;
- les fonctions historiques de l'hôte, de signature inchangée, qui
  décodent les options et appellent le cœur.

Les deux interfaces partagent donc les mêmes cœurs. Les stratégies
d'intégration (`StandardIntegrator`, `FiniteStrainIntegrator`,
`GreenLagrangeStrainIntegrator`, `LogarithmicStrainIntegrator`,
`ModelIntegrator`) sont des types choisis par le code généré.

Les fonctions d'initialisation et les post-traitements sont appelés par des
objets fonctions générés, jamais par des pointeurs vers membres : `nvc++`
ignore ces pointeurs sans message dans le code GPU.

### Backends

Chaque backend fournit une seule primitive : appliquer un objet fonction
aux points du lot et réduire les statistiques (nombre d'échecs, de résultats
peu fiables, plus petit facteur de pas, premier point en échec). Chacun
utilise la primitive idiomatique de son modèle de programmation.

| Backend | Primitive | File d'exécution | Asynchrone | Options |
|---|---|---|---|---|
| `cuda`, `hip` | un fil par point, vote par warp, opérations atomiques seulement en cas d'événement | flux (`cudaStream_t`, `hipStream_t`) | oui | `number_of_threads_per_block` (128 par défaut) |
| `stlpar` | `std::transform_reduce` sur un itérateur de comptage | aucune | non | `execution_policy`, `device_stack_size` |
| `kokkos` | `Kokkos::parallel_for`, opérations atomiques seulement en cas d'événement | instance de l'espace d'exécution | oui | `execution_space` |
| `sycl` | `parallel_for` et mises à jour atomiques | `sycl::queue`, obligatoire | oui, si aucune valeur n'est copiée et si la file est ordonnée | `data_transfer_policy`, `device_stack_size` |
| `openmp` | `target teams distribute parallel for` ou `parallel for` | aucune | non | `offload`, `device_stack_size` |

Quelques choix méritent une explication :

- SYCL n'utilise pas de réduction : leurs barrières font échouer
  AdaptiveCpp sur les gros noyaux exécutés sur CPU.
- Le noyau doit s'exécuter après les calculs déjà soumis à la file, qui
  peuvent produire ses entrées. Une file SYCL ordonnée le garantit. Avec
  une file non ordonnée, un calcul synchrone attend d'abord la fin de la
  file. Un calcul asynchrone y est refusé (arguments invalides) : rien ne
  garantirait l'ordre. La règle est la même pour les lois et les
  propriétés matériau.
- Le nombre de fils par bloc de CUDA et HIP est arrondi au multiple de la
  taille du warp, lue sur le GPU à l'exécution. Le vote par warp suppose
  des warps complets. La même règle vaut pour les lois et les propriétés
  matériau.
- Kokkos doit être une bibliothèque **partagée**, compilée **sans code
  relogeable** avec CUDA. Les bibliothèques générées sont chargées par
  `dlopen` : l'appelant et elles doivent partager la même instance de
  Kokkos. La fonction générée vérifie que Kokkos est initialisé, elle ne
  l'initialise jamais. mfem-mgis construit avec l'option `enable-kokkos`
  l'initialise à la demande (option `--kokkos` ou variable
  `MFEM_MGIS_KOKKOS`, voir `results/kokkos-polycristal-2026-10-11`). Les
  statistiques d'un calcul synchrone sont dans
  `Kokkos::SharedSpace`, comme pour les propriétés matériau. Les backends
  `OpenACC` et `OpenMPTarget` de Kokkos n'ont pas cet espace : une
  assertion statique le signale.
- Le backend `stlpar` sert aussi le GPU avec `nvc++ -stdpar=gpu`. Les
  bibliothèques liées à TBB (gcc, clang, icpx) restent chargées jusqu'à la
  fin du processus (`-Wl,-z,nodelete`) : les fils de TBB tournaient encore
  quand elles étaient déchargées à la sortie.

### Options d'exécution

Les options d'exécution (`mfront_gpb_ExecutionOptions`) ont le même sens
pour tous les backends :

- `queue` : file sur laquelle le calcul est soumis. Un pointeur nul désigne
  la file par défaut, quand elle existe.
- `statistics` : statistiques d'un calcul asynchrone, ou nul pour un calcul
  synchrone (voir [Synchrone et asynchrone](#synchrone-et-asynchrone)).
- `preserve_stack_size` : si non nul, la taille de pile des fils du GPU
  n'est jamais modifiée. Une erreur est rapportée si elle est trop petite.

Un lot vide est accepté sans aucune vérification : rien n'est lu, écrit
ou soumis, les champs et la file peuvent être nuls. La règle est la même
pour les propriétés matériau.

### Espace mémoire

Chaque bibliothèque déclare où ses calculs lisent les données. La valeur
est fixée à la compilation, d'après le backend, le compilateur et le mode
de mémoire (`MFront/GenericParallel/MemorySpace.hxx`).

| Backend et compilation | Espace mémoire |
|---|---|
| `cuda`, `hip` | `device` |
| `stlpar`, `openmp` sur l'hôte | `host` |
| `nvc++ -stdpar=gpu`, `nvc++ -mp=gpu` | `device`, ou `any` en mémoire unifiée |
| `kokkos` | `host` si l'espace d'exécution accède à la mémoire hôte, `device` sinon |
| `sycl` | `any` avec les transferts automatiques, `device` sinon |

`device` signifie que les calculs ne peuvent pas lire la mémoire ordinaire
de l'hôte. Les données doivent être sur le GPU ou en mémoire partagée.
`any` signifie que toute donnée convient. L'appelant n'a pas besoin de
connaître la liste des backends : ce symbole lui suffit.

### Erreurs et statuts

La fonction d'un lot renvoie un statut, aussi écrit dans la structure de
sortie (`mfront_gpb_OutputStatus`) :

| Statut | Sens |
|---|---|
| `1` | succès en tous les points |
| `0` | succès, mais le facteur de pas est inférieur à 0,99 en certains points |
| `-1` | échec du calcul en certains points |
| `-2` | arguments invalides, rien n'est calculé |
| `-3` | erreur du backend (allocation, lancement d'un noyau, etc.) |

En cas d'échec, la structure donne le nombre de points en échec et
**l'indice du premier**. Cet indice ne dépend pas de l'ordre des calculs
(minimum atomique ou réduction). Un tableau optionnel donne le statut de
chaque point.

Le message d'erreur détaillé d'un point n'est pas transmis : le code GPU
n'a ni exceptions ni tampon de message par fil. Pour l'obtenir, on rejoue le
point avec la bibliothèque `generic` de la même loi.

**Avec les opérations atomiques, le compte rendu des erreurs ne coûte rien
en cas de succès.** Les noyaux n'écrivent dans les statistiques que si un
événement survient en un point. Sur CUDA et HIP, un vote par warp évite
toute opération atomique dans les warps sans événement. En mode synchrone,
les backends `cuda`, `hip`, `sycl` et `kokkos` tiennent les statistiques
dans des blocs de mémoire accessibles par l'hôte et le GPU, alloués au
premier appel puis réutilisés
(`MFront/GenericParallel/StatusBlockPool.hxx`). Les backends `stlpar` et
`openmp` fusionnent les statuts des points par une réduction : son
résultat, quelques entiers, est copié vers l'hôte à chaque appel.
L'hôte initialise un bloc et le lit après la synchronisation : sans
événement, aucune donnée ne passe entre l'hôte et le GPU. Mesure, sur une
loi élastique et 1000 points : l'appel synchrone CUDA passe de 135 à 29 µs.

### Synchrone et asynchrone

Si les options d'exécution ne donnent pas de statistiques, le calcul est
synchrone. La fonction rend la main une fois le calcul fini et le statut
décrit le résultat.

Si elles en donnent, les statistiques sont cumulées dans cette structure,
qui doit être accessible par le GPU (et par l'hôte pour `stlpar` et
`openmp`, qui cumulent sur l'hôte). Les backends `cuda`, `hip`, `kokkos`
et `sycl` (sans copie) rendent alors la main dès le calcul soumis. L'appelant
synchronise la file, puis `mfront_gpb_report` (`BatchData.h`) remplit le
statut de sortie comme le fait un calcul synchrone, message compris.
C'est le pendant de `<nom>_generic_parallel_report` des propriétés
matériau. Plusieurs lots et l'assemblage peuvent ainsi être soumis avant
une seule synchronisation.

### Pile des fils du GPU

Certaines lois demandent une pile plus grande que la valeur par défaut des
GPU NVIDIA (1 Ko). Les backends `cuda` et `hip` calculent la taille
nécessaire à partir des attributs du noyau
(`MFront/GenericParallel/GPUStackSize.hxx`). Kokkos passe par un noyau
sonde jamais lancé (`MFront/GenericParallel/KokkosStackSize.hxx`). Ces
deux en-têtes servent aux lois et aux propriétés matériau. `nvc++` et
AdaptiveCpp ne donnent pas accès à cette taille : l'option
`device_stack_size` la fixe par le runtime CUDA chargé dans le processus
(`MFront/GenericParallel/DeviceStackSize.hxx`). La taille de pile est
commune à tous les calculs du processus : `preserve_stack_size` permet de
l'interdire.

Le noyau est un paramètre template de la fonction qui lit ses attributs.
Le runtime CUDA est lié statiquement par défaut : chaque bibliothèque a sa
propre instance, qui ne connaît que ses noyaux. Une fonction partagée qui
prendrait le noyau en argument serait résolue par l'éditeur de liens
dynamique vers l'instance d'une seule bibliothèque, et le noyau d'une
autre y serait inconnu (`invalid resource handle`). Le test CUDA des
propriétés matériau, qui charge trois bibliothèques, l'a montré.

### Couverture

L'interface couvre tout ce que l'interface `generic` sait faire :

- lois en petites déformations et lois générales (thermique, zones
  cohésives) ;
- grandes déformations, y compris les stratégies Green-Lagrange et Hencky
  et les conversions d'opérateur tangent ;
- fonctions d'initialisation et post-traitements ;
- rotations des lois orthotropes : gradients, forces thermodynamiques et
  tous les blocs de l'opérateur tangent, par un module partagé avec
  l'interface `generic` (`MFront/BehaviourRotationFunctions.hxx`) ;
- monocristaux (systèmes de glissement et matrices d'interaction construits
  à la compilation) ;
- variables de comportement et leurs paramètres, brique
  `StandardElastoViscoPlasticity`, solveur linéaire TDLS dans les lois
  implicites ;
- données partagées par tous les points (`@SharedData`, voir
  [Code GPU](#code-gpu-dans-tfel-et-mfront)) ;
- modèles autonomes : ils sont traités comme des lois sans gradient, comme
  par l'interface `generic`. Une bibliothèque de modèles exporte les mêmes
  fonctions et symboles qu'une bibliothèque de lois.

## Interface generic-parallel des propriétés matériau

Cette interface existait dans TFEL. Elle évalue une propriété matériau sur
un ensemble de points. Le prototype l'étend pour qu'elle offre les mêmes
possibilités que celle des lois.

### Ce qui a changé

- **Backends** : `kokkos` et `openmp` ajoutés aux backends existants
  (`cuda`, `hip`, `stlpar`, `sycl`).
- **Propriétés définies par des données ou appelant d'autres lois** :
  utilisables sur CUDA et HIP.
- **Options d'exécution en dernier argument** (version 2 des fonctions,
  `MFRONT_GPMP_API_VERSION`). Elles ont les membres des options des lois,
  avec le même sens : file, statistiques d'une évaluation asynchrone,
  préservation de la pile. Un pointeur nul donne les options par défaut.
  La file SYCL était un argument propre à ce backend, inséré au milieu des
  autres : la signature dépendait du backend. Elle passe maintenant par
  les options, et la signature est la même partout.
- **Statut de sortie et statistiques calqués sur les lois.** Le statut de
  sortie (`mfront_gpmp_OutputStatus`) donne, en plus du statut et du rang
  de la variable hors bornes, le nombre de points en échec, le premier
  d'entre eux et le nombre de points hors bornes avec la politique
  `Warning`. Chaque noyau calcule le statut de son point et ne le fusionne
  dans les statistiques (`mfront_gpmp_Statistics`) qu'en cas d'événement :
  opérations atomiques pour CUDA, HIP, SYCL et Kokkos, réduction pour les
  algorithmes parallèles et OpenMP, comme pour les lois. Coût nul en cas
  de succès avec les opérations atomiques.
- **Pile des fils** : CUDA, HIP et Kokkos l'adaptent au noyau, avec les
  en-têtes des lois. `preserve_stack_size` a donc le même effet partout.
- **Transferts SYCL** : la classe de copies des lois
  (`MFront/GenericParallel/SYCLDeviceCopies.hxx`) sert aussi aux
  propriétés matériau. Les statistiques d'une évaluation asynchrone sont
  copiées si besoin. La sortie est copiée à l'aller : un point en échec
  garde sa valeur. Une allocation USM de tout type est accessible, comme
  pour les lois. Les erreurs du noyau remontent (`wait_and_throw`).
- **Exceptions des runtimes** : rapportées par le statut `-6`, comme les
  erreurs des runtimes CUDA et HIP. Le statut `-2` n'est plus utilisé.
- **Ensemble vide** : accepté sans vérification, comme un lot vide.
  `errno` est restauré à chaque retour.
- **Nombre de points** : type `size_t` dans les noyaux CUDA, HIP et
  Kokkos et avec AdaptiveCpp (limite de 2^31 points levée).
- **Évaluation asynchrone.** Si les options donnent des statistiques,
  initialisées par `mfront_gpmp_initializeStatistics`, la fonction rend la
  main dès le calcul soumis (`cuda`, `hip`, `kokkos`, `sycl` sans copie).
  Une fois le calcul fini, la fonction exportée
  `<nom>_generic_parallel_report` analyse les statistiques et remplit le
  statut de sortie. L'évaluation synchrone appelle la même fonction :
  l'analyse n'est écrite qu'une fois.
- **Symboles** : backend, version des fonctions, espace mémoire, comme
  pour les lois.
- **Blocs réutilisés** : les statistiques d'une évaluation synchrone
  viennent d'une réserve de blocs, sans allocation par appel ; attente du
  flux au lieu de tout le GPU (CUDA, HIP) ; attente de l'instance utilisée
  au lieu de tout Kokkos.
- **Statut 1 avec la politique `Warning`** : une borne violée donne
  maintenant le statut 1 (résultat calculé, à utiliser avec précaution),
  comme l'interface `generic`. Le statut restait nul.

Le statut et le rang de la variable gardent le sens de l'interface
`generic` des propriétés matériau : `0` succès, `1` borne violée avec la
politique `Warning` (rang positif), `-1` borne physique violée ou borne
violée avec la politique `Strict` (rang négatif), `-4` résultat non fini.
L'interface ajoute `-5` (arguments invalides, file SYCL absente,
évaluation asynchrone sur une file SYCL non ordonnée), `-6` (erreur du
backend, y compris une exception de son runtime) et `-7` (sortie uniforme
incohérente). Quand
plusieurs vérifications échouent, la première dans l'ordre suivant est
rapportée : résultat non fini, bornes physiques, bornes, les variables
dans l'ordre, la sortie après les arguments.

## Comparaison entre lois et propriétés matériau

| | Lois et modèles | Propriétés matériau |
|---|---|---|
| Appel | structure du lot | arguments et pas |
| Options d'exécution | file, statistiques, pile | mêmes membres |
| Synchrone ou asynchrone | au choix de l'appelant | au choix de l'appelant |
| Compte rendu asynchrone | `mfront_gpb_report` (fonction C de `BatchData.h`) | `<nom>_generic_parallel_report` (symbole, il nomme la variable hors bornes) |
| Lot vide | accepté sans vérification | accepté sans vérification |
| Exception du runtime du backend | `-3` (erreur du backend) | `-6` (erreur du backend) |
| Transferts SYCL automatiques | `SYCLDeviceCopies.hxx` | même classe |
| Fusion des statuts des points | atomique (CUDA, HIP, SYCL, Kokkos), réduction (algorithmes parallèles, OpenMP) | identique |
| Nombre de points en échec, premier point en échec | oui | oui |
| Pile des fils | adaptée au noyau, `preserve_stack_size` | identique |
| Espace mémoire, backend, version | symboles exportés | symboles exportés |
| Tableau optionnel des statuts par point | oui | non |

Il reste une différence : le tableau optionnel des statuts par point. Il
n'est pas ajouté aux propriétés matériau (option écartée). MGIS ne
l'expose pour aucune des deux interfaces.

Les statistiques ne peuvent pas avoir exactement la même structure. Celles
des lois portent le facteur de pas de temps, celles des propriétés
matériau la première vérification en échec, qui sert à nommer la variable
dans le message.

## Code GPU dans TFEL et MFront

Le code exécuté sur GPU n'a ni exceptions, ni fonctions virtuelles appelées
sur des objets créés par l'hôte, ni accès aux variables globales de l'hôte.
Les changements suivants rendent TFEL et le code généré compatibles, sans
changer les résultats sur CPU.

- **Annotations.** Les en-têtes mathématiques et matériau sont annotés
  (`TFEL_HOST_DEVICE`). Chaque en-tête est compilé avec clang CUDA, qui
  refuse tout appel d'une fonction hôte depuis le code GPU
  (`scripts/clang-cuda-header-scan.sh`).
- **Code GPU et compilateurs à une passe.** `TFEL_DEVICE_CODE` est défini
  quand du code GPU est compilé. Les compilateurs à une passe (`nvc++`,
  AdaptiveCpp) le définissent dans toute l'unité de compilation : seules les
  bibliothèques `generic-parallel` qu'ils compilent suivent alors les règles
  du code GPU sur l'hôte.
- **Lois sans fonctions virtuelles.** `MechanicalBehaviour` ne déclare plus
  l'interface des lois en fonctions virtuelles pures. La table virtuelle
  demande des informations de type à l'exécution, absentes du code GPU. Le
  code généré ne l'utilisait pas.
- **Types trivialement copiables.** L'opérateur tangent des lois en grandes
  déformations utilise `TrivialGenTypeBase` au lieu de `std::variant`, dont
  la visite n'est pas supportée par `nvc++` dans le code GPU.
- **Données constantes.** Les données construites à la compilation
  (systèmes de glissement) passent par `tfel::utilities::getConstantObject`.
  nvcc remplit de zéros les variables statiques de fonction dans le code
  GPU.
- **Valeurs propres** sans récursion ni assertion dans le code GPU
  (algorithme de Jacobi pour les valeurs propres multiples).
- **Splines cubiques** : itérateurs avancés sans `std::advance`,
  `std::next` et `std::prev`, que nvcc compile en fonctions vides dans le
  code GPU.
- **Erreurs portables.** `tfel::fail<Exception>(message)` lève l'exception
  sur l'hôte et renvoie `false` sur le GPU. `abortResolution` (DSL
  implicites) et `abortIntegration` (Runge-Kutta) arrêtent l'intégration
  d'un point. Les 36 lois de test dont le code utilisateur n'était pas
  utilisable sur GPU ont été rendues portables avec ces outils, sans
  changer leurs résultats.
- **Données partagées** (`@SharedData`). Des données communes à tous les
  points, par exemple les orientations d'un polycristal lues dans un
  fichier, sont calculées sur l'hôte puis copiées une fois par GPU et par
  contexte.
- **Paramètres des variables de comportement.** Les valeurs des paramètres
  sont lues sur l'hôte et transmises par valeur aux noyaux
  (`mfront_parameters_values`).
- **Lois de matériau appelées par les lois** : une implémentation dans
  l'espace de noms `mfront::device`, générée dans l'en-tête de l'interface
  `mfront` des propriétés matériau.
- **Même classe de loi sur l'hôte et le GPU**, et comparaison à la
  jacobienne numérique faite seulement sur l'hôte.
- **Variables statiques.** MFront écrit une variable statique
  (`@StaticVariable`, `@IntegerConstant`, nombres de systèmes de
  glissement) `TFEL_STATIC_VARIABLE(v)`, dans les blocs de code et dans les
  entrées des propriétés matériau. Cette macro de `TFELConfig.hxx` donne la
  variable elle-même dans le code de l'hôte, qui ne change donc pas, et une
  copie initialisée par une expression constante, `(decltype(v){v})`, dans
  le code GPU (`TFEL_DEVICE_CODE`). Le code GPU ne peut pas faire
  référence à une variable statique constexpr de l'hôte, ce que fait un
  opérateur qui prend ses arguments par référence dans une expression
  évaluée à l'exécution, par exemple `T + dT - Tref` avec les quantités
  (`@UseQt true`, nvcc : « identifier is undefined in device code »).
  Les algorithmes de Runge-Kutta écrivent aussi ainsi la rigidité
  recalculée aux étapes intermédiaires et les constantes de rkCastem.
  Reproducteur : `reproducers/mfront-static-quantity-cuda` ; loi de test
  `StrainHardeningCreep3`.

## Corrections du code existant

Ces défauts ont été trouvés en chemin. Chacun est corrigé dans son propre
commit et couvert par un test existant ou nouveau. Certains changent les
résultats sur CPU : c'était le but.

Interface `generic` des lois :

- statut vérifié après une intégration en grandes déformations ;
- modules de l'opérateur `DTAU_DDF` avec la déformation de Hencky ;
- état de fin de pas donné aux post-traitements ;
- fonctions d'initialisation des stratégies Green-Lagrange et Hencky ;
- vitesse du son avec l'opérateur de prédiction ;
- pointeurs restaurés après une exception ;
- contrainte des fonctions d'initialisation pour toutes les mesures de
  contrainte ;
- exceptions du constructeur de la loi rapportées comme des échecs ;
- rotations : décalages, scalaires et vecteurs corrigés ; rotation de
  tous les blocs de l'opérateur tangent (tableaux de variables, couples de
  variables sans type TFEL).

Modèles et propriétés matériau :

- bornes des modèles `generic` vérifiées dans tous les modes de
  compilation ;
- classes générées pour les modèles `generic` en liaison interne : les
  paramètres d'un modèle étaient partagés par toutes ses bibliothèques
  chargées dans un processus ;
- bornes physiques des arguments des propriétés matériau exportées même
  sans bornes standard ;
- interface `generic` des propriétés matériau : message de la borne
  physique supérieure (« below » au lieu de « over »), `errno` restauré
  quand une borne supérieure est violée avec la politique `Strict`, et
  documentation du signe du rang de la variable (négatif quand
  l'évaluation échoue, y compris pour une borne avec la politique
  `Strict`, comme le code le fait et comme les tests de MGIS l'attendent) ;
- interface `generic-parallel` des propriétés matériau : code généré avec
  les vérifications désactivées, mémoire SYCL libérée après chaque appel,
  pointeurs des arguments copiés dans les noyaux (le tableau de l'appelant
  peut ne pas être lisible par le GPU), étendue des champs avec pas et
  valeurs entre les sorties avec pas préservées (SYCL), transferts dans le
  bloc `try`, backend dans la garde de l'en-tête, bornes et valeurs par
  défaut des paramètres écrites avec 14 chiffres, taille de grille des
  noyaux CUDA, `errno` remplacé par la détection des résultats non finis
  dans les noyaux (`errno` est par fil et absent sur GPU).

Divers :

- options de l'éditeur de liens après les fichiers objets dans les
  Makefile générés ;
- loi de test `DigitSeparatorTest` : dérivée de la loi d'écoulement non
  initialisée ;
- test CUDA des propriétés matériau qui écrivait en mémoire hôte depuis un
  noyau (masqué par la gestion hétérogène de la mémoire) ;
- algorithmes de Runge-Kutta : les modificateurs de la rigidité gardaient,
  dans une lambda statique, le nom de la première loi traitée par un appel
  de MFront ; avec plusieurs lois dans le même appel, les autres faisaient
  référence aux variables statiques de la première ;
- brique DDIF2 : test inversé sur le nombre de contraintes de rupture non
  constantes (une seule était refusée ; avec trois, seule la première
  était utilisée) et première énergie de rupture utilisée dans les trois
  directions. Vérifié par paires de lois (valeurs constantes et valeurs
  non constantes égales, MTest dans chaque direction) ;
- DSL isotropes avec `use_stress_update_algorithm` : coefficient de
  Poisson converti en `stress` au lieu de `real`, ce qui ne compilait pas
  avec `@UseQt true` (loi de test `Plasticity8`, mêmes résultats que
  `Plasticity6`) ;
- message de débogage de rkCastem (facteur affiché faux).

Défaut relevé mais pas corrigé : dans rkCastem, le facteur de
normalisation de l'erreur ne prend jamais `young` (le résultat de
`checkVariable("young")` est ignoré dans
`findStressErrorNormalizationFactor`). Le corriger changerait les pas de
temps des lois existantes : à décider avec les développeurs de TFEL.

## Incompatibilités

Le support du code GPU change quelques interfaces de programmation. La
liste est dans les notes de version de TFEL
(`tfel/docs/web/release-notes-5.2.md`, section « Known incompatibilities ») :
`MechanicalBehaviour` sans fonctions virtuelles, fonctions de
`BoundsCheck`, `LogarithmicStrainHandler` défini dans son en-tête, systèmes
de glissement construits à la compilation. Les fonctions des propriétés
matériau `generic-parallel` prennent un argument de plus (version 2).

## Configurations et compilateurs

TFEL est construit et testé avec 14 configurations de l'interface des lois
(`build-tfel.sh`). Les fichiers de configuration sont dans
`tfel/mfront/configuration-files/`.

| Configuration | Compilateur | Cible |
|---|---|---|
| `cuda-nvcc`, `cuda-clang` | nvcc 13.0, clang 20 | GPU NVIDIA |
| `hip-hipcc`, `hip-clang` | hipcc 7.2, clang 20 | GPU NVIDIA (plateforme NVIDIA de HIP) |
| `stlpar-parunseq-gcc`, `-clang`, `-icpx` | g++ 14, clang 20, icpx 2026.1 | CPU, TBB |
| `stlpar-parunseq-nvhpc-gpu` | nvc++ 25.11 | GPU NVIDIA |
| `kokkos` | Kokkos 5.2 | GPU NVIDIA et CPU (OpenMP) |
| `sycl-default-acpp`, `sycl-gpu-acpp` | AdaptiveCpp 25.10 | CPU, GPU NVIDIA |
| `sycl-default-icpx` | icpx 2026.1 | CPU |
| `openmp-host-gcc` | g++ 14 | CPU |
| `openmp-offload-nvhpc-gpu` | nvc++ 25.11 | GPU NVIDIA |

Options notables :

- `nvc++ -stdpar=gpu` utilise `-gpu=mem:managed,nordc`. Le mode mémoire
  unifiée par défaut est mal compilé par nvc++ 25.11. Sans `nordc`, des
  symboles d'initialisation statique de même nom entrent en conflit à
  l'édition de liens.
- Les configurations `stlpar` liées à TBB ajoutent `-Wl,-z,nodelete` (voir
  [Backends](#backends)).

Les défauts de compilateurs rencontrés sont contournés de façon ciblée,
jamais en retirant une fonctionnalité sur GPU. Chacun a un reproducteur :
`results/nvhpc-bugs/README.md` et `reproducers/README.md`.

## Documentation dans l'arbre de TFEL

La documentation de TFEL, en anglais, est à jour avec ces changements.
Chaque affirmation a été vérifiée dans le code.

- `tfel/docs/web/generic-parallel-behaviour-interface.md` : utilisation de
  l'interface des lois et des modèles (fonctions, contrat de données,
  états, statuts, files de chaque backend, fichiers de configuration
  livrés, écriture de lois compatibles GPU, tests). La page est
  construite et figure dans le menu du site.
- `tfel/docs/web/generic-parallel-material-property-interface.md` :
  utilisation de l'interface des propriétés matériau (prototypes, options
  d'exécution, évaluation asynchrone, statuts, backends).
- `tfel/docs/web/release-notes-5.2.md` : lois sur GPU, incompatibilités,
  nouveautés de `TFEL/Math`, `TFEL/Utilities` et `TFEL/Material`,
  corrections.
- Autres pages mises à jour : interface `generic` (rotations de tous les
  blocs, noms des fonctions des lois en grandes déformations), `TFEL/Config`
  (substitutions et configurations), `tfel-math` (`abortResolution`),
  variables de comportement, fonctions d'initialisation, post-traitements,
  tutoriels des polycristaux, aide des mots-clés (`@SharedData`,
  `@CompareToNumericalJacobian`, `@Profiling`). La section `@SharedData`
  est ajoutée aux pages des mots-clés des DSL de lois.
- Commentaires de code et documentation écrite dans le code généré :
  espace mémoire des lots, files SYCL, listes de backends, pas des
  tableaux, paramètres du constructeur.

Les pages de mots-clés ne sont pas régénérées en entier : TFEL a corrigé
des fautes directement dans ces pages, pas dans les fichiers d'aide dont
elles sont tirées. Une régénération complète les ferait revenir.

Un défaut de CMake a été trouvé pendant cette relecture : les tests
`sycl-gpu-acpp` utilisent les fichiers de configuration de
`sycl-default-acpp`, qui n'étaient installés que si cette configuration
était aussi choisie. Ils sont maintenant installés avec l'une ou l'autre.
