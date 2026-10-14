# Plan : choix de l'exécution sur CPU, sur GPU ou mixte

Statut : plan du 13/10/2026 (tâche 9a). Schéma d'options retenu le
13/10/2026. Réalisé le 13/10/2026 (étapes 1 à 5, voir la dernière
section) : mfem-mgis 75bfd92a54 à 1348a95a5e, exemples dbf659f et suivant,
mm-opera-hpc 7416eb8 et suivants, deux relectures indépendantes traitées.

## Objectif

Avant cette tâche, placer un calcul sur le GPU demandait plusieurs
réglages indépendants, et chaque exemple les exposait à sa façon :

- l'option `--device` (ou la variable `MFEM_MGIS_DEVICE`) configure le
  périphérique de MFEM ;
- les variables `MFEM_MGIS_BEHAVIOUR_INTEGRATION`,
  `MFEM_MGIS_GRADIENTS_COMPUTATION`, `MFEM_MGIS_INNER_FORCES_COMPUTATION`,
  `MFEM_MGIS_TANGENT_OPERATOR_COMPUTATION` et
  `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION` choisissent les calculs par
  noyaux ;
- `--hypre-on-host` (ou `MFEM_MGIS_HYPRE_ON_HOST`) garde hypre sur l'hôte ;
- la bibliothèque de lois (`-l`) décide où les lois sont intégrées : une
  bibliothèque `generic-parallel` d'un backend GPU les intègre sur le GPU,
  une bibliothèque `generic` sur l'hôte ; les exemples ne construisent que
  la seconde ;
- ex1 et ex3 ne lisent que les variables d'environnement, les variantes
  MUMPS et PETSc d'ex4 n'ont pas d'option `-l`, le test de Robin n'a
  aucune option.

Le but : un seul jeu d'options, déclaré par mfem-mgis et compris de la
même façon par tous les exemples, et toutes les combinaisons de placement
possibles, sauf impossibilité technique réelle, documentée.

## Options

Déclarées par `mfem_mgis::declareDefaultOptions`, lues par
`mfem_mgis::initialize` :

- `--execution cpu|gpu|mixed` (`-x`) :
  - `cpu` : aucun périphérique, tout est calculé sur l'hôte ;
  - `gpu` : périphérique `gpu` de MFEM (CUDA ou HIP selon sa
    construction), tout est calculé sur le GPU ;
  - `mixed` : périphérique `gpu`, et chaque phase est placée par son
    sous-interrupteur (GPU par défaut) ;
- en `mixed` seulement (refusés sinon, car contradictoires) :
  - `--finite-elements cpu|gpu` : données de MFEM, gradients, forces
    intérieures, tangente ;
  - `--behaviours cpu|gpu` : intégration des lois (aussi matériau par
    matériau, par l'API) ;
  - `--linear-solver-device cpu|gpu` : hypre et les solveurs linéaires ;
  - `--post-processings cpu|gpu` ;
- `--device` est gardé pour les experts (noms de MFEM : `debug`, `omp`,
  `ceed-cuda`...). Il est refusé s'il contredit `--execution` (par
  exemple `--execution cpu --device cuda`).

Sans `--execution`, rien ne change : le comportement actuel (options et
variables existantes) est gardé, ce qui laisse les tests de mfem-mgis
inchangés.

Ordre de priorité, du plus fort au plus faible : paramètre donné au
problème, variable d'environnement, réglage déduit de `--execution`,
valeur par défaut. Un paramètre donné explicitement garde donc la main
(un test peut toujours forcer un chemin).

Réalisé autrement (voir [Réalisation](#réalisation-13102026)) : avec
`--execution`, les variables d'environnement de placement
(`MFEM_MGIS_DEVICE`, `MFEM_MGIS_HYPRE_ON_HOST`,
`MFEM_MGIS_BEHAVIOUR_INTEGRATION`, `MFEM_MGIS_*_COMPUTATION`) sont
ignorées, avec un avertissement. L'ordre est donc : paramètre donné au
problème, réglage déduit de `--execution`, valeur par défaut. Sans
`--execution`, les variables s'appliquent comme avant.

## Correspondance avec les réglages existants

| Phase | GPU | CPU |
|---|---|---|
| éléments finis | `BehaviourIntegration=Batch`, gradients, forces intérieures et tangente `Kernels` sur le périphérique | `ElementByElement`, boucles de MFEM sur l'hôte (réalisé : en `mixed` seulement ; `Kernels` sur l'hôte en mode `cpu`) |
| lois | bibliothèque d'un backend GPU de `generic-parallel` | bibliothèque `generic` (ou `generic-parallel` sur l'hôte) |
| solveur linéaire | hypre sur le GPU | hypre sur l'hôte (`--hypre-on-host`) |
| post-traitements | `PostProcessingsComputation=Kernels` sur le périphérique | `ElementByElement` sur l'hôte |

En mode `cpu`, les éléments finis et les post-traitements gardent le
chemin par défaut de mfem-mgis (élément par élément) ; les noyaux restent
accessibles par les paramètres et les variables existants, pour les
mesures.

Réalisé autrement : en mode `cpu`, les éléments finis sont calculés par
noyaux sur l'hôte, plus rapides, avec repli élément par élément si les
intégrateurs ne les supportent pas. Ils restent élément par élément quand
ils sont placés sur l'hôte en mode `mixed`. Seuls les post-traitements
gardent le chemin élément par élément sur l'hôte.

## Combinaisons

En mode `mixed`, les 16 combinaisons des quatre phases sont visées. Ce
qui est déjà possible :

- éléments finis sur l'hôte et lois sur le GPU : les gradients calculés
  élément par élément sont écrits dans les tableaux des matériaux, dans la
  mémoire de MFEM, puis lus sur le GPU par l'intégration par lots ;
- lois sur l'hôte et éléments finis sur le GPU : les tableaux sont
  synchronisés sur l'hôte avant l'intégration point par point, puis
  relus sur le GPU ;
- solveur sur le GPU avec une matrice assemblée sur l'hôte : MFEM copie la
  matrice dans la mémoire de hypre ;
- post-traitements sur l'hôte avec un calcul sur le GPU : synchronisation
  des matériaux sur l'hôte, déjà faite par les post-traitements.

Limites réelles, à refuser avec un message clair :

- un solveur qui n'existe que d'un côté : MUMPS, UMFPACK, PETSc (sur
  l'hôte) ; AmgX et cuDSS (sur le GPU). Le solveur choisi par
  l'application l'emporte : un conflit avec `--linear-solver-device` est
  une erreur, pas un déplacement silencieux ;
- le placement de hypre est global au processus (fixé avant son
  initialisation) : tous les solveurs de hypre d'un calcul sont du même
  côté ;
- le périphérique de MFEM est global : en `mixed`, il est toujours
  configuré, et « éléments finis sur le CPU » veut dire boucles de MFEM
  sur l'hôte, les vecteurs restant dans la mémoire de MFEM.

## Transferts

Les transferts entre l'hôte et le GPU passent par la mémoire de MFEM
(`Read`, `Write`, `HostRead`...) et par les synchronisations des
matériaux. Pour que leur surcoût soit mesurable :

- une section de profilage par synchronisation des matériaux
  (`Material::synchronizeOnHost`, `...OnDevice`), avec le nombre d'octets ;
- une section autour des copies de la matrice vers hypre quand
  l'assemblage et le solveur ne sont pas du même côté ;
- un tableau de mesure (poutre r4, une combinaison par ligne : temps total,
  temps des transferts) dans `results/`.

Copies asynchrones : la mémoire de MFEM copie de façon synchrone. Une
copie asynchrone (flux CUDA ou HIP, mémoire épinglée) ne sera ajoutée que
si la mesure montre un gain, et seulement là où l'ordre des opérations est
garanti (par exemple la recopie des contraintes vers l'hôte pendant
l'assemblage de la tangente sur le GPU).

## Lois

- Une fonction CMake commune, installée par mfem-mgis :
  `mfem_mgis_add_behaviours_library(NAME <nom> SOURCES <.mfront...>
  [BACKENDS ...])`, qui construit la bibliothèque `generic` et, si MFEM
  utilise un GPU, celle du backend GPU de `generic-parallel` (`cuda` ou
  `hip`), avec un nommage fixe. La signature réalisée est différente
  (voir l'étape 3 de la [Réalisation](#réalisation-13102026)).
- Une fonction de mfem-mgis qui choisit la bibliothèque d'après le
  placement des lois : `mfem_mgis::getBehavioursLibrary(ctx, nom)`.
- L'option `-l` des exemples force une bibliothèque (pour les mesures).

## Étapes

1. mfem-mgis : options et réglages d'exécution (`ExecutionSettings`),
   valeurs par défaut des paramètres déduites du placement, refus des
   contradictions ; test de toutes les combinaisons sur un petit problème
   (mêmes résultats que tout sur l'hôte) et des refus.
2. Sections de profilage des transferts.
3. Fonction CMake des lois et choix automatique de la bibliothèque.
4. Exemples de mfem-mgis-examples, puis mm-opera-hpc, puis ex1 et ex3 :
   options communes, lois construites par la fonction commune, variables
   d'environnement retirées des scripts.
5. Mesure du surcoût de chaque combinaison (poutre), documentation.

## Réalisation (13/10/2026)

- Étape 1 : options et réglages (`ExecutionSettings`,
  `getExecutionSettings`), lus avant l'initialisation de MPI et de hypre ;
  avec `--execution`, les variables d'environnement de placement sont
  ignorées (avertissement), un paramètre donné au problème garde la main.
  En mode `cpu`, les éléments finis sont calculés par noyaux sur l'hôte
  (plus rapide : 12,2 s au lieu de 16,5 s sur la poutre, 8 processus),
  avec repli élément par élément si les intégrateurs ne les supportent
  pas ; l'intégration par lots par défaut a le même repli. Le paramètre
  `BehaviourExecutionLocation` de `addBehaviourIntegrator` place les lois
  d'un matériau.
- Étape 2 : copies des tableaux des matériaux comptées
  (`getMaterialArraysTransfers`) et affichées après la table de profilage.
  Les autres transferts (propriétés matériau, variables externes, vecteurs
  et matrices de MFEM) ne sont pas comptés : le surcoût total se lit dans
  les temps.
- Étape 3 : `mfem_mgis_add_behaviours_library` (module CMake installé) et
  `getBehavioursLibrary`. Signature réalisée, sans `NAME` ni `BACKENDS` :
  `mfem_mgis_add_behaviours_library(<cible> SOURCES <fichiers.mfront>...
  [LIBRARIES <noms>...] [INCLUDE_DIRECTORY <dossier>] [DEPENDS
  <fichiers>...] [EXCLUDE_FROM_ALL] [NO_GPU])`. Les bibliothèques sont
  rangées dans `<dossier>/cpu/` et `<dossier>/gpu/`
  (`cmake/modules/MFEMMGISBehaviours.cmake` de mfem-mgis).
- Étape 4 : tous les exemples, mm-opera-hpc, ex1 et ex3 ; tests sur GPU par
  `--execution`.
- Étape 5 : `results/execution-beam-2026-10-13-1432`. Les copies des
  tableaux coûtent au plus 0,42 s (de 0,2 % à 4,5 % du calcul selon le
  placement) ; le coût d'un placement mixte vient du processeur qui
  calcule. Les copies asynchrones ne sont donc pas faites : le gain
  possible est petit devant le risque (ordre des opérations, mémoire
  épinglée).
- Limites réelles (documentées dans le guide de mfem-mgis) : solveurs qui
  n'existent que d'un côté ; solveurs de MFEM sur le GPU dès qu'un GPU est
  utilisé ; placement de hypre global ; cuDSS sur le GPU avec hypre sur
  l'hôte seulement par les options expertes. Avec HIP, les noyaux de
  mfem-mgis sont compilés pour HIP depuis la tâche 9c : compilation et
  édition de liens vérifiées pour un GPU AMD (gfx90a), pas d'essai sur un
  tel GPU.
- Pistes : garder les états sur le GPU après la mise à jour de fin de pas
  (faite sur le GPU, mais suivie d'un retour sur l'hôte voulu par le
  contrat de `MaterialArrays`, voir `results/L7-transfers-2026-10-11`),
  avec des écritures sur l'hôte synchronisées explicitement (120 copies
  évitables sur la poutre, 0,07 s) ; compter les autres transferts.
