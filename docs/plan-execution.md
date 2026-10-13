# Plan : choix de l'exécution sur CPU, sur GPU ou mixte

Statut : plan du 13/10/2026 (tâche 9a). Schéma d'options retenu le
13/10/2026 ; ce document en tire la réalisation.

## Objectif

Aujourd'hui, placer un calcul sur le GPU demande plusieurs réglages
indépendants, et chaque exemple les expose à sa façon :

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

## Correspondance avec les réglages existants

| Phase | GPU | CPU |
|---|---|---|
| éléments finis | `BehaviourIntegration=Batch`, gradients, forces intérieures et tangente `Kernels` sur le périphérique | `ElementByElement`, boucles de MFEM sur l'hôte |
| lois | bibliothèque d'un backend GPU de `generic-parallel` | bibliothèque `generic` (ou `generic-parallel` sur l'hôte) |
| solveur linéaire | hypre sur le GPU | hypre sur l'hôte (`--hypre-on-host`) |
| post-traitements | `PostProcessingsComputation=Kernels` sur le périphérique | `ElementByElement` sur l'hôte |

En mode `cpu`, les éléments finis et les post-traitements gardent le
chemin par défaut de mfem-mgis (élément par élément) ; les noyaux restent
accessibles par les paramètres et les variables existants, pour les
mesures.

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
  `hip`), avec un nommage fixe.
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
