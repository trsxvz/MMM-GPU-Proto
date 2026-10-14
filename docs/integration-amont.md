# Intégration amont

Ce document propose un découpage des changements en contributions pour les
dépôts officiels de TFEL, de MGIS et de mfem-mgis, puis des exemples et de
mm-opera-hpc. Chaque contribution doit pouvoir être relue et testée seule.
Les corrections déjà proposées à MFEM et à AmgX sont listées à la fin.

Certains commits du prototype mêlent plusieurs sujets. Les contributions
seront donc reconstruites à partir de l'état final des branches, fichier par
fichier, plutôt que par sélection de commits.

## TFEL

L'ordre suit les dépendances : chaque contribution s'appuie sur les
précédentes. Les détails techniques sont dans [tfel.md](tfel.md).

1. **Corrections du code existant**, sans code GPU. Ce sont les corrections
   listées dans [tfel.md](tfel.md#corrections-du-code-existant). Elles
   peuvent changer les résultats sur CPU, puisque ce sont des corrections.
   Chacune est couverte par un test existant ou nouveau.
2. **Code GPU dans TFEL** (en-têtes mathématiques et matériau) :
   annotations hôte et GPU, `TFEL_DEVICE_CODE` pour les compilateurs à une
   passe, lois sans fonctions virtuelles, `TrivialGenTypeBase`,
   `getConstantObject`, valeurs propres sans récursion, contournements de
   compilateurs (par exemple les itérateurs des splines cubiques). Les
   résultats sur CPU ne changent pas. Les incompatibilités sont listées dans
   les notes de version.
3. **Code généré compatible GPU dans MFront** : annotations des classes
   générées, classe de loi compatible GPU, pilotes de l'interface
   `generic` coupés en un cœur hôte et GPU et les fonctions historiques,
   fonctions de rotation et appels de méthodes utilisables sur GPU, type de
   l'opérateur tangent vérifié avant sa conversion, paramètres transmis à la
   classe de loi et à ses variables de comportement, erreurs portables
   (`tfel::fail`, `abortResolution`, `abortIntegration`), données partagées
   (`@SharedData`), variables statiques écrites `TFEL_STATIC_VARIABLE(v)`
   (copie dans le code GPU, variable elle-même sur l'hôte). Les résultats
   des interfaces existantes ne changent pas.
4. **Interface `generic-parallel` des lois et des modèles**, avec un
   premier ensemble de backends : CUDA, HIP, algorithmes parallèles du C++
   et OpenMP sur l'hôte. Contrat de données versionné, options d'exécution,
   calcul synchrone ou asynchrone, erreurs sans coût en cas de succès,
   symbole d'espace mémoire, données partagées copiées sur chaque GPU, tests
   et documentation.
5. **Backends et compilateurs supplémentaires** : Kokkos, SYCL, déport de
   `nvc++` (`stdpar` et OpenMP), clang CUDA et HIP. La plupart des
   contournements de compilateurs vont dans cette contribution.
6. **Extensions de l'interface des propriétés matériau** : backends Kokkos
   et OpenMP, propriétés définies par des données ou appelant d'autres lois
   sur CUDA et HIP, options d'exécution en dernier argument avec le même sens
   pour tous les backends (version 2 des fonctions), évaluation asynchrone
   et fonction de compte rendu, symboles de backend, de version et d'espace
   mémoire, code des noyaux partagé par les backends, statuts dans des blocs
   réutilisés, statut 1 avec la politique `Warning`.
7. **Précision float** (option de DSL `real_type`, patchs 0112 à 0120) :
   type numérique des interfaces `generic` et `generic-parallel`, symbole
   `<f>_real_size` et versions d'API décalées, refus du float par les
   autres interfaces et les modèles du DSL `Model`, `IsAssignableTo` vers
   un flottant plus étroit, MTest qui contrôle le type, critères de
   convergence bornés à la génération avec avertissements, documentation.
   En double, le code généré ne change que par le symbole de taille.

Le contrat de données des lois
(`MFront/GenericParallelBehaviour/BatchData.h`) est versionné dès la
quatrième contribution, pour que MGIS puisse le vérifier.

### Vérifications pour chaque contribution

- la suite de tests complète de TFEL (`ctest`) ;
- les comparaisons des tests `generic-parallel` pour les configurations de
  la contribution, à `-O3` ;
- aucun changement des résultats des interfaces existantes, sauf les
  corrections de la première contribution ;
- formatage par `clang-format` limité aux lignes modifiées.

## MGIS

Les contributions de MGIS suivent la quatrième contribution de TFEL, dont
elles lisent le contrat. Les détails sont dans [mgis.md](mgis.md).

1. **Corrections** : axes des matrices de rotation, propriété matériau
   donnée par un `std::span` d'une valeur, statut d'échec des fonctions
   d'initialisation et des post-traitements point par point.
2. **Lots des lois et des modèles** : chargement et vérifications, copies
   des structures de TFEL et test de leur disposition, intégration sans
   copie, fonctions d'initialisation, post-traitements, rotations,
   basculement automatique des fonctions existantes, politique d'exécution,
   refus des fonctions d'un point.
3. **Mémoire** : emplacement des tableaux vérifié contre l'espace mémoire
   des bibliothèques, refus des accès de l'hôte aux tableaux placés sur un
   GPU, ressources mémoire et leurs en-têtes.
4. **Propriétés matériau** : chargement des propriétés `generic-parallel`
   et fonction `evaluate`.
5. **Rangements et opérations de copie** : rangement des tableaux d'un
   gestionnaire (`MaterialDataManagerInitializer::layout`), opérations de
   copie (`FieldOperations`), vue sans contrôle des valeurs gardées sur
   l'hôte (`make_unchecked_view`) et mise à jour des seules entrées d'un
   état (`updateInputValues`), dont dépendent les données des matériaux
   de mfem-mgis sur le GPU, fichiers HDF5 et liaisons C, Fortran et
   Python qui suivent le rangement.
6. **Profilage** : fonction de synchronisation appelée autour des
   sections mesurées, pour attendre les noyaux lancés sur un GPU.
7. **Documentation** : page web des lots et notes de version.
8. **Précision float** (option `MGIS_REAL_TYPE`, patchs 0064 à 0070) :
   type contrôlé au chargement de chaque point d'entrée, tests en float
   sur le CPU et avec CUDA, documentation. Elle suit la contribution 7 de
   TFEL.

Le test CUDA de MGIS demande un GPU NVIDIA. Il est désactivé par défaut
(option `enable-generic-parallel-cuda-tests`), pour que l'intégration
continue de MGIS reste sur CPU.

## mfem-mgis

La branche `gpu-proto` de mfem-mgis part du master officiel `8fac67336e`
(13/10/2026, fusion de la PR #485). L'ordre suit les dépendances :

- les contributions 1 et 2 ne demandent ni le TFEL ni le MGIS du
  prototype ;
- la contribution 3 demande la synchronisation du profilage de MGIS
  (sixième contribution de MGIS) ;
- la contribution 5 demande l'emplacement des tableaux et les rangements
  de MGIS (troisième et cinquième contributions de MGIS) ;
- à partir de la contribution 6, l'interface `generic-parallel` de TFEL
  (quatrième contribution de TFEL) et les lots de MGIS (deuxième
  contribution de MGIS).

Le plan d'ensemble est dans [plan-mfem-mgis.md](plan-mfem-mgis.md).

1. **Corrections du code existant**, sans code GPU :
   - préconditionneur ParaSails redonné à un nouvel opérateur après la
     destruction du précédent (défaut de `mfem::HypreParaSails`) ;
   - réglages de l'ILU de hypre réappliqués après chaque nouvel opérateur
     (MFEM les remet aux valeurs par défaut) ;
   - tolérance absolue des solveurs de Krylov de hypre (refusée, ou donnée
     comme tolérance relative à FGMRES) ;
   - noms courts des options déclarées par `declareDefaultOptions` (`-up`,
     `-no-up`, `petscrc_file`) ignorés par `initialize` ;
   - exception levée par l'opérateur tangent pendant la correction de
     Newton : échec de la résolution en séquentiel, message puis
     `MPI_Abort` en parallèle, au lieu de `std::terminate` ;
   - exceptions des projections L2, déclarées `noexcept`, rapportées dans
     le contexte d'exécution ;
   - prédiction séquentielle : diagonale mise à un aux degrés de liberté
     imposés, comme en parallèle ;
   - CMake : plusieurs dossiers rendus par `tfel-config`, une ligne de
     `env.sh` par dossier ;
   - documentation des constructeurs qui prennent un contexte
     d'exécution.
2. **Fonctionnalités sans GPU**, utiles sur CPU :
   - ordre des règles d'intégration choisi par l'utilisateur
     (`IntegrationRuleOrder` : avec l'ordre 3, 8 points au lieu de 27 sur
     les hexaèdres d'ordre 1) ;
   - options de BoomerAMG (seuil de connexion, grossissement,
     interpolation, relaxation, cycles), réappliquées après chaque nouvel
     opérateur ;
   - Newton inexact (tolérance du solveur linéaire adaptée au résidu,
     Eisenstat et Walker) ;
   - lisseurs de Jacobi et de Chebyshev de MFEM comme préconditionneurs ;
   - analyse de MUMPS réutilisée tant que le motif de la matrice ne change
     sur aucun processus (`ReuseAnalysis`, nombre d'analyses donné par
     `getNumberOfAnalyses`, pour MUMPS et cuDSS ; test
     `DirectSolverReuseAnalysisTest`) ;
   - paramètres du problème transmis par les constructeurs de
     `PeriodicNonLinearEvolutionProblem` ;
   - `TransientHeatTransfer` en déformations planes et en contraintes
     planes ;
   - sections de profilage (préparation du solveur linéaire, prédiction,
     fin de pas, post-traitements).
3. **Construction et configuration pour le GPU** :
   - MFEM compilé avec CUDA, sources contenant des noyaux de MFEM
     compilées comme du code CUDA, ou comme du code HIP quand MFEM utilise
     HIP (plateforme AMD, architectures de MFEM ; compilation vérifiée,
     sans essai sur un GPU AMD) ;
   - `FindMFEM` : `MFEM_DIR` prime sur `CMAKE_PREFIX_PATH` ;
   - option `--device` d'`initialize` (un GPU par processus, périphérique
     déjà configuré par le programme, périphérique créé avant PETSc) ;
   - réserve de mémoire de hypre sur le GPU ;
   - mémoire de l'hôte standard pour le périphérique configuré par
     `initialize` : avec MFEM compilé avec Umpire, les vecteurs
     temporaires des boucles de l'hôte sur les éléments passaient par
     Umpire (2 à 5 fois plus lents) ;
   - hypre sur l'hôte avec le périphérique sur le GPU
     (`--hypre-on-host`), pour les petits systèmes linéaires ;
   - solveur direct cuDSS sur le GPU (`CuDSSSolver`, enveloppe de
     `mfem::CuDSSSolver` qui le recrée quand la taille ou le motif change),
     avec la réutilisation de l'analyse de la contribution 2 ;
   - sections de profilage qui attendent les noyaux lancés.
4. **Chemin GPU du code existant** : toute la suite de tests sur le GPU
   (variable `MFEM_MGIS_DEVICE`), et les corrections trouvées ainsi :
   - accès de l'hôte déclarés à MFEM (`HostRead`, `HostReadWrite`) :
     composantes des projections L2, valeurs imposées des conditions de
     Dirichlet, export ParaView, interpolation des champs, erreur L2
     analytique ;
   - ILU de hypre limitée au niveau 0 et ParaSails refusé quand hypre
     calcule sur le GPU (défaillances de hypre).
5. **Données des matériaux dans la mémoire de MFEM** : tableaux des points
   d'intégration rangés composante par composante, synchronisés sur
   l'hôte avant chaque accès de l'hôte, éléments et décalages des espaces
   de quadrature lisibles sur le GPU, copies des états par des noyaux.
6. **Intégration des lois par lots**, sur l'hôte et sur le GPU : lots de
   MGIS, rotations des lois orthotropes, `update` et `revert` sur le GPU,
   modèles point par point (`PointWiseModel`). Sans l'option
   `--execution` (contribution 11), l'intégration point par point reste le
   défaut ; le paramètre `BehaviourIntegration` (ou la variable
   `MFEM_MGIS_BEHAVIOUR_INTEGRATION`) choisit les lots.
7. **Mécanique par noyaux** :
   - géométrie des éléments (jacobiens et leurs inverses calculés une
     fois), gradients, forces intérieures, matrices de rigidité ;
   - assemblage direct dans la matrice creuse (motif construit à partir
     des éléments, positions des entrées calculées une fois, échanges
     gardés), par additions atomiques ou dans un ordre fixe, qui rend les
     résultats identiques d'un passage à l'autre (paramètre
     `TangentOperatorAssembly`, variable
     `MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY`) ;
   - prédiction, action et diagonale de l'opérateur tangent sans
     assemblage (`MatrixFree`), préconditionneur sur l'opérateur assemblé
     (`AssembledOperator`) ;
   - autres intégrateurs de la forme (Robin, ressorts), conditions aux
     limites indépendantes des inconnues (pression imposée, source de
     chaleur) assemblées une fois et ajoutées aux forces intérieures
     calculées par noyaux.

   Sans l'option `--execution`, le calcul élément par élément reste le
   défaut ; les paramètres `GradientsComputation`,
   `InnerForcesComputation` et `TangentOperatorComputation` (ou les
   variables correspondantes) choisissent les noyaux.
8. **Autres formulations par noyaux** : thermique stationnaire et
   transitoire, endommagement micromorphe, FBar, régularisation de
   Faltus et al.
9. **Solveur de point fixe** (`FixedPoint`) : accélérations de MTest,
   méthode NewtonK0 (opérateur élastique factorisé, tolérance adaptée),
   prédiction élastique, solveurs itératifs de la fabrique.
10. **Post-traitements par noyaux** (variable
    `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION`) : intégrales, projections
    L2, valeurs aux nœuds.
11. **Réglages d'exécution** (voir [plan-execution.md](plan-execution.md)) :
    option `--execution cpu|gpu|mixed` et placement de chaque phase en
    mode `mixed` (`--finite-elements`, `--behaviours`,
    `--linear-solver-device`, `--post-processings`), valeurs par défaut
    déduites du placement (lots et noyaux, sauf avec PETSc), variables de
    placement ignorées avec un avertissement, placement de la loi d'un
    matériau (paramètre `BehaviourExecutionLocation`), copies des tableaux
    des matériaux comptées, bibliothèques de lois pour le CPU et le GPU
    (`mfem_mgis_add_behaviours_library`, `getBehavioursLibrary`). Elle
    dépend des lots, des noyaux et des post-traitements.
12. **Précision float** (voir [plan-precision.md](plan-precision.md)) :
    type réel donné par MFEM, configuration (refus d'un MGIS d'un autre
    type et d'AmgX en float), paramètres stockés en `double`, types MPI,
    tolérances des tests adaptées au float, régularisation d'Anderson
    relevée à 10 ε, et, en float, noyaux de raideur du GPU compilés sans
    contraction en multiplications-additions fusionnées. Elle suit les
    contributions de précision de TFEL et de MGIS.
13. **Facultatif, à proposer plus tard** :
    - Kokkos initialisé par mfem-mgis (option `enable-kokkos`), pour les
      lois compilées avec le backend Kokkos ;
    - solveurs AmgX, qui utilisent directement l'interface C d'AmgX :
      après la fusion de la correction proposée à AmgX et de sa recette
      spack.

Beaucoup de commits de `gpu-proto` touchent plusieurs contributions. Par
exemple, le passage de toute la suite sur le GPU porte des corrections
(contribution 1), de la configuration (3) et le chemin GPU (4), et
`LinearSolverFactory.cxx` porte des morceaux de six contributions. Comme
pour TFEL, les contributions seront reconstruites à partir de l'état
final, par morceaux de fichier.

### Vérifications pour chaque contribution

- la suite de tests complète de mfem-mgis sur CPU, en séquentiel et sur
  deux processus ;
- à partir de la contribution 4, la suite sur le GPU
  (`MFEM_MGIS_DEVICE=cuda`) ;
- à partir de la contribution 6, la suite avec les variables qui forcent
  les lots et les noyaux (`MFEM_MGIS_BEHAVIOUR_INTEGRATION=Batch`,
  `MFEM_MGIS_*_COMPUTATION=Kernels`), sur CPU et sur GPU ;
- les exemples de mfem-mgis-examples et les cas de mm-opera-hpc suivis par
  le prototype, avec les mêmes résultats sur CPU et sur GPU à l'arrondi
  près ;
- un test pour chaque correction et chaque fonctionnalité (certaines n'en
  ont pas encore, par exemple les paramètres de
  `PeriodicNonLinearEvolutionProblem`) ;
- formatage par `clang-format` limité aux lignes modifiées (le dépôt n'est
  pas formaté dans son ensemble).

## mfem-mgis-examples et mm-opera-hpc

Les changements des exemples (31 patchs depuis `600da70`) et de
mm-opera-hpc (33 patchs depuis `49147f2`) suivent les contributions de
mfem-mgis dont ils se servent. Ils seront proposés après elles, par
thème :

- options de mfem-mgis et de la loi dans chaque exemple et chaque cas
  (bibliothèque de lois, périphérique, solveurs, profilage), puis
  `--execution` et les bibliothèques de lois pour le CPU et le GPU
  (contribution 11) ;
- synchronisations sur l'hôte avant les lectures de l'hôte (contribution
  5) ;
- choix du solveur non linéaire (Newton, NewtonK0) et de la tolérance
  adaptée (contributions 2 et 9), solveurs directs et réutilisation de
  l'analyse (contributions 2 et 3), AmgX (contribution 13) ;
- exemples et cas en float : types MPI, options typées, critères et
  tolérances adaptés (contribution 12) ;
- tests sur le GPU, enregistrés seulement quand la loi pour le GPU est
  construite ;
- dans mm-opera-hpc, tests comparés aux solutions de référence (Cast3M,
  solution analytique de la bulle) et corrections relevées par la
  relecture : aide `-h` sans échec, unité de la pression de la bulle,
  fichiers de sortie ouverts après les options et écrits par le seul
  processus principal.

## Déjà proposé en amont

Corrections de MFEM et d'AmgX proposées pendant le prototype à
leurs dépôts officiels (ouvertes au 13/10/2026) :

- MFEM #5604 : ressources partagées d'`AmgXSolver` (durée de vie,
  initialisation séquentielle) ;
- MFEM #5605 : ordre des sommets des éléments de bord d'un `SubMesh`,
  repris de l'élément de bord parent ;
- MFEM #5615 : `HypreParaSails::SetOperator` prend le communicateur du
  nouvel opérateur. mfem-mgis contourne ce défaut dans la contribution 1 ;
- AmgX #386 : multigrille classique avec CUDA 13 (« invalid device
  function », reproducteur `reproducers/amgx-classical-cuda13`).

La reconstruction de MFEM et d'AmgX dans l'environnement attend leur
fusion (voir `deps/spack/README.md`).
