# Plan détaillé pour mfem-mgis

Statut : plan établi le 09/10/2026, relu de façon indépendante le même
jour, puis confronté à tout le code de mfem-mgis (82 fonctionnalités
inventoriées ; oublis et couvertures partielles intégrés). Avancement
au 09/10/2026 au soir : L0, L1 et L2 faits et relus (voir le
[journal](../JOURNAL.md)). L3 avancé : gradients, résidu et tangente
assemblée par noyaux, par les points d'entrée de
`NonLinearEvolutionProblemImplementation` (`Mult`, `GetGradient`), action
et diagonale de la tangente sans assemblage (option `MatrixFree`, avec un
préconditionneur de Jacobi seulement) ; restent les blocs nodaux,
l'opérateur propre complet et l'assemblage partiel de MFEM. Tranche
minimale de L6 faite (loi `cuda` intégrée sur les tableaux du GPU, fin de
pas sur le GPU, validation par le périphérique `debug`, profil sans copie
cachée), prédiction de L4 par noyaux faite. Sur un GPU : poutre r4 en
3,3 s contre 17,5 s sur 8 cœurs, polycristal en 13,2 s contre 31,7 s
(`results/L6-*`). Au 10/10/2026 : action de la tangente sans
assemblage par un noyau fusionné et dérivées des inconnues par un noyau
propre (`results/L4-*`), préconditionneur construit sur la tangente
assemblée (`AssembledOperator`) ; L5 vague 1 en grande partie faite :
solveur `FixedPoint` (méthodes `ElasticOperator` et `ModifiedNewton`,
toutes les accélérations de MTest, prédiction par l'AMG de K0), options
de BoomerAMG ; sur 8 cœurs le point fixe bat Newton (×1,5 à ×1,7), sur le
GPU égalité avec le lisseur de Chebyshev (`results/L5-*`). `NewtonK0` fait
(méthode du solveur `FixedPoint`, tout solveur itératif de MFEM de la
fabrique préconditionné par le solveur de K0) : sur le GPU, ×1,2 sur la
poutre et ×1,06 sur le polycristal, sur 8 cœurs ×1,19 sur le
polycristal, mêmes itérations que Newton. Newton inexact fait (option
`AdaptiveLinearSolverTolerance`, Eisenstat-Walker) : avec NewtonK0, ×1,59
contre Newton sur le GPU (poutre). Restent en L5 : précision mixte,
recherche linéaire, lisseurs sans matrice, vagues 2 et 3. Les choix
de fond sont dans [plan-mfem-mgis.md](plan-mfem-mgis.md) ; la liste des
solveurs dans [solveurs-gpu.md](solveurs-gpu.md). Ce document découpe le
travail en lots, dit quels fichiers changent, ce qui est testé et quand
un lot est fini.

Base : mfem-mgis master `a011d97d6f` (06/10/2026), dépôt `mfem-mgis/`.
Les références `fichier:ligne` renvoient à cette base. TFEL et MGIS sont
ceux du prototype (MGIS installé dans `install/mgis` par `build-mgis.sh`),
MFEM 4.10 avec CUDA (et les backends OpenMP et libCEED) vient de
l'environnement `mmm-gpu`. Le plan s'applique à mfem-mgis, à
mfem-mgis-examples et à mm-opera-hpc.

## Sommaire

1. [Le code aujourd'hui](#le-code-aujourdhui)
2. [Principes](#principes)
3. [Paramètres proposés](#paramètres-proposés)
4. [Vue d'ensemble des lots](#vue-densemble-des-lots)
5. [L0. Base de travail](#l0-base-de-travail)
6. [L1. Données aux points dans la mémoire de MFEM](#l1-données-aux-points-dans-la-mémoire-de-mfem)
7. [L2. Lois intégrées par lots sur CPU](#l2-lois-intégrées-par-lots-sur-cpu)
8. [L3. Noyaux par lots et opérateur propre](#l3-noyaux-par-lots-et-opérateur-propre)
9. [L4. Opérateur tangent et prédiction](#l4-opérateur-tangent-et-prédiction)
10. [L5. Solveurs](#l5-solveurs)
11. [L6. GPU](#l6-gpu)
12. [L7. Tout mfem-mgis](#l7-tout-mfem-mgis)
13. [L8. Sorties et post-traitements](#l8-sorties-et-post-traitements)
14. [L9. Exemples et mesures](#l9-exemples-et-mesures)
15. [L10. Documentation et propositions amont](#l10-documentation-et-propositions-amont)
16. [Stratégie de test](#stratégie-de-test)
17. [Campagne de mesures](#campagne-de-mesures)
18. [Questions à trancher](#questions-à-trancher)
19. [Risques](#risques)

## Le code aujourd'hui

Ce que le chemin par lots doit remplacer ou contourner, relevé dans le
code (fichiers de `src/` et `include/MFEMMGIS/`).

### Boucle de résolution

- `NonLinearEvolutionProblemImplementation<parallel>` hérite de
  `(Par)NonlinearForm` et de `NonLinearEvolutionProblemImplementationBase`.
  Il porte un `MultiMaterialNonLinearIntegrator* const` (classe `final`),
  ajouté au formulaire par `AddDomainIntegrator`
  (`NonLinearEvolutionProblemImplementation.cxx:335` et `:557`).
- Une itération de Newton fait **trois passes élément par élément sur
  l'hôte** :
  1. `integrate(u, it, dt)` (`...Implementation.cxx:485-514` et
     `:715-745`) : pour chaque élément, valeurs élémentaires, puis
     intégration de la loi point par point ;
  2. `Mult` du formulaire : `AssembleElementVector`, qui lit les forces
     thermodynamiques ;
  3. `GetGradient` du formulaire : `AssembleElementGrad`, qui lit `K`.
- **`Mult` n'intègre pas la loi.** L'intégration est déclenchée par une
  action du solveur non linéaire de mfem-mgis, appelée à chaque nouvelle
  estimation (`NonLinearSolverBase.cxx:39-58`). Un solveur non linéaire
  extérieur (MFEM, SUNDIALS, PETSc), qui n'appelle que `Mult`, calculerait
  le résidu avec des contraintes périmées. Seule exception, le mode PETSc :
  `AssembleElementVector` intègre la loi, et un échec de la loi arrête le
  programme (`MFEM_VERIFY`, `MultiMaterialNonLinearIntegrator.cxx:290-297`).
- Le résidu passe par `computeResidual` (`NonLinearSolverBase.cxx:135-140`),
  la tangente par `getJacobian` (`:142-152`), c'est-à-dire par
  `oper->Mult` et `oper->GetGradient` : un opérateur propre peut s'y
  brancher.
- Le Newton (`mfem_mgis::NewtonSolver`) est propre à mfem-mgis : il dérive
  de `mfem::IterativeSolver` par `NonLinearSolverBase`, pas de
  `mfem::NewtonSolver`. Il divise la correction par 2 quand la loi échoue
  (`NewtonSolver.cxx:166-195`), sans autre recherche linéaire.
- Le solveur linéaire reçoit la tangente à chaque itération
  (`NewtonSolver.cxx:245`, `prec` y est le solveur linéaire) ; son
  `SetOperator` reconstruit le préconditionneur (Krylov de MFEM et de
  hypre).
- Le seul solveur non linéaire de la fabrique est `Newton`
  (`NonLinearSolverFactory.cxx`), plus `PetscNonlinearSolver` quand PETSc
  est demandé.
- La prédiction (`...Implementation.cxx:98-233`) passe par une
  `BilinearForm` et une `LinearForm` (`StiffnessMatrixIntegrator`,
  `InternalForcesIntegrator`), élément par élément. Elle se règle en C++
  seulement (`setPredictionPolicy`), sur quatre axes : stratégie (dont
  `CONSTANT_GRADIENTS_INTEGRATION_PREDICTION`), opérateur, opérateur
  d'intégration, incrément de temps nul
  (`AbstractNonLinearEvolutionProblem.hxx:39-205`).
- Les conditions de Dirichlet écrivent les valeurs imposées dans `u` sur
  l'hôte (`UniformDirichletBoundaryCondition.cxx:64-79`).
- Fin de pas : `update` et `revert` de MGIS sur l'hôte, matériau par
  matériau (`BehaviourIntegratorBase.cxx:779-785`).
- Profilage : des sections de mesure existent déjà (`CatchTimeSection`,
  `ProfilingSection`, sections « NS::... »).

### Intégrateurs, matériaux, données aux points

- Un intégrateur de comportement **est** le matériau :
  `BehaviourIntegratorBase` dérive de `Material`, qui dérive de
  `mgis::behaviour::MaterialDataManager`.
- Un attribut d'élément peut porter **plusieurs** intégrateurs, dont les
  contributions sont sommées (`MultiMaterialNonLinearIntegrator.cxx:272-276`
  et `:301-318`).
- `Material` construit le gestionnaire de MGIS sans initialiseur
  (`Material.cxx:48-61`) : tableaux alloués par MGIS sur l'hôte, rangés
  point par point.
- L'intégration est faite point par point
  (`BehaviourIntegratorBase::performsLocalBehaviourIntegration`,
  `BehaviourIntegratorBase.cxx:672-743`), avec une
  `BehaviourDataView` et des tampons partagés. Elle demande une
  bibliothèque `generic`.
- Le vrai code des intégrateurs est dans
  `StandardBehaviourIntegratorCRTPBase.ixx` : gradients (`:41-142`),
  forces internes (`:154-219`), tangente (`:222-290`), avec
  `CalcPhysDShape` et des `DenseMatrix` membres. Les matrices `B` des
  classes concrètes sont **générées** (`code-generation/behaviour-integrator.cxx`,
  GiNaC). Indexation élémentaire `u[ni + c * nnodes]`, ordre natif des
  degrés de liberté de MFEM.
- **Règle d'intégration** : `IntRules.Get(geom, 2 * t.OrderGrad(&e))`
  (classes générées, par exemple
  `IsotropicTridimensionalStandardSmallStrainMechanicsBehaviourIntegrator.cxx:16-20`),
  soit 27 points sur l'hexaèdre d'ordre 1. Elle est choisie élément par
  élément : un maillage mixte a plusieurs règles.
- **Rotations des lois orthotropes** : les gradients sont tournés avant
  l'intégration et les blocs de la tangente après ; les forces
  thermodynamiques restent dans le repère du matériau (elles servent
  d'état au pas suivant) et ne sont tournées qu'à l'assemblage, dans une
  copie (`StandardBehaviourIntegratorCRTPBase.ixx:134-139` et
  `:208-210`). Matrices par point données par un `variant`
  (`Material.hxx:99-106`), rotations par pointeurs de fonction de MGIS.
- `PartialQuadratureSpace` : points contigus par élément, éléments d'un
  matériau dans l'ordre croissant, mais décalages rangés dans une
  `std::unordered_map` (`PartialQuadratureSpace.hxx:230-236`).
- `PartialQuadratureFunction` : `std::vector` ou mémoire empruntée, vues
  `std::span`, accès `values[o * stride + begin]`.
  `getIntegrationPointValues` rend un `std::span` contigu
  (`PartialQuadratureFunction.hxx:198-228` et `:323-352`), appelé une
  quarantaine de fois ; les intégrateurs lisent les données par
  `o * stride`. **Tout le chemin actuel suppose le rangement point par
  point.**
- Propriétés matériau et variables externes : évaluateurs
  (`AbstractQPEvaluator`), évalués sur l'hôte dans `setup`
  (`BehaviourIntegratorBase.cxx:360-670`).

### Solveurs linéaires

- Fabrique `LinearSolverFactory.cxx` : Krylov de MFEM (séquentiel et
  parallèle), `HyprePCG`, `HypreGMRES`, `HypreFGMRES`, `MUMPSSolver`,
  `UMFPackSolver` ; préconditionneurs `HypreBoomerAMG` (stratégies
  `Elasticity`, `System`, `None` ; `System` par défaut quand l'espace a
  plusieurs composantes), `HypreEuclid`, `HypreILU`, `HypreParaSails`,
  `HypreDiagScale`. Tous les préconditionneurs et solveurs directs
  demandent une matrice assemblée.
- `SystemHypreBoomerAMG` rappelle `SetSystemsOptions` à chaque
  `SetOperator` (`LinearSolverFactory.cxx:46-50`), ce qui remet deux
  réglages de l'AMG à leur valeur par défaut.
- Deux défauts relevés au passage : `HypreGMRES` teste un paramètre
  `AbsoluteTolerance` qui n'est pas dans la liste autorisée (code mort),
  `HypreFGMRES` aussi, et appelle `SetTol` au lieu d'une tolérance
  absolue (`mfem::HypreFGMRES` n'a pas de `SetAbsTol` en MFEM 4.9 ni 4.10 : passer
  par `HYPRE_ParCSRFlexGMRESSetAbsoluteTol`).

### GPU

- Aucun usage de `mfem::Device`, `AssemblyLevel`, `forall`,
  `ElementRestriction`, `QuadratureInterpolator`, `GeometricFactors` ou
  `DofToQuad`. `mfem_mgis::initialize` ne crée pas de périphérique.
- `mfem::forall` ne lance un noyau CUDA que si le fichier est compilé
  par nvcc (`__CUDACC__`, `general/forall.hpp:789-794`) ; sinon le corps
  s'exécute sur l'hôte. Les fichiers de mfem-mgis qui contiennent des
  noyaux doivent donc être compilés comme du CUDA quand MFEM l'est.

### Points d'extension et comportements à préserver

Relevés en confrontant le plan à tout le code :

- **Intégrateurs du formulaire autres que celui de mfem-mgis.** Le
  problème hérite publiquement de `(Par)NonlinearForm` : l'utilisateur
  peut y ajouter ses intégrateurs. mfem-mgis y met lui-même la pression
  imposée (intégrateur de bord) et la source de chaleur (intégrateur de
  domaine à marqueurs, `UniformHeatSourceBoundaryCondition.cxx:208` et
  `:223`) ; l'exemple ex8 de mfem-mgis-examples y met une condition de
  Robin (intégrateur de faces de bord dont la tangente dépend de `u`,
  `ex8/src/RobinBC.cxx`). Le mode
  `UseMultiMaterialNonLinearIntegrator=false`
  (`NonLinearEvolutionProblemImplementationBase.cxx:18-42`) se passe de
  l'intégrateur de mfem-mgis.
- **Degrés de liberté imposés** : ceux du formulaire
  (`SetEssentialTrueDofs`), pas seulement ceux des conditions de
  Dirichlet. Le problème périodique les pose directement
  (`PeriodicNonLinearEvolutionProblem.cxx:100`, `:192`, `:235`, `:291`) :
  `getEssentialDegreesOfFreedom()` y est vide. L'interface des conditions
  de Dirichlet est publique et écrit `u` sur l'hôte.
- **Données des matériaux lues ou écrites sur l'hôte par l'API publique** :
  `getGradient`, `getThermodynamicForce`, `getInternalStateVariable`
  (`Material.cxx:269-370`), `MechanicalPostProcessings`,
  `makeGridFunction`, `computeL2Projection`, `StateManager`,
  `makeGradientEvaluator`, post-traitements donnés par `std::function`,
  initialisation des variables internes par l'utilisateur. Toutes ces vues
  reposent sur un pas entre points (`data_stride`).
- **Propriétés matériau et variables externes posées par l'utilisateur
  dans sa propre mémoire hôte**, par l'API de MGIS (tests
  `StationaryNonLinearHeatTransferTest`, `TwoLayerCube` ; mm-opera-hpc
  `cermet/cermet.cxx:118-135`).
- **Convergence des solveurs linéaires** : `hasConverged` et
  `getNumberOfIterationsAtConvergence` (`SolverUtilities.cxx:103-163`)
  déclarent convergé tout solveur qu'elles ne connaissent pas.
- **Échecs** : un attribut sans intégrateur arrête aujourd'hui le calcul
  (`MultiMaterialNonLinearIntegrator.cxx:30-53`) ; l'intégration s'arrête
  au premier élément en échec ; le statut 0 (« résultats peu fiables »)
  compte comme un succès.
- **Prédiction** : elle ajoute la contribution des conditions aux limites
  (`addLinearFormIntegrators`, `...Implementation.cxx:160`) et respecte
  `DiscardLinearSolverFailure`. `getLinearizedOperators` est une API
  publique (utilisée par `PredictionTest`).
- **Séquentiel** : la tangente assemblée est une `SparseMatrix` (UMFPACK,
  Krylov de MFEM) ; les préconditionneurs de hypre y sont refusés
  (`LinearSolverFactory.cxx:245-300`).
- **Dépendances entre problèmes** :
  `resolveBehaviourIntegratorsDependencies`
  (`NonLinearEvolutionProblem.cxx:795-854`) relie les intégrateurs de deux
  problèmes, qui doivent garder la même règle d'intégration.
- **Régularisation par gradient implicite** (`L2Projection.cxx:624-702`)
  et intégrales (`computeMeasure`, `computeIntegral`, énergies,
  résultante) : boucles hôte sur les points.
- **Ce qui n'a rien à faire** : la fabrique d'évaluateurs (presque
  entièrement en commentaire), les dépendances entre modèles (non
  implémentées), la reprise de l'état des matériaux et les liaisons
  Python (absentes), la lecture des maillages, la gestion du temps.

### Ce que l'assemblage partiel de MFEM 4.10 impose

Vérifié dans les sources de MFEM 4.9 et 4.10 (`fem/nonlinearform.cpp`,
`fem/nonlinearform_ext.cpp`, `fem/pnonlinearform.cpp`) :

- `NonlinearForm::Mult` rend la main juste après l'extension : les
  intégrateurs de bord (`AddBoundaryIntegrator`) sont **ignorés sans
  message** (pression imposée, `UniformImposedPressureBoundaryCondition`).
  Un intégrateur de domaine à marqueurs sans assemblage partiel arrête le
  calcul (`AssemblePA` non implémenté) : c'est le cas de la source de
  chaleur (`UniformHeatSourceBoundaryCondition`).
- Les vecteurs élémentaires de l'extension sont en ordre lexicographique,
  pas dans l'ordre natif des noyaux actuels.
- `GetGradient` rend un opérateur sans matrice, conditions de Dirichlet
  comprises ; `ParNonlinearForm::GetLocalGradient` refuse l'extension :
  pas de `HypreParMatrix` par ce chemin.

Conséquence : l'opérateur propre (point d'entrée B du plan) est le chemin
principal. L'assemblage partiel de MFEM (point d'entrée A) reste une
option, refusée avec un message clair si le problème a des intégrateurs
de bord ou à marqueurs.

## Principes

- **Le chemin actuel ne change pas** et reste le chemin par défaut. Les
  tests existants ne changent pas.
- **Travail dans `mfem-mgis/`, `mfem-mgis-examples/` et `mm-opera-hpc/`,
  branche `gpu-proto` de chacun** (bases `a011d97d6f`, `8397eec`,
  `8bd56eb`), copiés dans ce dépôt par `scripts/sync-repository.sh`. Messages de commit
  en anglais. Un lot fini est un ensemble de commits relus.
- **Code nouveau regroupé** : `include/MFEMMGIS/Batch/` et `src/Batch/`
  (nom à discuter). Les classes existantes gagnent des points
  d'extension, pas de code GPU.
- **Noyaux écrits une fois**, avec `mfem::forall`, génériques sur
  l'hypothèse et le type de gradient. Ils servent sur CPU et sur GPU, et
  aux deux points d'entrée. Les fichiers de `src/Batch/` sont compilés
  comme du CUDA (langage CUDA de CMake) quand MFEM est compilé avec CUDA,
  comme MFEM le fait pour ses propres sources.
- **Toutes les variantes, choisies par option**, comparées sur les mêmes
  cas (rangement, point d'entrée, opérateur tangent, règle d'intégration,
  solveur).
- **Chaque option a un test qui échoue si l'option est ignorée.**
- **Pas d'échec silencieux** : une configuration non prise en charge
  (intégrateur de bord avec l'assemblage partiel, stratégie `Elasticity`
  sur GPU, loi `generic` sur un périphérique GPU, rangement composante par
  composante sur le chemin actuel...) est refusée avec un message.

## Paramètres proposés

Noms à discuter, dans le style des paramètres existants.

| Où | Paramètre | Valeurs |
|---|---|---|
| problème | `AssemblyStrategy` | `Legacy` (défaut, chemin actuel), `BatchedOperator` (opérateur propre), `PartialAssembly` (assemblage partiel de MFEM) |
| problème | `QuadratureDataLayout` | `PointByPoint` (défaut), `ComponentByComponent` (refusé avec `Legacy`) |
| problème | `TangentOperator` | `Assembled` (option c, défaut), `MatrixFree` (option a), `ElementMatrices` (assemblage par élément), `ElasticOperator` (option b) |
| problème | `IntegrationRuleOrder` | ordre de la règle d'intégration ; par défaut celui d'aujourd'hui (`2 * OrderGrad`) |
| problème | `Prediction` | dictionnaire reprenant les quatre axes actuels (`Strategy`, `Operator`, `IntegrationOperator`, `NullTimeIncrement`), plus la stratégie nouvelle `Extrapolation` |
| solveur non linéaire | `Newton` | options nouvelles : `AdaptiveLinearSolverTolerance` (Eisenstat-Walker), `LineSearch`, `PreconditionerUpdate` |
| solveur non linéaire | `FixedPoint` | `Method` (`ElasticOperator`, `ModifiedNewton`), `ElasticOperatorUpdate`, `LinearSolver` (exact ou nombre de cycles d'AMG), `Acceleration` (`None`, `Anderson`, `Cast3M`, sécantes, `IronsTuck`, `Steffensen`, delta carrés), `AccelerationDepth`, `AccelerationDamping`, `AccelerationTrigger`, `AccelerationPeriod`, `OperatorPrecision` (voir [L5](#point-fixe-et-méthodes-voisines)) |
| solveur non linéaire | `Newton` | `PreconditionerOperator` : `Tangent` (défaut) ou `Elastic` (AMG de `K0` gelé, méthode `NewtonK0`) |
| solveur non linéaire | `LBFGS` (préconditionné, à écrire), `KINSOL` | options propres |
| préconditionneur | `OperatorJacobi`, `OperatorChebyshev`, `NodalBlockJacobi`, `LOR` | sans matrice |
| préconditionneur | `HypreBoomerAMG` | options nouvelles : grossissement, interpolation, lisseur, balayages, seuil, niveaux |
| préconditionneur | `ElasticOperator` | enveloppe : préconditionneur construit une fois sur l'opérateur élastique assemblé, quel que soit l'opérateur du Krylov |
| solveur, préconditionneur | `Ginkgo*`, `AmgX`, `CuDSS` | selon [solveurs-gpu.md](solveurs-gpu.md) |
| application | `--device` | configuration de `mfem::Device` (`cpu`, `cuda`, `debug`...), ajoutée à `declareDefaultOptions` |

`PreconditionerUpdate` : `EachIteration` (comportement actuel),
`EachTimeStep`, `Never`, `EveryNIterations`. Avec `Never`, le
préconditionneur est enveloppé dans un objet dont `SetOperator` ne
reconstruit rien.

**Les options de l'opérateur tangent, précisément :**

- **option c** : `TangentOperator=Assembled`, Newton, solveur et
  préconditionneur quelconques sur la tangente assemblée ;
- **option a** : `TangentOperator=MatrixFree`, Newton, Krylov de MFEM,
  préconditionneur `ElasticOperator` (AMG sur l'opérateur élastique) ou
  sans matrice (`OperatorJacobi`, `OperatorChebyshev`,
  `NodalBlockJacobi`) ;
- **option b** : `TangentOperator=ElasticOperator` et solveur
  `FixedPoint` : l'opérateur linéaire est l'opérateur élastique, assemblé
  une fois ; chaque itération applique une résolution approchée (un
  cycle d'AMG, ou un Krylov court) à `F(u)`, avec l'accélération
  d'Anderson. Le Newton avec `TangentOperator=ElasticOperator` est la même
  méthode sans Anderson et avec un Krylov complet : gardé comme variante ;
- `TangentOperator=ElementMatrices` : variante de l'option a, avec les
  matrices élémentaires gardées.

## Vue d'ensemble des lots

| Lot | Contenu | Dépend de | Taille | Fin |
|---|---|---|---|---|
| L0 | construction des trois dépôts, compatibilité MFEM 4.10, référence, périphérique, petites corrections | environnement `mmm-gpu` | S | tests existants verts avec MFEM 4.10 et CUDA, temps de référence |
| L1 | données aux points dans la mémoire de MFEM, deux rangements | L0 | M | suite verte, rangement composante par composante sur les données |
| L2 | lois intégrées par lots sur CPU | L1 | M | mêmes résultats, gain mesuré avec `stlpar` et `openmp` |
| L3 | noyaux par lots, opérateur propre, assemblage partiel | L2 | L | résidu et tangente égaux au chemin actuel, suite verte |
| L4 | opérateur tangent (a, b, c, éléments), prédiction, règle d'intégration | L3 | M | les options sur les cas de référence |
| L5 | solveurs en trois vagues, point fixe et ses variantes | L3, L4 | L | vague 1 et point fixe complets |
| L6 | GPU | L3 (tranche minimale), puis L4, L5 | L | mêmes résultats qu'en CPU, aucune copie cachée en un processus |
| L7 | toutes les hypothèses, cinématiques et modèles | L3 à L6 | L | suite complète de mfem-mgis, complétée des tests manquants, sur le chemin par lots |
| L8 | sorties et post-traitements | L1, L6 | S | sorties identiques, coût mesuré |
| L9 | exemples et mesures | L6 | M | campagne de mesures publiée |
| L10 | documentation, propositions amont | tous | M | documents et branches prêts |

Tailles : S (quelques jours), M (une à deux semaines), L (plusieurs
semaines). Chemin critique : L0, L1, L2, L3, tranche GPU de L6, L4, L5.
La tranche GPU minimale (petites déformations 3D, une loi, option c puis
a) arrive tôt pour lever les risques du GPU avant de généraliser.

## L0. Base de travail

Objectif : construire mfem-mgis dans l'environnement du prototype, avec
la référence de temps du chemin actuel.

Travail :

- Branche `gpu-proto` dans `mfem-mgis/`. Script `build-mfem-mgis.sh` à
  la racine (comme `build-mgis.sh`) : MFEM 4.10 avec CUDA de `mmm-gpu`,
  MGIS de `install/mgis`, TFEL de `install/tfel-mmm-gpu`, construction,
  `make check`.
- CMake de mfem-mgis : langage CUDA activé quand MFEM est compilé avec
  CUDA, pour les fichiers qui contiendront des noyaux.
- **Compatibilité avec MFEM 4.10** : mfem-mgis n'a été utilisé qu'avec
  4.8 et 4.9 ; corriger ce que 4.10 casse (API retirées ou dépréciées),
  en commits séparés, proposables en amont.
- mfem-mgis-examples et mm-opera-hpc construits de la même façon
  (`build-mfem-mgis.sh` construit les trois), leurs tests verts.
- Vérifier que la suite actuelle passe avec MFEM compilé pour CUDA mais
  périphérique `cpu`. MFEM de l'environnement est compilé avec Umpire :
  ses types de mémoire par défaut sont alors ceux d'Umpire.
- Entiers de hypre : la recette spack de mfem-mgis demande hypre `+int64`
  par défaut, que hypre refuse avec CUDA (« use +mixedint ») ;
  l'environnement utilise hypre en entiers 32 bits. À documenter.
- Option `--device` : `mfem::Device` créé par `mfem_mgis::initialize` (ou
  par l'application) à partir d'une chaîne (`Config.cxx:174-222`). Tester
  l'ordre avec `mfem::Hypre::Init`. Tout préconditionneur de hypre doit
  être créé après le périphérique (les valeurs par défaut de BoomerAMG
  dépendent du GPU, `linalg/hypre.cpp:5029`).
- Mesures par phase : étendre les sections de profilage existantes
  (`CatchTimeSection`, `ProfilingSection`) à toutes les phases
  (intégration, résidu, tangente, préconditionneur, résolution, fin de
  pas), avec une synchronisation du GPU avant chaque lecture d'un temps.
- Corriger les deux défauts de `HypreGMRES` et `HypreFGMRES` (petits,
  indépendants, proposables en amont tout de suite).

Tests : suite existante ; un test du paramètre `--device` (le
périphérique est celui demandé).

Fin : suite verte, temps de référence notés pour les cas de la
[campagne de mesures](#campagne-de-mesures).

## L1. Données aux points dans la mémoire de MFEM

Objectif : les tableaux des matériaux vivent dans une mémoire choisie par
mfem-mgis (`mfem::Memory`, hôte ou GPU), dans l'un des deux rangements.
MGIS les utilise sans copie.

Travail :

- `Material` : construire le `MaterialDataManager` avec un
  `MaterialDataManagerInitializer` (`Material.cxx:48-61`). `Material`
  dérive directement de `MaterialDataManager`, construit avant les
  membres : la mémoire doit être tenue par une base déclarée avant
  `MaterialDataManager` (idiome de la base qui porte les membres), ou
  créée par une fabrique avant le matériau. Les tableaux (gradients,
  forces thermodynamiques, variables internes, énergies des deux états,
  blocs de la tangente) sont des `mfem::Vector` ; leurs pointeurs sont
  donnés à MGIS avec l'emplacement (`HOST` ou `DEVICE`, d'après
  `mfem::Device`) et le rangement (`QuadratureDataLayout`). Variante à
  comparer : tableaux alloués par MGIS avec une ressource mémoire
  (`MGIS/MemoryResources`), vus par MFEM.
- **Validité des copies** : un `mfem::Memory` a une copie hôte et une
  copie périphérique, avec des drapeaux de validité que MFEM tient à jour.
  MGIS écrit par des pointeurs bruts : les pointeurs donnés à MGIS
  doivent être obtenus par `Write` ou `ReadWrite` (sur le périphérique)
  juste avant chaque appel, et tout accès hôte doit passer par `HostRead`
  ou `HostReadWrite`. Cette politique est écrite et testée dans ce lot.
- `PartialQuadratureSpace` : ajouter, à côté de la table de hachage, des
  tableaux `mfem::Array<int>` des éléments du matériau et des décalages,
  lisibles sur le périphérique. La table reste pour le code actuel. Une
  refonte de cette classe est en cours en amont (voir
  [Risques](#risques)).
- **Rangement composante par composante** : le chemin actuel suppose le
  rangement point par point partout (accès contigus, `o * stride`). Ce
  rangement est donc refusé avec `AssemblyStrategy=Legacy`, et réservé
  au chemin par lots. Fait dans L1 : le rangement se choisit à la construction
  d'un matériau ; les intégrateurs du chemin actuel construisent toujours
  leurs matériaux point par point, et les accesseurs et évaluateurs à pas
  refusent l'autre rangement avec un message. Le paramètre du problème
  (`QuadratureDataLayout`) et son refus avec `AssemblyStrategy=Legacy`
  viennent avec le chemin par lots (L2, L3). Dans ce lot, seules les données des matériaux et
  leur description changent ; les accès du chemin actuel restent point
  par point.
- **Accès hôte par l'API publique** (accesseurs de `Material`,
  post-traitements, projections, `StateManager`, évaluateurs de
  gradients, initialisation des variables internes) : chaque accesseur
  synchronise vers l'hôte (`HostRead`) ; chaque écriture hôte est
  déclarée (`HostReadWrite`) pour que la copie du périphérique soit
  invalidée. Avec le rangement composante par composante, ces vues
  (fondées sur un pas entre points) sont refusées, ou servies par une
  copie transposée.
- **Propriétés matériau et variables externes dans la mémoire hôte de
  l'utilisateur** (API de MGIS) : copiées vers la mémoire de la loi à
  chaque `setup` quand la loi calcule sur GPU, ou refusées avec un message
  qui dit comment faire.
- Fin de pas : `update` et `revert` passent par les surcharges de MGIS qui
  prennent des `FieldOperations`. MGIS fournit `HostFieldOperations` et
  des opérations CUDA, HIP, SYCL, Kokkos et OpenMP en en-têtes, à compiler
  par mfem-mgis avec la chaîne correspondante. Variante à comparer :
  `MFEMFieldOperations`, écrites avec `mfem::forall` et les copies de
  `mfem::Memory` (compilées comme du CUDA quand MFEM l'est).

Tests : suite complète inchangée (rangement point par point, données
dans `mfem::Memory`) ; refus du rangement composante par composante par
les accesseurs et les évaluateurs fondés sur un pas entre points (le
refus avec `Legacy` vient avec le paramètre du problème, en L2 ou L3) ;
accesseurs hôte après un calcul sur le périphérique (valeurs à jour) ;
tests unitaires des opérations de fin de pas dans les deux rangements ;
test de la politique de validité sous le périphérique `debug` de MFEM.

Fin : suite verte, données dans la mémoire de MFEM.

## L2. Lois intégrées par lots sur CPU

Objectif : remplacer l'intégration point par point par un appel par
matériau, sans toucher à l'assemblage. Premier gain mesurable : les
backends `stlpar` et `openmp` de TFEL sur CPU.

Travail :

- Découper `integrate` (`...Implementation.cxx:485-514`) :
  1. passe élément par élément, qui calcule seulement les gradients (code
     actuel de `implementIntegrate` sans l'appel de la loi), gradients
     macroscopiques et mise à jour des variables externes à partir des
     inconnues compris (trait `updateExternalStateVariablesFromUnknownsValues`),
     puis rotation des gradients ;
  2. pour chaque matériau, un appel `mgis::behaviour::integrate` sur une
     vue du matériau (`make_view`), avec la politique d'exécution choisie ;
  3. rotation des seuls blocs de la tangente, sur vues. **Les forces
     thermodynamiques ne sont pas tournées** : elles restent dans le
     repère du matériau, l'assemblage actuel les tourne dans une copie.
- Plusieurs intégrateurs sur un même attribut : un appel par intégrateur
  (chacun est un matériau).
- Les bibliothèques `generic` restent acceptées : avec un espace de
  travail, MGIS les intègre point par point (référence). Les bibliothèques
  `generic-parallel` sont choisies par le chemin de la bibliothèque.
- Propriétés matériau et variables externes : les valeurs évaluées dans
  `setup` sont données à la vue (uniformes ou tableaux).
- Échecs : MGIS rend le premier point en échec et son statut ; le
  comportement actuel (division de la correction de Newton, sous-pas)
  doit être conservé. Changement de comportement à documenter : par lots,
  tous les points sont intégrés (aujourd'hui l'intégration s'arrête au
  premier élément en échec), donc l'état change au-delà du point en
  échec ; l'itération est de toute façon rejetée. Le statut 0 (« résultats
  peu fiables ») et la réduction du statut entre processus MPI sont gardés
  tels quels.
- **Propriétés matériau de début de pas** : les lois `generic` ne lisent
  que celles de fin de pas, mais le traitement par lots de MGIS contrôle
  les deux états (`mgis/src/Integrate.cxx`, `checkMaterialStateView`). Un
  problème qui ne pose que celles de fin de pas serait refusé : mfem-mgis
  pose aussi celles de début de pas (mêmes valeurs à `setup`), ou MGIS
  accepte leur absence quand la loi ne les lit pas (à décider, voir
  [Questions](#questions-à-trancher)).

Tests : suite complète avec une bibliothèque `generic` (identique au
chemin actuel), puis avec `stlpar` et `openmp` (identique à l'arrondi
près) ; loi orthotrope (rotations) ; plusieurs intégrateurs sur un
attribut ; test d'un point en échec.

Fin : mêmes résultats, gain mesuré sur la poutre et le polycristal.

Réalisé (09/10/2026), écarts au plan :

- paramètre provisoire `BehaviourIntegration` (`PointByPoint`, `Batch`),
  valeur par défaut donnée par la variable d'environnement
  `MFEM_MGIS_BEHAVIOUR_INTEGRATION` ; il rejoindra `AssemblyStrategy` ;
- propriétés matériau de début de pas (question 6) : mfem-mgis donne à la
  loi, pour l'appel, des vues sur celles de fin de pas, comme le chemin
  point par point ; MGIS n'est pas modifié ;
- forces thermodynamiques : toujours tournées point par point par le
  résidu quand la loi fournit la rotation d'un point ; pour une loi
  orthotrope `generic-parallel` (rotations de lots seulement), tournées
  pour tous les points après l'intégration, dans un tableau à part (les
  forces du matériau restent dans son repère). Un cache invalidé par des
  versions des tableaux a été essayé puis retiré après la relecture : les
  écritures de MGIS par ses vues lui échappent ;
- échecs : MGIS s'arrête au premier point en échec pour une loi `generic`
  et intègre tous les points d'une loi `generic-parallel` ;
- le gain dépend du cas : négligeable sur la poutre (loi peu coûteuse),
  important sur le polycristal (voir les mesures).

## L3. Noyaux par lots et opérateur propre

Objectif : gradients, résidu, action de la tangente et matrices
élémentaires calculés par des noyaux par lots, sans passe élément par
élément sur l'hôte.

Noyaux (`src/Batch/`, `mfem::forall`, petites déformations 3D d'abord) :

- **géométrie** : facteurs géométriques (`Mesh::GetGeometricFactors`,
  jacobiens et déterminants) et fonctions de base (`DofToQuad`) pour la
  règle d'intégration de chaque matériau. Par défaut, la même règle que
  le chemin actuel (`2 * OrderGrad`). Les noyaux de MFEM demandent une
  règle par type de géométrie : un maillage mixte est traité géométrie
  par géométrie. Fait : pour un maillage sans nœuds, les jacobiens sont
  calculés par un noyau à partir des sommets (les facteurs géométriques de
  MFEM demandaient un espace d'éléments finis d'ordre un et une
  projection des coordonnées sur l'hôte, 1,8 s pour 1,6 million de
  tétraèdres) ;
- **gradients** : valeurs élémentaires (`ElementRestriction`), dérivées
  aux points, conversion vers les gradients de la loi (déformation en
  notation de Voigt avec √2, comme les matrices `B` actuelles), gradients
  macroscopiques ajoutés, écriture dans les gradients de fin de pas du
  matériau, dans son rangement, pour les seuls éléments du matériau ;
- **résidu** : `Bᵀ σ w det J` par point (avec la rotation des forces dans
  une copie pour les lois orthotropes), somme par élément, puis
  transposée de la restriction ;
- **action de la tangente** : `Σ Bᵀ K B v` ;
- **diagonale** et **blocs nodaux 3×3** de la tangente (lisseurs) ;
- **matrices élémentaires** `Ke = Σ Bᵀ K B w det J`.

Les matrices `B` sont aujourd'hui générées par GiNaC. Deux variantes à
comparer : noyaux écrits à la main, génériques sur l'hypothèse, ou
générés par une extension de `code-generation/behaviour-integrator.cxx`.

Points d'entrée :

- **B, opérateur propre** (`AssemblyStrategy=BatchedOperator`, en
  premier). `Mult` et `GetGradient` de
  `NonLinearEvolutionProblemImplementation<p>` redéfinis. L'opérateur
  refait ce que fait aujourd'hui `(Par)NonlinearForm` :
  - prolongation `P` et `Pᵀ` (sur GPU,
    `DeviceConformingProlongationOperator`) ;
  - degrés de liberté imposés du formulaire (`SetEssentialTrueDofs`,
    périodique compris) : résidu nul, tangente éliminée ou
    `ConstrainedOperator`, diagonale à 1 pour les lisseurs ;
  - **toutes les autres listes d'intégrateurs du formulaire** (domaine
    avec marqueurs, bord, faces intérieures, faces de bord : pression,
    source de chaleur, Robin d'ex8, intégrateurs de l'utilisateur), sauf
    l'intégrateur de mfem-mgis : appelées sur le chemin actuel pour le
    résidu et la tangente assemblée ; pour l'option a, action sans
    matrice de ces termes, ou assemblage à part ; refus avec un message
    sinon, jamais d'omission silencieuse ;
  - l'intégrateur de domaine actuel n'est plus ajouté au formulaire
    (sinon le domaine serait compté deux fois) ;
  - si `Mult` intègre la loi, l'action « nouvelle estimation » du solveur
    de mfem-mgis ne l'intègre plus (sinon deux intégrations par
    itération).

  `Mult` intègre la loi puis assemble le résidu (option), ce qui rend
  possibles les solveurs non linéaires extérieurs ; un échec de la loi
  est signalé sans arrêter le programme. `mgis_integrator` (type `final`)
  devient une interface commune aux deux chemins, ou le chemin par lots a
  son propre objet. Refus explicites sur le chemin par lots : attribut
  d'élément sans intégrateur (sinon contribution nulle et matrice
  singulière), mode `UseMultiMaterialNonLinearIntegrator=false`,
  maillages NURBS et collections autres que H1, intégrateurs autres que
  les 18 intégrateurs standard du générateur tant qu'ils n'ont pas de
  noyau.
- **A, assemblage partiel de MFEM** (`AssemblyStrategy=PartialAssembly`) :
  un `NonlinearFormIntegrator` qui implémente `AssemblePA`, `AddMultPA`,
  `AssembleGradPA`, `AddMultGradPA`, `AssembleGradDiagonalPA`, avec les
  vecteurs élémentaires en ordre lexicographique. Refus explicite si le
  problème a des intégrateurs de bord ou à marqueurs.

Tests :

- par élément : résidu et matrice élémentaires des noyaux contre
  `AssembleElementVector` et `AssembleElementGrad`, déplacements
  aléatoires, maillages hexaédriques et tétraédriques, ordres 1 et 2,
  plusieurs matériaux, numérotation byNODES et byVDIM, loi orthotrope ;
- global : résidu, action de la tangente, diagonale contre le chemin
  actuel, en séquentiel et sur 2 processus ;
- problème : `UniaxialTensileTest`, `TwoLayerCube`, `PredictionTest`,
  `ImposedPressureTest`, `ImposedDirichletBoundaryConditionAtClosestNode`,
  un problème périodique, avec les deux points d'entrée (le second doit
  refuser la pression) : mêmes itérations, mêmes résultats ; un test de
  chaque refus explicite.

Fin : suite des petites déformations 3D verte sur les deux points
d'entrée et les deux rangements.

## L4. Opérateur tangent et prédiction

Objectif : les options de l'opérateur tangent (voir
[Paramètres proposés](#paramètres-proposés)), le choix de la règle
d'intégration, et une prédiction qui n'utilise plus de formes bilinéaires
élément par élément.

Travail :

- `TangentOperator=Assembled` (c) : matrices élémentaires des noyaux,
  puis `HypreParMatrix`. Sur CPU d'abord. Motif creux, matrice hypre et
  schéma de communication gardés d'une itération à l'autre (seules les
  valeurs changent). Fait (11/10) : positions des entrées des matrices
  élémentaires dans les lignes de la matrice creuse calculées une fois,
  sans GPU (`computeStiffnessMatricesPositions`, assemblage ×1,1 à ×2,9 sur
  un cœur, sans gain sur le GPU), et produits `K B` et `Bᵀ K B` écrits sans
  former `B` (table constexpr des entrées non nulles : 254 puis 64
  registres en grandes transformations 3D, assemblage ×1,1 à ×3,1 sur un
  cœur et ×1,4 à ×1,9 sur le GPU, résultats identiques). Inverses des
  jacobiennes calculées une fois par le noyau de géométrie : assemblage
  sur le GPU ×1,7 à ×5,9. Matrice globale avec motif et échanges gardés
  (`batch::GlobalStiffnessMatrix`, espaces conformes, hypre sur l'hôte) :
  produit `Pᵀ A P` de 15 à 37 ms remplacé par 2 ms sur 8 processus.
  Restent la prédiction (même classe) et le cas de hypre sur le GPU.
- `TangentOperator=MatrixFree` (a) : action des noyaux, diagonale et
  blocs nodaux pour les lisseurs, AMG construit une fois sur l'opérateur
  élastique assemblé.
- `TangentOperator=ElementMatrices` : matrices élémentaires gardées,
  action élément par élément.
- `TangentOperator=ElasticOperator` (b) : avec `FixedPoint` (L5).
- `IntegrationRuleOrder` : règle choisie par l'utilisateur, par exemple
  8 points sur l'hexaèdre d'ordre 1 au lieu de 27 ; vérification de
  l'absence de modes parasites sur les cas de référence. Fait (11/10) :
  paramètre de la discrétisation, même ordre pour toutes les géométries,
  règles à poids non positif refusées, espaces de quadrature et
  intégrations par la même règle (`applyIntegrationRuleOrder`).
- Prédiction (`...Implementation.cxx:98-233`) réécrite avec les noyaux :
  opérateur élastique, sécant ou tangent par l'intégration de prédiction
  de MGIS, second membre par le noyau de résidu plus la contribution des
  conditions aux limites (`addLinearFormIntegrators`), respect de
  `DiscardLinearSolverFailure`, et les quatre axes actuels conservés.
  Nouvelle stratégie : extrapolation de l'incrément du pas précédent.
  Paramètre `Prediction`. `getLinearizedOperators` (API publique) : rendu
  par les noyaux, ou refusé sur le chemin par lots avec un message.
- Séquentiel : l'option c donne une `SparseMatrix` (UMFPACK, Krylov de
  MFEM), le parallèle une `HypreParMatrix`.
- Solveur ou préconditionneur fourni par l'utilisateur qui demande une
  matrice avec `MatrixFree` ou `ElasticOperator` : refus explicite.
- `IntegrationRuleOrder` change la taille des espaces de quadrature, les
  fonctions aux points de l'utilisateur, les liens entre problèmes
  (`resolveBehaviourIntegratorsDependencies`) et les références des tests :
  les intégrales (`computeMeasure`, `computeIntegral`, énergies,
  résultante) suivent la règle choisie ; deux problèmes liés gardent la
  même règle.

Tests : chaque option sur les cas de référence (mêmes solutions à la
tolérance de Newton près) ; un test par option qui échoue si elle est
ignorée (nombre d'itérations, type de l'opérateur rendu, nombre de points
d'intégration) ; prédictions contre `PredictionTest`, `PredictionTest2` et
`ImposedPressureTest2` (seul test qui combine prédiction et pression) ;
option c en séquentiel avec UMFPACK.

Fin : les options et les prédictions sur les cas de référence.

## L5. Solveurs

Objectif : brancher les solveurs de [solveurs-gpu.md](solveurs-gpu.md),
par vagues, dans les fabriques de mfem-mgis.

**Vague 1** (nécessaire aux options a et b et au GPU) :

- `PreconditionerUpdate` dans le Newton de mfem-mgis
  (`NewtonSolver.cxx:245`) : le solveur linéaire reçoit toujours la
  nouvelle tangente, mais son préconditionneur est enveloppé dans un
  objet dont `SetOperator` ne reconstruit rien quand il ne le faut pas
  (avec les Krylov de hypre, la phase de construction du préconditionneur
  doit aussi être neutralisée) ; enveloppe `ElasticOperator` ;
- solveur `FixedPoint` avec accélération d'Anderson ;
- Newton inexact (Eisenstat-Walker, fait : option
  `AdaptiveLinearSolverTolerance` du Newton et de `NewtonK0`) et recherche
  linéaire, écrits dans le
  Newton de mfem-mgis (les crochets de `mfem::NewtonSolver` ne s'y
  appliquent pas) ;
- lisseurs sans matrice : `OperatorJacobi`, `OperatorChebyshev`,
  `NodalBlockJacobi` ;
- options de BoomerAMG exposées, sans être écrasées par
  `SystemHypreBoomerAMG::SetOperator` ; stratégie `Elasticity` refusée sur
  GPU (MFEM arrête le programme en byVDIM), `System` proposée à la place ;
- `hasConverged` et `getNumberOfIterationsAtConvergence` étendus à tous
  les solveurs ajoutés (aujourd'hui tout solveur inconnu est déclaré
  convergé : un échec passerait en silence), avec un test qui fait
  échouer chaque solveur ; `DiscardLinearSolverFailure` et les
  informations d'itérations gardés dans les solveurs nouveaux ;
- point fixe et méthodes voisines : voir la section suivante.

**Vague 2** (comparaisons sur GPU) :

- L-BFGS préconditionné, à écrire (`LBFGSSolver` de MFEM ignore le
  préconditionneur ; `SNESQN` de PETSc comme référence) ;
- Ginkgo dans la fabrique : solveurs et préconditionneurs distribués de
  MFEM 4.10 (Schwarz), préconditionneurs locaux (Jacobi par blocs à
  précision adaptative, ILU, ISAI) ; cuDSS par `CuDSSSolver` (réutilisation
  de l'analyse symbolique ; fait le 12/10, solveur `CuDSSSolver` de la
  fabrique, ex8 sur le GPU 9,8 -> 4,4 s) ;
- AmgX par `AmgXSolver`, avec `UpdateOperator` pour la réutilisation
  complète ;
- KINSOL (`KINSolver`), avec l'opérateur qui intègre la loi (L3) ;
- réparation du chemin PETSc (même opérateur, plus d'arrêt sur échec de
  la loi), puis GAMG avec espace quasi nul ;
- précision mixte (Ginkgo, puis AmgX `dDFI` par une enveloppe propre).

**Vague 3** (selon les mesures) : p-multigrille et LOR (ordre 2 et plus),
SuperLU_DIST et STRUMPACK sur GPU, MGR de hypre pour les
problèmes couplés, FSAI et ILU de hypre, Krylov pipeliné, multigrille
de Ginkgo (par son API), réutilisation partielle de la hiérarchie AMG
(GAMG, AmgX par une enveloppe propre).

Tests : chaque solveur sur un cas où il converge, contre une solution de
référence ; chaque option avec un test qui échoue si elle est ignorée
(nombre de constructions du préconditionneur pour `PreconditionerUpdate`,
nombre d'itérations pour Anderson, réglages de BoomerAMG relus après
plusieurs itérations...).

Avec un schéma de couplage itératif (plusieurs résolutions dans un même
pas), la politique de réutilisation du préconditionneur vaut aussi entre
itérations de couplage.

Fin : vague 1 complète et mesurée ; vagues suivantes au fil de L9.

### Point fixe et méthodes voisines

Étude préalable (04/10/2026, hors de ce dépôt, prototype branché dans la
fabrique de solveurs de mfem-mgis) : sur une poutre en plasticité J2
(112 000 inconnues), le point fixe sur l'opérateur élastique `K0`,
accéléré par Anderson, avec un seul cycle d'AMG gelé par itération, est
2,15 fois plus rapide que le meilleur Newton sur CPU, et environ 4 fois
sur GPU (2 fois plus rapide qu'un Newton à AMG de `K0` gelé). Sur un
polycristal (loi chère, écrouissage faible), il est 3,3 fois plus lent
sur CPU : 8 fois plus d'itérations, donc d'intégrations de la loi. Le GPU
déplace exactement ce rapport. D'où toutes les variantes ci-dessous,
comparées sur les mêmes cas.

**Méthodes** (solveur non linéaire `FixedPoint`, paramètre `Method`) :

- `ElasticOperator` : point fixe sur `K0` (rigidité initiale, comme
  Cast3M). `K0` est assemblé une fois, pour tout le calcul ou une fois par
  pas (`ElasticOperatorUpdate`). Chaque itération : intégration de la loi,
  résidu, une résolution approchée avec `K0`, accélération ;
- `ModifiedNewton` : tangente assemblée et préconditionnée (ou factorisée)
  à la première itération du pas seulement ;
- `NewtonK0` (fait comme méthode du solveur `FixedPoint`) : tangente
  cohérente, Krylov préconditionné par l'AMG de `K0` construit une fois.
  Mêmes itérations que Newton, sans construction d'AMG par itération ;
  c'est le vrai concurrent du point fixe sur GPU. Le Krylov est un
  solveur itératif de MFEM de la fabrique (`CGSolver` par défaut,
  `GMRESSolver`, `FGMRESSolver`, `BiCGSTABSolver`, `MINRESSolver`,
  `SLISolver`) : le gradient conjugué suppose une tangente symétrique
  définie positive.

**Résolution avec `K0`** (paramètre `LinearSolver`) :

- exacte : factorisation unique (MUMPS ou UMFPACK sur CPU, cuDSS sur GPU) ;
- inexacte : `n` cycles d'AMG de `K0` gelé par itération (un seul suffit
  avec Anderson), ou un Krylov court préconditionné par cet AMG.

**Accélérations** (paramètre `Acceleration`), toutes portées des
accélérateurs de MTest (TFEL, `mtest/src/*AccelerationAlgorithm.cxx`,
cadre unifié publié dans Computers & Mathematics with Applications
70(9), 2015) vers les vecteurs de MFEM : produits scalaires et combinaisons de
vecteurs, donc sur GPU sans code propre :

- `None` (référence : en J2 à écrouissage linéaire `H`, le taux est
  `3μ/(3μ+H)` dans les zones plastiques, prévisible) ;
- `Anderson` (sur les itérés `UAnderson` et sur les résidus `FAnderson` ;
  profondeur, amortissement, déclenchement, période) ;
- `Cast3M` ;
- sécantes : `Secant`, `CrossedSecant`, `AlternateSecant` ;
- `IronsTuck`, `Steffensen` ;
- delta carrés : `CrossedDelta2`, `Crossed2Delta`, `Crossed2Deltabis`,
  `AlternateDelta2`, `Alternate2Delta`.

Le prototype n'avait porté que Cast3M et Anderson. La revue
bibliographique de l'étude montre qu'Anderson n'est pas toujours le
meilleur : d'où le portage de toute la famille.

**Précision** (paramètre `OperatorPrecision`) : `K0` ou son AMG en simple
précision, résidu et accélération en double. Sur la RTX 5080, l'étude
mesure environ 30 % de gain à nombre d'itérations égal.

**Prédiction** : la prédiction élastique (indispensable en plasticité)
réutilise `K0` déjà factorisé ou son AMG : plus de factorisation par pas
(elle coûtait 35 s sur 44 avec MUMPS dans l'étude).

**Mesures propres à ces méthodes** : temps par phase (intégration, résidu,
assemblage, construction du préconditionneur, résolutions,
accélération), comme dans le prototype ; itérations ; et la comparaison
« meilleur contre meilleur » (jamais contre un Newton à solveur direct
seul, qui fausse les gains).

**Tests** : taux de convergence du point fixe non accéléré égal à
`3μ/(3μ+H)` sur un cas J2 à écrouissage linéaire (test exact) ; chaque
accélérateur sur un problème linéaire (Anderson y équivaut à GMRES
préconditionné par `K0`) ; mêmes solutions que Newton à la tolérance près
sur les cas de référence ; un test par paramètre qui échoue s'il est
ignoré (nombre de constructions de l'AMG, nombre de factorisations).

## L6. GPU

Objectif : tout le calcul d'une itération sur le GPU, données comprises.

Tranche minimale (dès la fin de L3) : petites déformations 3D, une loi J2
compilée pour `cuda`, option c puis a.

Travail :

- périphérique `cuda` de MFEM, tableaux des matériaux sur le GPU (L1),
  noyaux de L3 sur le GPU (fichiers compilés par nvcc), intégration par
  la bibliothèque `cuda` (ou `kokkos`) de la loi ;
- refus d'une bibliothèque `generic` sur un périphérique GPU, sauf
  configuration mixte explicite (matériau en mémoire partagée) ;
- conditions de Dirichlet appliquées sur le GPU
  (`Vector::SetSubVector` avec une valeur s'exécute sur le périphérique
  depuis MFEM 4.9) ; plus d'écriture hôte dans `u` pour les conditions de
  mfem-mgis. Conditions écrites par l'utilisateur (interface publique qui
  écrit `u` sur l'hôte) : `u` est synchronisé vers l'hôte avant l'appel et
  vers le périphérique après (correct, mais coûteux), et une interface
  pour le périphérique est proposée ;
- solveurs qui ne marchent que sur l'hôte (UMFPACK, MUMPS) refusés avec un
  périphérique GPU, sauf configuration mixte explicite (copie de la
  matrice) ;
- fin de pas sur le GPU par les opérations CUDA de MGIS (ou la variante
  `MFEMFieldOperations`) ;
- option c sur GPU : matrice globale assemblée sur le GPU à partir des
  matrices élémentaires. Le remplissage de matrice creuse de MFEM sur le
  périphérique (`ElementRestriction::FillJAndData`) travaille sur des
  degrés de liberté scalaires : pour 3 composantes, à vérifier, ou à faire
  par blocs de composantes ou par l'IJ de hypre sur GPU. Élimination des
  degrés de liberté imposés sur le périphérique ;
- plusieurs matériaux : appels asynchrones de MGIS puis une seule
  synchronisation ;
- vérification qu'aucune copie cachée n'a lieu : GPU sans HMM et profil
  Nsight Systems sur le vrai GPU ; politique de validité des copies
  contrôlée par le périphérique `debug` de MFEM, qui exécute les noyaux
  sur l'hôte (avec une loi `openmp` ou `stlpar`). En plusieurs
  processus, sans MPI compatible GPU, les échanges de la prolongation et
  de hypre passent par l'hôte : le critère ne vaut qu'en un processus.

Tests : les tests de L3 et L4 sur le périphérique `cuda`, comparés au CPU
à l'arrondi près ; les mêmes sous le périphérique `debug` avec une loi
`openmp`.

Fin : mêmes résultats qu'en CPU sur les cas de référence, aucune copie
cachée en un processus, premières mesures.

## L7. Tout mfem-mgis

Objectif : le chemin par lots couvre ce que couvre le chemin actuel.

Ordre proposé :

1. déformations planes et contraintes planes ;
2. grandes déformations (gradient de la transformation, contrainte de
   Piola-Kirchhoff, conversions de MGIS sur vues) ;
3. thermique non linéaire stationnaire et transitoire (la transitoire
   demande un noyau sur les valeurs des fonctions de forme, pas seulement
   sur leurs dérivées ; la température lue depuis les inconnues comme
   variable externe doit vivre sur le GPU). Stationnaire fait (gradient de
   température, flux, bloc dérivée du flux par rapport à la température,
   valeurs des fonctions de forme, température aux points recopiée dans la
   variable externe sur l'hôte). Transitoire faite (gradient = température
   au point, opérateur B = valeur de la fonction de forme, intégrateur
   `TransientHeatTransfer` déclaré aussi en 2D) ; reste la variable externe
   gardée sur le GPU ;
4. lois orthotropes : matrices de rotation par point en tableau, rotations
   de MGIS sur vues (plus de pointeurs de fonction) ;
5. propriétés matériau et variables externes sous toutes leurs formes,
   dont les propriétés `generic-parallel` évaluées par MGIS sur le GPU ;
   évaluateurs non uniformes copiés vers le GPU à chaque `setup` après
   synchronisation de leur source ;
6. conditions aux limites de bord (pression, source de chaleur) par des
   noyaux de bord, pour le GPU. Fait autrement (11/10) : ces conditions ne
   dépendent pas des inconnues, leur contribution est assemblée une fois
   pour une valeur unité, puis multipliée par la valeur à la fin du pas et
   ajoutée sur le GPU aux forces intérieures calculées par noyaux
   (`AbstractBoundaryCondition::isIndependentOfTheUnknowns`). La forme non
   linéaire ne contient plus que l'intégrateur de mfem-mgis : plus de
   résidu assemblé sur l'hôte, tangente assemblée par noyaux sans
   l'assemblage de MFEM, tangente sans assemblage possible avec une
   pression. Restent les conditions qui dépendent des inconnues (Robin
   d'ex8), qui demanderaient de vrais noyaux de bord ;
7. modèles et schémas de couplage (`Simulation`, `IterativeCouplingScheme`,
   `LoopCouplingScheme`, `NonLinearModel`, `PointWiseModel`,
   `QPEvaluators`), dépendances entre intégrateurs de deux problèmes
   (`resolveBehaviourIntegratorsDependencies`) : l'exemple
   thermomécanique ex8 en dépend, et les échanges de champs aux points
   entre modèles se font sur l'hôte.
   `PointWiseModel` (loi aux points sans éléments finis) est un premier
   cas simple de loi par lots sur GPU ;
8. micromorphe ;
9. régularisations (FBar, Faltus) : elles demandent la solution courante
   à l'assemblage, à traiter en dernier. FBar fait (11/10) : gradient au
   centre des éléments par une règle d'ordre 0, `F̄ = α F`, tangente
   standard modifiée plus un terme de rang un par élément (voir
   `batch::FBarData`), en 3D et en déformations planes, matrice assemblée,
   diagonale et action sans matrice. Faltus fait (11/10) : gradient au
   centre de la même règle, résidu `Σ w Dᵀ α (F - F0)`, tangente
   `α Σ w Dᵀ D` constante calculée à la volée par paires de nœuds (bloc
   `d_a · d_b` fois l'identité, `d` différence des dérivées des fonctions
   de forme au point et au centre), en 3D et en déformations planes
   (contraintes planes refusées) ;
10. problèmes périodiques.

Tests : la suite de mfem-mgis passée sur le chemin par lots, cas par cas.
La suite actuelle ne teste ni FBar, ni Faltus, ni les contraintes planes,
ni la thermique transitoire, ni la source de chaleur, ni le couplage
itératif, ni PETSc : ces tests sont à écrire dans ce lot (sur le chemin
actuel d'abord, comme référence), ou les exemples correspondants de
mfem-mgis-examples entrent dans la validation.

Fin : suite complète, complétée de ces tests, sur le chemin par lots, sur
CPU et sur GPU.

## L8. Sorties et post-traitements

- D'abord : copie vers l'hôte avant chaque post-traitement (`HostRead`),
  comportement actuel inchangé. Avec le rangement composante par
  composante, les accès des post-traitements sont à adapter ou les
  données à transposer vers une copie point par point.
- Ensuite, selon les mesures : forces thermodynamiques moyennes, résultante
  sur un bord, énergies et intégrales (`computeMeasure`,
  `computeIntegral`) par des réductions sur le GPU, avec la règle
  d'intégration choisie ; projections L2 et régularisation par gradient
  implicite (`L2Projection.cxx:624-702`) par l'assemblage partiel de
  MFEM.

Tests : sorties identiques sur CPU et sur GPU, dans les deux rangements.

## L9. Exemples et mesures

- mfem-mgis-examples et mm-opera-hpc : options `--device`,
  `AssemblyStrategy`, `TangentOperator`, `QuadratureDataLayout`,
  `IntegrationRuleOrder`, solveur, passées par la ligne de commande.
- Cas particuliers de ces dépôts à traiter : la condition de Robin d'ex8
  (intégrateur de faces de bord dont la tangente dépend de `u`) et son
  couplage thermomécanique ; les propriétés matériau posées dans la
  mémoire hôte de l'utilisateur (`cermet`) ; les polycristaux périodiques
  (degrés de liberté imposés par le problème périodique) ; la bulle
  (`bubble`).
- Campagne de mesures (section suivante), résultats dans `results/`,
  scripts dans `scripts/`.

## L10. Documentation et propositions amont

- Documentation de mfem-mgis (en anglais, dans son arbre) : chemin par
  lots, paramètres, configurations testées.
- Découpage en propositions amont indépendantes : corrections (L0),
  données aux points (L1), lois par lots (L2), noyaux et opérateur (L3),
  solveurs (L5), GPU (L6). Voir [integration-amont.md](integration-amont.md).

## Stratégie de test

- Référence unique : le chemin actuel, avec une bibliothèque `generic`.
- Trois niveaux de comparaison : élément (résidu et matrice
  élémentaires), global (résidu, action de la tangente, diagonale),
  problème (itérations de Newton, déplacements, contraintes, variables
  internes). Tolérance relative 1e-12 aux deux premiers niveaux, à
  l'arrondi près au troisième.
- Chaque option nouvelle a un test qui échoue si l'option est ignorée.
- GPU : mêmes tests sur le périphérique `cuda`, comparés au CPU ; mêmes
  tests sous le périphérique `debug` de MFEM avec une loi `openmp`, pour
  la validité des copies.
- Passe ciblée (une vingtaine de minutes) à chaque commit, suite complète
  à la fin d'un lot, relecture indépendante avant de déclarer un lot fini.

## Campagne de mesures

Cas :

1. poutre en plasticité J2, 3D, petites déformations (cas des essais
   préalables) ;
2. ssna303 (mfem-mgis-examples ex2, grandes déformations) ;
3. polycristal (mm-opera-hpc, plusieurs matériaux, loi de monocristal) ;
4. thermique non linéaire ;
5. cube raffiné jusqu'à remplir la mémoire du GPU (passage à l'échelle).

Mesures : temps par pas et par phase (gradients, loi, résidu, tangente,
construction du préconditionneur, résolution, fin de pas, sorties),
itérations de Newton et de Krylov, mémoire du GPU, transferts entre hôte
et GPU (dont ceux des échanges MPI sans MPI compatible GPU).

Comparaisons :

- CPU (1 à 16 cœurs, MPI) contre GPU ;
- chemin actuel contre chemin par lots sur CPU ;
- rangement point par point contre composante par composante ;
- options a, b, c et éléments de l'opérateur tangent ;
- règle d'intégration (27 points contre 8 sur l'hexaèdre d'ordre 1) ;
- solveurs, vague par vague ;
- bibliothèque `cuda` contre `kokkos` de la même loi, et `generic` sur
  l'hôte.

## Questions à trancher

1. **Noms des paramètres** proposés plus haut.
2. **Ordre des points d'entrée** : opérateur propre d'abord (proposé), à
   cause des intégrateurs de bord ignorés ou refusés par l'assemblage
   partiel de MFEM.
3. **MFEM 4.10** : tranché le 09/10/2026, tout le prototype passe en 4.10
   (voir [solveurs-gpu.md](solveurs-gpu.md#mfem-410)).
4. **Noyaux écrits à la main ou générés** (extension du générateur GiNaC
   de mfem-mgis).
5. **Propriétaire de la mémoire** : mfem-mgis (proposé) ou MGIS avec une
   ressource mémoire ; les deux sont prévus en L1, à comparer.
6. **Propriétés matériau de début de pas** (voir L2) : tranché en L2,
   mfem-mgis donne des vues sur celles de fin de pas, comme le chemin
   point par point.
7. **Base de travail** : `a011d97d6f`, ou la refonte de
   `PartialQuadratureSpace` en cours en amont une fois fusionnée.
8. **Propositions amont** : à quel moment et dans quel ordre.

## Risques

- **Taille** : beaucoup d'intégrateurs et de cas particuliers. Le chemin
  par lots commence par les cas les plus utilisés et reste optionnel.
- **Branches amont en cours** : la refonte de `PartialQuadratureSpace`
  (branche `466-refactor-partial-quadrature-space-to-better-support-spaces-on-boundaries`,
  89 fichiers) touche les classes que L1 modifie ; la branche
  `452-post-processings-only-consider-the-first-behaviour-integrator-of-a-material`
  touche les post-traitements de L1 et L8. Conflits à prévoir.
- **Assemblage partiel de MFEM** : limites décrites plus haut ; le
  point d'entrée B ne dépend pas de l'extension.
- **Compilation CUDA dans mfem-mgis** : les noyaux demandent nvcc ;
  temps de compilation et mémoire à surveiller, comme pour MFEM.
- **Coût de l'AMG sur GPU** : construction plus chère que la résolution ;
  d'où la réutilisation du préconditionneur en vague 1.
- **Double précision sur GPU grand public** : la RTX 5080 calcule la
  double précision 64 fois plus lentement que la simple ; les mesures
  doivent être confirmées sur un GPU de calcul (A100 ou H100).
- **MFEM 4.10** : mfem-mgis n'a été utilisé qu'avec 4.8 et 4.9 ;
  compatibilité vérifiée en L0. Compilation avec CUDA 13 : `ptxas` prend
  22 Go pour un fichier de MFEM (constructions encadrées par
  `scripts/memory-governor.sh`).
- **Échecs locaux de l'intégration** : sur GPU, un point en échec doit
  remonter sans copie de tout l'état ; MGIS rend le premier point en échec
  et un statut global.
- **Plusieurs processus sans MPI compatible GPU** : copies vers l'hôte à
  chaque échange, inévitables avec l'environnement local.

- **Plusieurs GPU par nœud et deux runtimes CUDA** : les noyaux de
  mfem-mgis passent par le runtime CUDA partagé, ceux de MFEM par la copie
  statique du runtime que contient `libmfem.so`. Les deux suivent le
  contexte courant du pilote, donc a priori le même GPU, mais ce n'est pas
  vérifiable sur une machine à un seul GPU. À tester sur un nœud à
  plusieurs GPU (un GPU par processus, opérations de MFEM et de mfem-mgis
  sur les mêmes tableaux).
- **Stockage externe entre matériaux** : un matériau peut lire ses
  variables externes dans les tableaux d'un autre (mode `EXTERNAL_STORAGE`,
  exemple ex8 : gonflement lu dans les variables internes du modèle de
  gonflement). La synchronisation d'un matériau ne couvre pas les tableaux
  de l'autre : à traiter quand ces tableaux seront écrits sur le GPU.
