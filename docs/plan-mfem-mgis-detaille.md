# Plan détaillé pour mfem-mgis

Statut : plan établi le 09/10/2026, rien n'est commencé dans mfem-mgis.
Les choix de fond sont dans [plan-mfem-mgis.md](plan-mfem-mgis.md) ; la
liste des solveurs dans [solveurs-gpu.md](solveurs-gpu.md). Ce document
découpe le travail en lots, dit quels fichiers changent, ce qui est testé
et quand un lot est fini.

Base : mfem-mgis master `a011d97d6f` (06/10/2026), dépôt `mfem-mgis/`.
Les références `fichier:ligne` renvoient à cette base. TFEL et MGIS sont
ceux du prototype (MGIS installé dans `install/mgis` par `build-mgis.sh`),
MFEM 4.9 avec CUDA vient de l'environnement `mmm-gpu`.

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
- L'intégration est déclenchée par une action du solveur non linéaire
  (`NonLinearSolverBase.cxx:39-58`), le résidu par
  `computeResidual` (`:135-140`), la tangente par `getJacobian`
  (`:142-152`). Ces deux fonctions passent par `oper->Mult` et
  `oper->GetGradient` : un opérateur propre peut s'y brancher.
- Le préconditionneur est reconstruit à chaque itération :
  `prec->SetOperator(getJacobian(u))` (`NewtonSolver.cxx:245`).
- Le seul solveur non linéaire est `Newton` (`NonLinearSolverFactory.cxx`),
  plus `PetscNonlinearSolver` quand PETSc est demandé.
- La prédiction (`...Implementation.cxx:98-233`) passe par une
  `BilinearForm` et une `LinearForm` (`StiffnessMatrixIntegrator`,
  `InternalForcesIntegrator`), élément par élément. Elle se règle en C++
  seulement (`setPredictionPolicy`) : il n'y a pas de paramètre.
- Les conditions de Dirichlet écrivent les valeurs imposées dans `u` sur
  l'hôte (`UniformDirichletBoundaryCondition.cxx:64-79`).
- Fin de pas : `update` et `revert` de MGIS sur l'hôte, matériau par
  matériau (`BehaviourIntegratorBase.cxx:779-785`).

### Intégrateurs, matériaux, données aux points

- Un intégrateur de comportement **est** le matériau :
  `BehaviourIntegratorBase` dérive de `Material`, qui dérive de
  `mgis::behaviour::MaterialDataManager`.
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
- Rotations des lois orthotropes : pointeurs de fonction de MGIS
  (`rotate_gradients_ptr`...), matrice par point donnée par un `variant`
  (`Material.hxx:99-106`).
- `PartialQuadratureSpace` : points contigus par élément, éléments d'un
  matériau dans l'ordre croissant, mais décalages rangés dans une
  `std::unordered_map` (`PartialQuadratureSpace.hxx:230-236`).
- `PartialQuadratureFunction` : `std::vector` ou mémoire empruntée, vues
  `std::span`, accès `values[o * stride + begin]`.
- Propriétés matériau et variables externes : évaluateurs
  (`AbstractQPEvaluator`), évalués sur l'hôte dans `setup`
  (`BehaviourIntegratorBase.cxx:360-670`).

### Solveurs linéaires

- Fabrique `LinearSolverFactory.cxx` : Krylov de MFEM (séquentiel et
  parallèle), `HyprePCG`, `HypreGMRES`, `HypreFGMRES`, `MUMPSSolver`,
  `UMFPackSolver` ; préconditionneurs `HypreBoomerAMG` (stratégies
  `Elasticity`, `System`, `None`), `HypreEuclid`, `HypreILU`,
  `HypreParaSails`, `HypreDiagScale`. Tous les préconditionneurs et
  solveurs directs demandent une matrice assemblée.
- Deux défauts relevés au passage : `HypreGMRES` teste un paramètre
  `AbsoluteTolerance` qui n'est pas dans la liste autorisée (code mort),
  `HypreFGMRES` aussi, et appelle `SetTol` au lieu de `SetAbsTol`.

### GPU

- Aucun usage de `mfem::Device`, `AssemblyLevel`, `forall`,
  `ElementRestriction`, `QuadratureInterpolator`, `GeometricFactors` ou
  `DofToQuad`. `mfem_mgis::initialize` ne crée pas de périphérique.

### Ce que l'assemblage partiel de MFEM 4.9 impose

Vérifié dans les sources de MFEM 4.9 (`fem/nonlinearform.cpp`,
`fem/nonlinearform_ext.cpp`, `fem/pnonlinearform.cpp`) :

- `NonlinearForm::Mult` rend la main juste après l'extension : les
  intégrateurs de bord (`AddBoundaryIntegrator`) et les marqueurs des
  intégrateurs de domaine sont **ignorés sans message**. Dans mfem-mgis,
  cela touche la pression imposée (`UniformImposedPressureBoundaryCondition`)
  et la source de chaleur (`UniformHeatSourceBoundaryCondition`).
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
- **Travail dans `mfem-mgis/`, branche `gpu-proto`**, base `a011d97d6f`,
  copié dans ce dépôt par `scripts/sync-repository.sh`. Messages de commit
  en anglais. Un lot fini est un ensemble de commits relus.
- **Code nouveau regroupé** : `include/MFEMMGIS/Batch/` et `src/Batch/`
  (nom à discuter). Les classes existantes gagnent des points
  d'extension, pas de code GPU.
- **Noyaux écrits une fois**, avec `mfem::forall`, génériques sur
  l'hypothèse et le type de gradient. Ils servent sur CPU et sur GPU, et
  aux deux points d'entrée.
- **Toutes les variantes, choisies par option**, comparées sur les mêmes
  cas (rangement, point d'entrée, opérateur tangent, solveur).
- **Chaque option a un test qui échoue si l'option est ignorée.**
- **Pas d'échec silencieux** : une configuration non prise en charge
  (intégrateur de bord avec l'assemblage partiel, stratégie `Elasticity`
  sur GPU, loi `generic` sur un périphérique GPU...) est refusée avec un
  message.

## Paramètres proposés

Noms à discuter, dans le style des paramètres existants.

| Où | Paramètre | Valeurs |
|---|---|---|
| problème | `AssemblyStrategy` | `Legacy` (défaut, chemin actuel), `BatchedOperator` (opérateur propre), `PartialAssembly` (assemblage partiel de MFEM) |
| problème | `QuadratureDataLayout` | `PointByPoint` (défaut), `ComponentByComponent` |
| problème | `TangentOperator` | `Assembled` (option c, défaut), `MatrixFree` (option a), `ElementMatrices` (assemblage par élément), `ElasticOperator` (option b) |
| problème | `Prediction` | `None`, `ElasticOperator`, `SecantOperator`, `TangentOperator`, `Extrapolation` (aujourd'hui réglé en C++ seulement) |
| solveur non linéaire | `Newton` | options nouvelles : `AdaptiveLinearSolverTolerance` (Eisenstat-Walker), `LineSearch`, `PreconditionerUpdate` |
| solveur non linéaire | `FixedPoint` | `AndersonDepth`, `AndersonDamping`, `PreconditionerUpdate` |
| solveur non linéaire | `LBFGS`, `KINSOL` | options propres |
| préconditionneur | `OperatorJacobi`, `OperatorChebyshev`, `NodalBlockJacobi`, `LOR` | sans matrice |
| préconditionneur | `HypreBoomerAMG` | options nouvelles : grossissement, interpolation, lisseur, balayages, seuil, niveaux |
| préconditionneur | `ElasticOperator` | enveloppe : préconditionneur construit une fois sur l'opérateur élastique |
| solveur, préconditionneur | `Ginkgo*`, `AmgX`, `CuDSS` | selon [solveurs-gpu.md](solveurs-gpu.md) |
| application | `--device` | configuration de `mfem::Device` (`cpu`, `cuda`, `debug`...), ajoutée à `declareDefaultOptions` |

`PreconditionerUpdate` : `EachIteration` (comportement actuel),
`EachTimeStep`, `Never`, `EveryNIterations`.

## Vue d'ensemble des lots

| Lot | Contenu | Dépend de | Taille | Fin |
|---|---|---|---|---|
| L0 | construction, référence, périphérique | environnement `mmm-gpu` | S | tests existants verts avec MFEM CUDA, temps de référence |
| L1 | données aux points dans la mémoire de MFEM, deux rangements | L0 | M | suite verte dans les deux rangements |
| L2 | lois intégrées par lots sur CPU | L1 | M | mêmes résultats, gain mesuré avec `stlpar` et `openmp` |
| L3 | noyaux par lots, opérateur propre, assemblage partiel | L2 | L | résidu et tangente égaux au chemin actuel, suite verte |
| L4 | opérateur tangent (a, b, c, éléments), prédiction | L3 | M | les quatre options sur les cas de référence |
| L5 | solveurs, en trois vagues | L3, L4 | L | vague 1 complète |
| L6 | GPU | L3 (tranche minimale), puis L4, L5 | L | mêmes résultats qu'en CPU, aucune copie cachée |
| L7 | toutes les hypothèses et cinématiques | L3 à L6 | L | suite complète de mfem-mgis sur le chemin par lots |
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
  la racine (comme `build-mgis.sh`) : MFEM 4.9 avec CUDA de `mmm-gpu`,
  MGIS de `install/mgis`, TFEL de `install/tfel-mmm-gpu`, construction,
  `make check`.
- Vérifier que la suite actuelle passe avec MFEM compilé pour CUDA mais
  périphérique `cpu` (rien n'est attendu sur GPU à ce stade).
- Option `--device` : `mfem::Device` créé par `mfem_mgis::initialize` (ou
  par l'application) à partir d'une chaîne ; `Config.cxx:174-222`.
- Chronomètres par phase dans le chemin actuel (intégration, résidu,
  tangente, préconditionneur, résolution, fin de pas), affichés au niveau
  de verbosité adéquat. Ils servent de référence à toutes les mesures.
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
  `MaterialDataManagerInitializer` (`Material.cxx:48-61`). Les tableaux
  (gradients, forces thermodynamiques, variables internes, énergies des
  deux états, blocs de la tangente) sont des `mfem::Vector` possédés par
  le matériau ; leurs pointeurs sont donnés à MGIS avec l'emplacement
  (`HOST` ou `DEVICE`, d'après `mfem::Device`) et le rangement
  (`QuadratureDataLayout`). Variante à comparer : tableaux alloués par
  MGIS avec une ressource mémoire (`MGIS/MemoryResources`), vus par MFEM.
- `PartialQuadratureSpace` : ajouter, à côté de la table de hachage, des
  tableaux `mfem::Array<int>` des éléments du matériau et des décalages,
  lisibles sur le périphérique. La table reste pour le code actuel.
- `PartialQuadratureFunction` : décrire le rangement composante par
  composante (pas entre composantes) en plus du pas entre points ; les
  accès actuels (`getIntegrationPointValues`...) doivent donner les mêmes
  valeurs dans les deux rangements.
- Les fonctions qui lisent les données (post-traitements, projections,
  sorties, `computeMeanThermodynamicForcesValues`...) passent par ces
  accès : à vérifier une par une dans les deux rangements.
- Fin de pas : `update` et `revert` passent par les surcharges de MGIS qui
  prennent des `FieldOperations`. Ajouter `MFEMFieldOperations`, qui
  implémente `mgis::FieldOperations` avec `mfem::forall` et les copies de
  `mfem::Memory` : portable CPU, CUDA et HIP, sans compiler mfem-mgis avec
  nvcc.

Tests : toute la suite avec `QuadratureDataLayout=ComponentByComponent`
(résultats identiques bit à bit au rangement point par point, puisque le
chemin actuel calcule les mêmes valeurs) ; un test unitaire de
`MFEMFieldOperations` (copies, transpositions) ; un test de chaque accès
de `PartialQuadratureFunction` dans les deux rangements.

Fin : suite verte dans les deux rangements.

## L2. Lois intégrées par lots sur CPU

Objectif : remplacer l'intégration point par point par un appel par
matériau, sans toucher à l'assemblage. Premier gain mesurable : les
backends `stlpar` et `openmp` de TFEL sur CPU.

Travail :

- Découper `integrate` (`...Implementation.cxx:485-514`) :
  1. passe élément par élément, qui calcule seulement les gradients (code
     actuel de `implementIntegrate` sans l'appel de la loi), rotation
     comprise ;
  2. pour chaque matériau, un appel `mgis::behaviour::integrate` sur une
     vue du matériau (`make_view`), avec la politique d'exécution choisie ;
  3. rotation des forces et des blocs de la tangente, sur vues.
- Les bibliothèques `generic` restent acceptées : avec un espace de
  travail, MGIS les intègre point par point (référence). Les bibliothèques
  `generic-parallel` sont choisies par le chemin de la bibliothèque.
- Propriétés matériau et variables externes : les valeurs évaluées dans
  `setup` sont données à la vue (uniformes ou tableaux).
- Échecs : MGIS rend le premier point en échec et son statut ; le
  comportement actuel (recherche linéaire de Newton, sous-pas) doit être
  conservé.

Tests : suite complète avec une bibliothèque `generic` (identique au
chemin actuel), puis avec `stlpar` et `openmp` (identique à l'arrondi
près) ; test d'un point en échec.

Fin : mêmes résultats, gain mesuré sur la poutre et le polycristal.

## L3. Noyaux par lots et opérateur propre

Objectif : gradients, résidu, action de la tangente et matrices
élémentaires calculés par des noyaux par lots, sans passe élément par
élément sur l'hôte.

Noyaux (`src/Batch/`, `mfem::forall`, petites déformations 3D d'abord) :

- **géométrie** : facteurs géométriques (`Mesh::GetGeometricFactors`,
  jacobiens et déterminants) et fonctions de base (`DofToQuad`) pour la
  règle d'intégration de chaque matériau (`IntRules.Get(geom, 2 * order)`,
  la même que le chemin actuel) ;
- **gradients** : valeurs élémentaires (`ElementRestriction`), dérivées
  aux points, conversion vers les gradients de la loi (déformation en
  notation de Voigt avec √2, comme les matrices `B` actuelles), écriture
  dans les gradients de fin de pas du matériau, dans son rangement,
  pour les seuls éléments du matériau ;
- **résidu** : `Bᵀ σ w det J` par point, somme par élément, puis
  transposée de la restriction ;
- **action de la tangente** : `Σ Bᵀ K B v` ;
- **diagonale** et **blocs nodaux 3×3** de la tangente (lisseurs) ;
- **matrices élémentaires** `Ke = Σ Bᵀ K B w det J`.

Les matrices `B` sont aujourd'hui générées par GiNaC. Deux variantes à
comparer : noyaux écrits à la main, génériques sur l'hypothèse, ou
générés par une extension de `code-generation/behaviour-integrator.cxx`.

Points d'entrée :

- **B, opérateur propre** (`AssemblyStrategy=BatchedOperator`, en
  premier) : `Mult` et `GetGradient` de
  `NonLinearEvolutionProblemImplementation<p>` redéfinis. Le domaine passe
  par les noyaux ; les intégrateurs de bord restent sur le chemin actuel
  et sont ajoutés (le formulaire ne les ignore pas sur ce chemin).
  `mgis_integrator` (type `final`) devient une interface commune aux deux
  chemins, ou le chemin par lots a son propre objet.
- **A, assemblage partiel de MFEM** (`AssemblyStrategy=PartialAssembly`) :
  un `NonlinearFormIntegrator` qui implémente `AssemblePA`, `AddMultPA`,
  `AssembleGradPA`, `AddMultGradPA`, `AssembleGradDiagonalPA`, avec les
  vecteurs élémentaires en ordre lexicographique. Refus explicite si le
  problème a des intégrateurs de bord ou à marqueurs.

Tests :

- par élément : résidu et matrice élémentaires des noyaux contre
  `AssembleElementVector` et `AssembleElementGrad`, déplacements
  aléatoires, maillages hexaédriques et tétraédriques, ordres 1 et 2,
  plusieurs matériaux, numérotation byNODES et byVDIM ;
- global : résidu, action de la tangente, diagonale contre le chemin
  actuel, en séquentiel et sur 2 processus ;
- problème : `UniaxialTensileTest`, `TwoLayerCube`, `PredictionTest`
  avec les deux points d'entrée : mêmes itérations, mêmes résultats.

Fin : suite des petites déformations 3D verte sur les deux points
d'entrée et les deux rangements.

## L4. Opérateur tangent et prédiction

Objectif : les options de l'opérateur tangent du plan, et une prédiction
qui n'utilise plus de formes bilinéaires élément par élément.

Travail :

- `TangentOperator=Assembled` (c) : matrices élémentaires des noyaux,
  puis matrice `HypreParMatrix`. Sur CPU d'abord : assemblage de MFEM ou
  IJ de hypre. Sur GPU : voir L6.
- `TangentOperator=MatrixFree` (a) : action des noyaux, diagonale et
  blocs nodaux pour les lisseurs, AMG construit sur l'opérateur élastique
  (assemblé une fois).
- `TangentOperator=ElementMatrices` : matrices élémentaires gardées,
  action élément par élément.
- `TangentOperator=ElasticOperator` (b) : point fixe avec l'opérateur
  élastique (voir L5, `FixedPoint`).
- Prédiction (`...Implementation.cxx:98-233`) réécrite avec les noyaux :
  opérateur élastique, sécant ou tangent par l'intégration de prédiction
  de MGIS, second membre par le noyau de résidu. Nouvelle prédiction par
  extrapolation de l'incrément du pas précédent. Paramètre `Prediction`.

Tests : chaque option sur les cas de référence (mêmes solutions à la
tolérance de Newton près) ; un test par option qui échoue si elle est
ignorée (nombre d'itérations, type de l'opérateur rendu) ; prédictions
contre `PredictionTest` et `PredictionTest2`.

Fin : les quatre options et les prédictions sur les cas de référence.

## L5. Solveurs

Objectif : brancher les solveurs de [solveurs-gpu.md](solveurs-gpu.md),
par vagues, dans les fabriques de mfem-mgis.

**Vague 1** (nécessaire aux options a et b et au GPU) :

- `PreconditionerUpdate` dans `NewtonSolver`
  (`NewtonSolver.cxx:245`) ; enveloppe `ElasticOperator` ;
- solveur `FixedPoint` avec accélération d'Anderson ;
- Newton inexact (`SetAdaptiveLinRtol`) et recherche linéaire
  (`ComputeScalingFactor`) ;
- `L-BFGS` (`LBFGSSolver`) ;
- lisseurs sans matrice : `OperatorJacobi`, `OperatorChebyshev`,
  `NodalBlockJacobi` ;
- options de BoomerAMG exposées ; stratégie `Elasticity` refusée sur GPU
  (MFEM arrête le programme), `System` proposée à la place.

**Vague 2** (comparaisons sur GPU) :

- Ginkgo dans la fabrique (Krylov, Jacobi par blocs à précision
  adaptative, ILU, ISAI ; local à chaque processus en MFEM 4.9) ;
- AmgX par `AmgXSolver`, avec `UpdateOperator` pour la réutilisation ;
- KINSOL (`KINSolver`) ;
- réparation du chemin PETSc, puis GAMG avec espace quasi nul ;
- précision mixte (Ginkgo, puis AmgX `dDFI` par une enveloppe propre).

**Vague 3** (selon les mesures) : p-multigrille et LOR (ordre 2 et plus),
cuDSS, SuperLU_DIST et STRUMPACK sur GPU, MGR de hypre pour les
problèmes couplés, FSAI et ILU de hypre, Krylov pipeliné, Ginkgo
distribué (MFEM 4.10).

Tests : chaque solveur sur un cas où il converge, contre une solution de
référence ; chaque option avec un test qui échoue si elle est ignorée
(nombre de constructions du préconditionneur pour `PreconditionerUpdate`,
nombre d'itérations pour Anderson...).

Fin : vague 1 complète et mesurée ; vagues suivantes au fil de L9.

## L6. GPU

Objectif : tout le calcul d'une itération sur le GPU, données comprises.

Tranche minimale (dès la fin de L3) : petites déformations 3D, une loi J2
compilée pour `cuda`, option c puis a.

Travail :

- périphérique `cuda` de MFEM, tableaux des matériaux sur le GPU (L1),
  noyaux de L3 sur le GPU, intégration par la bibliothèque `cuda` (ou
  `kokkos`) de la loi ;
- refus d'une bibliothèque `generic` sur un périphérique GPU, sauf
  configuration mixte explicite (matériau en mémoire partagée) ;
- conditions de Dirichlet appliquées sur le GPU
  (`Vector::SetSubVector` s'exécute sur le périphérique depuis MFEM 4.9) ;
  plus d'écriture hôte dans `u` ;
- fin de pas sur le GPU par `MFEMFieldOperations` ;
- option c sur GPU : matrice globale assemblée sur le GPU à partir des
  matrices élémentaires (remplissage de matrice creuse de MFEM sur le
  périphérique, ou IJ de hypre sur GPU) ;
- plusieurs matériaux : appels asynchrones de MGIS puis une seule
  synchronisation ;
- vérification qu'aucune copie cachée n'a lieu : périphérique `debug` de
  MFEM, GPU sans HMM, profil Nsight Systems.

Tests : les tests de L3 et L4 sur le périphérique `cuda`, comparés au CPU
à l'arrondi près ; un test sous le périphérique `debug`.

Fin : mêmes résultats qu'en CPU sur les cas de référence, aucune copie
cachée, premières mesures.

## L7. Tout mfem-mgis

Objectif : le chemin par lots couvre ce que couvre le chemin actuel.

Ordre proposé :

1. déformations planes et contraintes planes ;
2. grandes déformations (gradient de la transformation, contrainte de
   Piola-Kirchhoff, conversions de MGIS sur vues) ;
3. thermique non linéaire stationnaire et transitoire ;
4. lois orthotropes : matrices de rotation par point en tableau, rotations
   de MGIS sur vues (plus de pointeurs de fonction) ;
5. propriétés matériau et variables externes sous toutes leurs formes,
   dont les propriétés `generic-parallel` évaluées par MGIS sur le GPU ;
6. conditions aux limites de bord (pression, source de chaleur) par des
   noyaux de bord, pour le GPU ;
7. micromorphe ;
8. régularisations (FBar, Faltus) : elles demandent la solution courante
   à l'assemblage, à traiter en dernier ;
9. problèmes périodiques.

Tests : la suite de mfem-mgis passée sur le chemin par lots, cas par cas.

Fin : suite complète sur le chemin par lots, sur CPU et sur GPU.

## L8. Sorties et post-traitements

- D'abord : copie vers l'hôte avant chaque post-traitement (`HostRead`),
  comportement actuel inchangé.
- Ensuite, selon les mesures : forces thermodynamiques moyennes, résultante
  sur un bord et énergies par des réductions sur le GPU ; projections L2
  par l'assemblage partiel de MFEM.

Tests : sorties identiques sur CPU et sur GPU.

## L9. Exemples et mesures

- mfem-mgis-examples et mm-opera-hpc : options `--device`,
  `AssemblyStrategy`, `TangentOperator`, `QuadratureDataLayout`, solveur,
  passées par la ligne de commande.
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
- GPU : mêmes tests sur le périphérique `cuda`, comparés au CPU ; un test
  sous le périphérique `debug` de MFEM.
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
et GPU.

Comparaisons :

- CPU (1 à 16 cœurs, MPI) contre GPU ;
- chemin actuel contre chemin par lots sur CPU ;
- rangement point par point contre composante par composante ;
- options a, b, c et éléments de l'opérateur tangent ;
- solveurs, vague par vague ;
- bibliothèque `cuda` contre `kokkos` de la même loi, et `generic` sur
  l'hôte.

## Questions à trancher

1. **Noms des paramètres** proposés plus haut.
2. **Ordre des points d'entrée** : opérateur propre d'abord (proposé), à
   cause des intégrateurs de bord ignorés par l'assemblage partiel de
   MFEM 4.9.
3. **MFEM 4.10** : second environnement pour cuDSS et Ginkgo distribué
   (voir [solveurs-gpu.md](solveurs-gpu.md#mfem-49-ou-410)), ou passage de
   tout le prototype à 4.10.
4. **Noyaux écrits à la main ou générés** (extension du générateur GiNaC
   de mfem-mgis).
5. **Propriétaire de la mémoire** : mfem-mgis (proposé) ou MGIS avec une
   ressource mémoire ; les deux sont prévus en L1, à comparer.
6. **Propositions amont** : à quel moment et dans quel ordre.

## Risques

- **Taille** : beaucoup d'intégrateurs et de cas particuliers. Le chemin
  par lots commence par les cas les plus utilisés et reste optionnel.
- **Assemblage partiel de MFEM 4.9** : limites décrites plus haut ; le
  point d'entrée B ne dépend pas de l'extension.
- **Coût de l'AMG sur GPU** : construction plus chère que la résolution ;
  d'où la réutilisation du préconditionneur en vague 1.
- **Double précision sur GPU grand public** : la RTX 5080 calcule la
  double précision 64 fois plus lentement que la simple ; les mesures
  doivent être confirmées sur un GPU de calcul (A100 ou H100).
- **Compilation de MFEM avec CUDA 13** : environnement `mmm-gpu` en cours
  de construction ; MFEM 4.9 précède CUDA 13 de peu.
- **Échecs locaux de l'intégration** : sur GPU, un point en échec doit
  remonter sans copie de tout l'état ; MGIS rend le premier point en échec
  et un statut global.
