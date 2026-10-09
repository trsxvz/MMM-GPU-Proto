# Plan détaillé pour mfem-mgis

Statut : plan établi le 09/10/2026, relu de façon indépendante le même
jour (constats intégrés). Rien n'est commencé dans mfem-mgis. Les choix
de fond sont dans [plan-mfem-mgis.md](plan-mfem-mgis.md) ; la liste des
solveurs dans [solveurs-gpu.md](solveurs-gpu.md). Ce document découpe le
travail en lots, dit quels fichiers changent, ce qui est testé et quand
un lot est fini.

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
  absolue (`mfem::HypreFGMRES` n'a pas de `SetAbsTol` en MFEM 4.9 : passer
  par `HYPRE_ParCSRFlexGMRESSetAbsoluteTol`).

### GPU

- Aucun usage de `mfem::Device`, `AssemblyLevel`, `forall`,
  `ElementRestriction`, `QuadratureInterpolator`, `GeometricFactors` ou
  `DofToQuad`. `mfem_mgis::initialize` ne crée pas de périphérique.
- `mfem::forall` ne lance un noyau CUDA que si le fichier est compilé
  par nvcc (`__CUDACC__`, `general/forall.hpp:789-794`) ; sinon le corps
  s'exécute sur l'hôte. Les fichiers de mfem-mgis qui contiennent des
  noyaux doivent donc être compilés comme du CUDA quand MFEM l'est.

### Ce que l'assemblage partiel de MFEM 4.9 impose

Vérifié dans les sources de MFEM 4.9 (`fem/nonlinearform.cpp`,
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
- **Travail dans `mfem-mgis/`, branche `gpu-proto`**, base `a011d97d6f`,
  copié dans ce dépôt par `scripts/sync-repository.sh`. Messages de commit
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
| solveur non linéaire | `FixedPoint` | `AndersonDepth`, `AndersonDamping`, `PreconditionerUpdate` |
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
| L0 | construction, référence, périphérique, petites corrections | environnement `mmm-gpu` | S | tests existants verts avec MFEM CUDA, temps de référence |
| L1 | données aux points dans la mémoire de MFEM, deux rangements | L0 | M | suite verte, rangement composante par composante sur les données |
| L2 | lois intégrées par lots sur CPU | L1 | M | mêmes résultats, gain mesuré avec `stlpar` et `openmp` |
| L3 | noyaux par lots, opérateur propre, assemblage partiel | L2 | L | résidu et tangente égaux au chemin actuel, suite verte |
| L4 | opérateur tangent (a, b, c, éléments), prédiction, règle d'intégration | L3 | M | les options sur les cas de référence |
| L5 | solveurs, en trois vagues | L3, L4 | L | vague 1 complète |
| L6 | GPU | L3 (tranche minimale), puis L4, L5 | L | mêmes résultats qu'en CPU, aucune copie cachée en un processus |
| L7 | toutes les hypothèses, cinématiques et modèles | L3 à L6 | L | suite complète de mfem-mgis sur le chemin par lots |
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
- CMake de mfem-mgis : langage CUDA activé quand MFEM est compilé avec
  CUDA, pour les fichiers qui contiendront des noyaux.
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
  au chemin par lots. Dans ce lot, seules les données des matériaux et
  leur description changent ; les accès du chemin actuel restent point
  par point.
- Fin de pas : `update` et `revert` passent par les surcharges de MGIS qui
  prennent des `FieldOperations`. MGIS fournit `HostFieldOperations` et
  des opérations CUDA, HIP, SYCL, Kokkos et OpenMP en en-têtes, à compiler
  par mfem-mgis avec la chaîne correspondante. Variante à comparer :
  `MFEMFieldOperations`, écrites avec `mfem::forall` et les copies de
  `mfem::Memory` (compilées comme du CUDA quand MFEM l'est).

Tests : suite complète inchangée (rangement point par point, données
dans `mfem::Memory`) ; refus du rangement composante par composante avec
`Legacy` ; tests unitaires des opérations de fin de pas dans les deux
rangements ; test de la politique de validité sous le périphérique
`debug` de MFEM.

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
  doit être conservé.

Tests : suite complète avec une bibliothèque `generic` (identique au
chemin actuel), puis avec `stlpar` et `openmp` (identique à l'arrondi
près) ; loi orthotrope (rotations) ; plusieurs intégrateurs sur un
attribut ; test d'un point en échec.

Fin : mêmes résultats, gain mesuré sur la poutre et le polycristal.

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
  par géométrie ;
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
  - degrés de liberté imposés : résidu nul, tangente éliminée ou
    `ConstrainedOperator`, diagonale à 1 pour les lisseurs ;
  - intégrateurs de bord : gardés sur le chemin actuel pour le résidu et
    la tangente assemblée ; une action sans matrice est à écrire pour
    l'option a ;
  - l'intégrateur de domaine actuel n'est plus ajouté au formulaire
    (sinon le domaine serait compté deux fois).

  `Mult` intègre la loi puis assemble le résidu (option), ce qui rend
  possibles les solveurs non linéaires extérieurs ; un échec de la loi
  est signalé sans arrêter le programme. `mgis_integrator` (type `final`)
  devient une interface commune aux deux chemins, ou le chemin par lots a
  son propre objet.
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
  `ImposedPressureTest` avec les deux points d'entrée (le second doit
  refuser la pression) : mêmes itérations, mêmes résultats.

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
  valeurs changent).
- `TangentOperator=MatrixFree` (a) : action des noyaux, diagonale et
  blocs nodaux pour les lisseurs, AMG construit une fois sur l'opérateur
  élastique assemblé.
- `TangentOperator=ElementMatrices` : matrices élémentaires gardées,
  action élément par élément.
- `TangentOperator=ElasticOperator` (b) : avec `FixedPoint` (L5).
- `IntegrationRuleOrder` : règle choisie par l'utilisateur, par exemple
  8 points sur l'hexaèdre d'ordre 1 au lieu de 27 ; vérification de
  l'absence de modes parasites sur les cas de référence.
- Prédiction (`...Implementation.cxx:98-233`) réécrite avec les noyaux :
  opérateur élastique, sécant ou tangent par l'intégration de prédiction
  de MGIS, second membre par le noyau de résidu, et les quatre axes
  actuels conservés. Nouvelle stratégie : extrapolation de l'incrément du
  pas précédent. Paramètre `Prediction`.

Tests : chaque option sur les cas de référence (mêmes solutions à la
tolérance de Newton près) ; un test par option qui échoue si elle est
ignorée (nombre d'itérations, type de l'opérateur rendu, nombre de points
d'intégration) ; prédictions contre `PredictionTest` et `PredictionTest2`.

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
- Newton inexact (Eisenstat-Walker) et recherche linéaire, écrits dans le
  Newton de mfem-mgis (les crochets de `mfem::NewtonSolver` ne s'y
  appliquent pas) ;
- lisseurs sans matrice : `OperatorJacobi`, `OperatorChebyshev`,
  `NodalBlockJacobi` ;
- options de BoomerAMG exposées, sans être écrasées par
  `SystemHypreBoomerAMG::SetOperator` ; stratégie `Elasticity` refusée sur
  GPU (MFEM arrête le programme en byVDIM), `System` proposée à la place.

**Vague 2** (comparaisons sur GPU) :

- L-BFGS préconditionné, à écrire (`LBFGSSolver` de MFEM ignore le
  préconditionneur ; `SNESQN` de PETSc comme référence) ;
- Ginkgo dans la fabrique : préconditionneurs (Jacobi par blocs à
  précision adaptative, ILU, ISAI) appliqués au bloc local de chaque
  processus en MFEM 4.9 ; solveurs de Ginkgo en séquentiel seulement ;
- AmgX par `AmgXSolver`, avec `UpdateOperator` pour la réutilisation
  complète ;
- KINSOL (`KINSolver`), avec l'opérateur qui intègre la loi (L3) ;
- réparation du chemin PETSc (même opérateur, plus d'arrêt sur échec de
  la loi), puis GAMG avec espace quasi nul ;
- précision mixte (Ginkgo, puis AmgX `dDFI` par une enveloppe propre).

**Vague 3** (selon les mesures) : p-multigrille et LOR (ordre 2 et plus),
cuDSS, SuperLU_DIST et STRUMPACK sur GPU, MGR de hypre pour les
problèmes couplés, FSAI et ILU de hypre, Krylov pipeliné, Ginkgo
distribué (MFEM 4.10), réutilisation partielle de la hiérarchie AMG
(GAMG, AmgX par une enveloppe propre).

Tests : chaque solveur sur un cas où il converge, contre une solution de
référence ; chaque option avec un test qui échoue si elle est ignorée
(nombre de constructions du préconditionneur pour `PreconditionerUpdate`,
nombre d'itérations pour Anderson, réglages de BoomerAMG relus après
plusieurs itérations...).

Fin : vague 1 complète et mesurée ; vagues suivantes au fil de L9.

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
  depuis MFEM 4.9) ; plus d'écriture hôte dans `u` ;
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
3. thermique non linéaire stationnaire et transitoire ;
4. lois orthotropes : matrices de rotation par point en tableau, rotations
   de MGIS sur vues (plus de pointeurs de fonction) ;
5. propriétés matériau et variables externes sous toutes leurs formes,
   dont les propriétés `generic-parallel` évaluées par MGIS sur le GPU ;
6. conditions aux limites de bord (pression, source de chaleur) par des
   noyaux de bord, pour le GPU ;
7. modèles et schémas de couplage (`Simulation`, `IterativeCouplingScheme`,
   `LoopCouplingScheme`, `NonLinearModel`, `PointWiseModel`,
   `QPEvaluators`) : l'exemple thermomécanique ex8 en dépend, et les
   échanges de champs aux points entre modèles se font sur l'hôte.
   `PointWiseModel` (loi aux points sans éléments finis) est un premier
   cas simple de loi par lots sur GPU ;
8. micromorphe ;
9. régularisations (FBar, Faltus) : elles demandent la solution courante
   à l'assemblage, à traiter en dernier ;
10. problèmes périodiques.

Tests : la suite de mfem-mgis passée sur le chemin par lots, cas par cas.

Fin : suite complète sur le chemin par lots, sur CPU et sur GPU.

## L8. Sorties et post-traitements

- D'abord : copie vers l'hôte avant chaque post-traitement (`HostRead`),
  comportement actuel inchangé. Avec le rangement composante par
  composante, les accès des post-traitements sont à adapter ou les
  données à transposer vers une copie point par point.
- Ensuite, selon les mesures : forces thermodynamiques moyennes, résultante
  sur un bord et énergies par des réductions sur le GPU ; projections L2
  par l'assemblage partiel de MFEM.

Tests : sorties identiques sur CPU et sur GPU, dans les deux rangements.

## L9. Exemples et mesures

- mfem-mgis-examples et mm-opera-hpc : options `--device`,
  `AssemblyStrategy`, `TangentOperator`, `QuadratureDataLayout`,
  `IntegrationRuleOrder`, solveur, passées par la ligne de commande.
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
   partiel de MFEM 4.9.
3. **MFEM 4.10** : second environnement pour cuDSS et Ginkgo distribué
   (voir [solveurs-gpu.md](solveurs-gpu.md#mfem-49-ou-410)), ou passage de
   tout le prototype à 4.10.
4. **Noyaux écrits à la main ou générés** (extension du générateur GiNaC
   de mfem-mgis).
5. **Propriétaire de la mémoire** : mfem-mgis (proposé) ou MGIS avec une
   ressource mémoire ; les deux sont prévus en L1, à comparer.
6. **Base de travail** : `a011d97d6f`, ou la refonte de
   `PartialQuadratureSpace` en cours en amont une fois fusionnée.
7. **Propositions amont** : à quel moment et dans quel ordre.

## Risques

- **Taille** : beaucoup d'intégrateurs et de cas particuliers. Le chemin
  par lots commence par les cas les plus utilisés et reste optionnel.
- **Branches amont en cours** : la refonte de `PartialQuadratureSpace`
  (branche `466-refactor-partial-quadrature-space-to-better-support-spaces-on-boundaries`,
  89 fichiers) touche les classes que L1 modifie ; la branche
  `452-post-processings-only-consider-the-first-behaviour-integrator-of-a-material`
  touche les post-traitements de L1 et L8. Conflits à prévoir.
- **Assemblage partiel de MFEM 4.9** : limites décrites plus haut ; le
  point d'entrée B ne dépend pas de l'extension.
- **Compilation CUDA dans mfem-mgis** : les noyaux demandent nvcc ;
  temps de compilation et mémoire à surveiller, comme pour MFEM.
- **Coût de l'AMG sur GPU** : construction plus chère que la résolution ;
  d'où la réutilisation du préconditionneur en vague 1.
- **Double précision sur GPU grand public** : la RTX 5080 calcule la
  double précision 64 fois plus lentement que la simple ; les mesures
  doivent être confirmées sur un GPU de calcul (A100 ou H100).
- **Compilation de MFEM avec CUDA 13** : environnement `mmm-gpu` en cours
  de construction ; MFEM 4.9 (décembre 2025) est postérieur à CUDA 13.0
  (août 2025), mais CUDA 13 a changé des en-têtes (CCCL 3).
- **Échecs locaux de l'intégration** : sur GPU, un point en échec doit
  remonter sans copie de tout l'état ; MGIS rend le premier point en échec
  et un statut global.
- **Plusieurs processus sans MPI compatible GPU** : copies vers l'hôte à
  chaque échange, inévitables avec l'environnement local.
