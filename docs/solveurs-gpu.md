# Solveurs et techniques pour l'étude CPU et GPU de mfem-mgis

Statut : liste établie le 09/10/2026. Rien n'est codé. Cette liste sert à
choisir les études de performance de mfem-mgis et les paquets de
l'environnement spack `mmm-gpu` ([deps/spack](../deps/spack/README.md)).

Périmètre : seulement ce qui tourne sur CPU et sur GPU, ou sur GPU seul.
Ce qui ne tourne que sur CPU est écarté (voir [Écartés](#écartés-cpu-seul)).

Chaque point indique où la fonctionnalité existe, ce qui a été vérifié et
comment, ce qu'il faudrait dans mfem-mgis, et un intérêt sur 10 pour ce
projet. « Vérifié » veut dire lu dans les sources ou les notes de version
de la version citée (voir [Sources](#sources)), pas encore essayé sur la
machine. « À vérifier » veut dire non confirmé.

## Sommaire

1. [Tableau d'ensemble](#tableau-densemble)
2. [Réutiliser les préconditionneurs](#réutiliser-les-préconditionneurs)
3. [Solveurs non linéaires](#solveurs-non-linéaires)
4. [Préconditionneurs sans matrice](#préconditionneurs-sans-matrice)
5. [AMG sur GPU](#amg-sur-gpu)
6. [Précision mixte](#précision-mixte)
7. [Solveurs directs sur GPU](#solveurs-directs-sur-gpu)
8. [Plusieurs GPU](#plusieurs-gpu)
9. [Assemblage et opérateur tangent](#assemblage-et-opérateur-tangent)
10. [Problèmes couplés](#problèmes-couplés)
11. [Outils de mesure](#outils-de-mesure)
12. [Écartés (CPU seul)](#écartés-cpu-seul)
13. [MFEM 4.10](#mfem-410)
14. [Correspondance avec l'environnement spack](#correspondance-avec-lenvironnement-spack)
15. [Sources](#sources)

## Tableau d'ensemble

Les lignes sans marque viennent de la liste de départ. Les lignes marquées
« ajout » complètent cette liste. Ginkgo est déjà retenu : il apparaît dans
plusieurs lignes.

| Technique | CPU/GPU | Intérêt | Où |
|---|---|---|---|
| Réutiliser le préconditionneur (AMG construit une fois, sur l'opérateur élastique, ou gardé plusieurs itérations ou pas de temps) | les deux | 9 | mfem-mgis : enveloppe qui ignore `SetOperator` (voir plus bas) |
| Ajout : réutilisation partielle de la hiérarchie AMG (prolongateurs gardés, opérateurs grossiers recalculés) | les deux | 7 | GAMG de PETSc ; AmgX (`AMGX_solver_resetup`, pas exposé par MFEM) |
| Point fixe accéléré par Anderson | les deux | 9 | mfem-mgis, KINSOL |
| Ajout : assemblage de la matrice tangente globale sur GPU (option c sur GPU), motif creux et matrice hypre réutilisés d'une itération à l'autre | les deux | 9 | mfem-mgis, hypre (IJ sur GPU), MFEM |
| Ajout : choix de la règle d'intégration (8 points au lieu de 27 sur l'hexaèdre d'ordre 1) | les deux | 8 | mfem-mgis |
| Newton inexact (Eisenstat-Walker, `NewtonSolver::SetAdaptiveLinRtol` de MFEM) | les deux | 8 | fait dans le Newton de mfem-mgis et la méthode `NewtonK0` (option `AdaptiveLinearSolverTolerance`) |
| Lisseurs sans matrice Jacobi et Chebyshev sur la diagonale (`OperatorJacobiSmoother`, `OperatorChebyshevSmoother`), indispensables à l'option a entièrement sans matrice | les deux | 8 | MFEM |
| Ajout : Jacobi par blocs nodaux 3×3 sans matrice (diagonale par blocs de la tangente) | les deux | 7 | mfem-mgis, comme la mini-application `hooke` de MFEM |
| Ajout : l1-Jacobi « en valeur absolue » sans matrice | les deux | 5 | MFEM 4.9, `miniapps/diag-smoothers` |
| Réglages GPU de BoomerAMG exposés dans mfem-mgis (PMIS, ext+i, l1-Jacobi, Chebyshev, nombre de niveaux, seuils) | les deux | 8 | hypre, MFEM |
| Ajout : BoomerAMG pour les systèmes sur GPU (stratégie « System », réglages préservés à chaque `SetOperator` ; « Elasticity » arrête le programme avec hypre GPU) | les deux | 8 | MFEM, mfem-mgis |
| Précision mixte : préconditionneur en simple précision, Krylov en double ; gain fort sur la RTX 5080 (FP64 lent), plus faible sur A100 ou H100 | surtout GPU | 6 (8 sur GPU grand public) | Ginkgo, hypre 3, AmgX |
| L-BFGS préconditionné par un cycle AMG ou p-MG (`LBFGSSolver` de MFEM), comme Ratel | les deux | 7 | à écrire : `LBFGSSolver` de MFEM ignore le préconditionneur ; `SNESQN` de PETSc en a un |
| p-multigrille sans matrice avec AMG sur le Q1 (ordre 2 et plus), ou LOR (`LORSolver`) | les deux | 7 | MFEM |
| Ajout : Ginkgo distribué (matrices MPI, Schwarz à deux niveaux, multigrille PGM distribuée) | les deux | 7 | Ginkgo 1.8 à 1.11 ; enveloppes MFEM 4.10 |
| Ajout : AmgX en précision mixte (mode `dDFI`) et matrices par blocs 3×3 | GPU NVIDIA | 6 | AmgX ; pas dans l'enveloppe de MFEM |
| Recherche linéaire pour Newton (crochet `ComputeScalingFactor`) | les deux | 6 | à écrire dans le Newton de mfem-mgis (il ne divise aujourd'hui la correction par 2 qu'en cas d'échec de la loi) |
| AMG élasticité avec modes rigides sur GPU : GAMG de PETSc avec espace quasi nul (Kokkos), nouveau chemin par blocs sur GPU (Adams, arXiv 2606.24748) ; demande de réparer le chemin PETSc de mfem-mgis | les deux | 6 | PETSc 3.26 (`aijkokkos` ; le chemin par blocs de l'article n'y est pas) |
| KINSOL de SUNDIALS par MFEM (`KINSolver` : point fixe avec Anderson, JFNK, recherche linéaire) | les deux | 6 | SUNDIALS 7.9, MFEM ; demande un opérateur dont `Mult` intègre la loi |
| Ajout : prédiction par extrapolation de l'incrément de déplacement | les deux | 5 | mfem-mgis |
| AmgX par `AmgXSolver` de MFEM (point de comparaison avec BoomerAMG, modes MPI) | GPU NVIDIA | 5 | AmgX, MFEM |
| cuDSS (solveur direct NVIDIA, un ou plusieurs GPU, MPI, refactorisation) | GPU NVIDIA | 5 | cuDSS ; enveloppe `CuDSSSolver` de MFEM 4.10 ; fait : solveur `CuDSSSolver` de mfem-mgis |
| Ajout : Krylov pipeliné (CG pipeliné de Ginkgo, `KSPPIPECG` de PETSc) | les deux | 5 | Ginkgo 1.10, PETSc |
| Ajout : Ratel comme calcul de référence sur le même cas | les deux | 5 | Ratel, PETSc, libCEED |
| SuperLU_DIST et STRUMPACK sur GPU (enveloppes MFEM existantes) | les deux | 4 | MFEM |
| MGR de hypre (GPU depuis 2.29) pour les problèmes couplés monolithiques | les deux | 4 | hypre (construction sur GPU depuis 2.24, sans mémoire unifiée depuis 2.29) |
| FSAI (GPU depuis 2.30) et ILU de hypre sur GPU | les deux | 4 | hypre |
| Ajout : hypre récent (filtrage des couplages entre composantes, cycle flexible, Schwarz avec recouvrement, solveur de Chebyshev) | les deux | 4 | hypre 2.32 à 3.2 |
| Ajout : pas de temps adaptatif guidé par la convergence de Newton | les deux | 4 | mfem-mgis |
| Ajout : solveurs non linéaires de PETSc sur vecteurs GPU (NGMRES, Anderson, quasi-Newton, NASM) | les deux | 4 | PETSc SNES |
| Ajout : MueLu de Trilinos (agrégation lissée avec modes rigides, Kokkos) | les deux | 3 | Trilinos ; pas d'interface MFEM |
| Ajout : recyclage de sous-espaces de Krylov entre itérations de Newton | les deux | 3 | PETSc (HPDDM, GCRO-DR) ; GPU à vérifier |
| Ajout : lois de comportement en simple précision sur GPU grand public | GPU | 3 | demande TFEL et MGIS |
| Ajout : algèbre dense par lots de MFEM (`BatchedLinAlg`, blocs 3×3 et matrices élémentaires) et graphes CUDA (petits noyaux) | surtout GPU | 4 | MFEM 4.8 et plus, CUDA |
| Ajout (outils) : profils NVTX de hypre, Caliper, Nsight Systems | les deux | 7 | hypre `+gpu-profiling`, MFEM |
| Ajout (outils) : pools mémoire Umpire pour MFEM et hypre | surtout GPU | 6 | Umpire |
| Ajout (outils) : MPI compatible GPU | GPU | 5 | MPI de la grappe, MFEM, hypre |

## Réutiliser les préconditionneurs

Sur GPU, la construction d'un AMG coûte souvent plus que plusieurs
résolutions. mfem-mgis la refait à chaque itération de Newton : le
solveur linéaire reçoit la nouvelle tangente
(`NewtonSolver.cxx:245`), et le `SetOperator` d'un Krylov rappelle
toujours celui de son préconditionneur (Krylov de MFEM,
`linalg/solvers.cpp:184-193` ; `HyprePCG::SetOperator` ; la phase de
construction des Krylov de hypre reconstruit aussi le préconditionneur).

- **Réutilisation complète** : construire l'AMG une fois et le garder.
  Variantes : sur l'opérateur élastique (garde toute la simulation), sur
  la tangente du début du pas (garde un pas), sur la tangente gardée
  `n` itérations. Le Krylov reste exact, seul le préconditionneur
  vieillit. Il faut une enveloppe de préconditionneur dont `SetOperator`
  ne reconstruit rien (et, avec les Krylov de hypre, un préconditionneur
  dont la phase de construction ne fait rien), plus la politique de mise
  à jour dans le Newton de mfem-mgis. Des mécanismes existent ailleurs :
  `KSPSetReusePreconditioner` et `-snes_lag_preconditioner` de PETSc,
  `SetMaxSetupCalls` de KINSOL. Intérêt 9.
- **Réutilisation partielle** (ajout) : garder la structure de la
  hiérarchie (agrégats, prolongateurs) et recalculer seulement les
  opérateurs grossiers `Pᵀ A P` avec les nouvelles valeurs. Coût de
  construction bien plus faible, qualité proche d'une construction
  complète tant que la tangente ne change pas trop. Existe dans :
  - GAMG de PETSc (`PCGAMGSetReuseInterpolation`) ;
  - AmgX (`AMGX_solver_resetup`, option `structure_reuse_levels`), pas
    exposé par MFEM : `AmgXSolver::UpdateOperator` remplace seulement les
    coefficients de la matrice, sans nouvelle construction, ce qui est
    une réutilisation complète (vérifié dans `linalg/amgxsolver.cpp` de
    MFEM 4.9) ;
  - MGR de hypre 3.2 (« component reuse », notes de version) ;
  - BoomerAMG : pas d'équivalent trouvé.
  Intérêt 7.

## Solveurs non linéaires

### Prérequis : un opérateur qui intègre la loi

Dans mfem-mgis, `Mult` (le résidu) n'intègre pas la loi : l'intégration
est déclenchée par le solveur non linéaire de mfem-mgis, par une action
appelée à chaque nouvelle estimation (`NonLinearSolverBase.cxx:39-58`).
Le `Mult` du formulaire ne fait qu'assembler les forces déjà calculées.
Un solveur non linéaire extérieur (`LBFGSSolver`, `KINSolver`,
`mfem::NewtonSolver`, SNES de PETSc) appelle seulement `Mult` : il
calculerait le résidu avec des contraintes périmées. Seule exception, le
mode PETSc, où `AssembleElementVector` intègre la loi lui-même ; un échec
de la loi y arrête le programme (`MFEM_VERIFY`,
`MultiMaterialNonLinearIntegrator.cxx:290-297`), ce qui explique
peut-être que le chemin PETSc plantait.

Il faut donc, pour tous les solveurs extérieurs, un opérateur dont `Mult`
intègre la loi puis assemble, et qui signale un échec de la loi sans
arrêter le programme. C'est le rôle de l'opérateur propre du plan
détaillé (point d'entrée B).

Le Newton de mfem-mgis (`mfem_mgis::NewtonSolver`) n'est pas celui de
MFEM : il dérive de `mfem::IterativeSolver` par `NonLinearSolverBase`,
et réécrit `Mult`. Les crochets de `mfem::NewtonSolver`
(`SetAdaptiveLinRtol`, `ComputeScalingFactor`) n'y servent pas.

### Techniques

Tous ces solveurs n'utilisent que des vecteurs, des actions d'opérateurs
et des produits scalaires : ils tournent sur GPU dès que les vecteurs et
les opérateurs y sont.

- **Point fixe accéléré par Anderson** (intérêt 9). À écrire dans
  mfem-mgis (le solveur de point fixe du plan, option b), ou par KINSOL
  (`KINSolver::EnableAndersonAcc`, pour `KIN_FP` et `KIN_PICARD`, vérifié
  dans `linalg/sundials.hpp` de MFEM 4.9). Avec `KIN_FP`, KINSOL attend la
  fonction de point fixe `G(u) = u - K0⁻¹ F(u)`, pas le résidu : il faut
  un opérateur qui l'applique, ou `KIN_PICARD`. Le point fixe seul
  convergeait mal sur le polycristal des essais préalables : Anderson est
  la première chose à essayer.
- **Newton inexact** (intérêt 8) : tolérance du solveur linéaire adaptée
  au résidu non linéaire (Eisenstat-Walker). À écrire dans le Newton de
  mfem-mgis, sur le modèle de `mfem::NewtonSolver::SetAdaptiveLinRtol`.
- **L-BFGS** (intérêt 7) : avec un cycle d'AMG ou de p-multigrille comme
  inverse initial, comme Ratel. `LBFGSSolver` de MFEM ne convient pas :
  il ignore le préconditionneur (`SetPreconditioner` affiche « L-BFGS
  won't use the given preconditioner », `linalg/solvers.hpp:917-920`, en
  4.9 comme en 4.10). À écrire ; `SNESQN` de PETSc
  (`-snes_qn_type lbfgs`, `-snes_qn_scale_type jacobian`) est la
  référence existante d'un L-BFGS à inverse initial préconditionné.
- **Recherche linéaire** (intérêt 6) : le Newton de mfem-mgis divise
  aujourd'hui la correction par 2 seulement quand la loi échoue
  (`NewtonSolver.cxx:166-195`). Une vraie recherche linéaire (sur la
  norme du résidu ou l'énergie) est à y écrire.
- **KINSOL** (intérêt 6) : point fixe avec Anderson, Picard, Newton-Krylov
  sans matrice (`SetJFNK`), recherche linéaire, jacobienne réutilisée
  (`SetMaxSetupCalls`). L'enveloppe de MFEM utilise les vecteurs CUDA de
  SUNDIALS quand MFEM est compilé avec CUDA (vérifié dans
  `linalg/sundials.hpp`). Demande l'opérateur qui intègre la loi. À
  essayer sur la machine.
- **Prédiction par extrapolation** (ajout, intérêt 5) : début de chaque pas
  à partir de l'incrément du pas précédent (linéaire), ou des deux
  précédents (quadratique). Moins chère que la prédiction par l'opérateur
  tangent de mfem-mgis : aucune résolution. Réduit le nombre
  d'itérations de Newton quand le chargement est régulier.
- **Pas de temps adaptatif** (ajout, intérêt 4) : mfem-mgis a le cadre
  (`AbstractTimeIncrementComputer`, sous-pas en cas d'échec), mais le
  calcul par défaut va directement à la fin de la séquence. Un calcul du
  pas guidé par le nombre d'itérations de Newton influence directement le
  coût total, donc la comparaison des variantes.
- **Solveurs non linéaires de PETSc** (ajout, intérêt 4) : NGMRES,
  Anderson, quasi-Newton (L-BFGS, `SNESQN`), NASM, avec des vecteurs sur
  GPU. Utiles comme implémentations de référence pour les points
  ci-dessus, une fois le chemin PETSc de mfem-mgis réparé (opérateur qui
  intègre la loi).

## Préconditionneurs sans matrice

Utiles à l'option a (tangente sans matrice). L'option a du plan
préconditionne par un AMG construit sur l'opérateur élastique assemblé ;
les lisseurs ci-dessous, qui ne demandent que la diagonale ou des blocs,
permettent une option a entièrement sans matrice (seuls, ou dans une
p-multigrille).

- **Jacobi et Chebyshev sur la diagonale** (intérêt 8) :
  `OperatorJacobiSmoother` et `OperatorChebyshevSmoother` de MFEM
  (vérifiés dans 4.9). La diagonale vient de `AssembleGradDiagonalPA`
  (point d'entrée de l'assemblage partiel) ou d'un noyau de mfem-mgis.
  Chebyshev demande une estimation de la plus grande valeur propre.
- **Jacobi par blocs nodaux 3×3** (ajout, intérêt 7) : inverser, nœud par
  nœud, le bloc 3×3 de la diagonale de la tangente. Plus robuste que le
  Jacobi scalaire en élasticité, où les composantes d'un même nœud sont
  fortement couplées. La mini-application `hooke` de MFEM 4.9 (élasticité
  non linéaire sans matrice, sur GPU) le fait
  (`ElasticityDiagonalPreconditioner`, type `BlockDiagonal`) : c'est un
  modèle direct pour mfem-mgis. L'algèbre dense par lots de MFEM
  (`BatchedLinAlg`, depuis 4.8) peut inverser ces blocs sur GPU.
- **l1-Jacobi en valeur absolue** (ajout, intérêt 5) : famille de lisseurs
  sans matrice convergents par construction (MFEM 4.9,
  `miniapps/diag-smoothers`, exemples en élasticité et en multigrille).
- **p-multigrille et LOR** (intérêt 7) : `GeometricMultigrid` de MFEM avec
  des lisseurs sans matrice et un AMG sur l'espace d'ordre 1, ou
  `LORSolver` (AMG sur la discrétisation d'ordre bas). Seulement utile
  pour l'ordre 2 et plus. MFEM 4.10 améliore l'assemblage LOR par lots.

## AMG sur GPU

### BoomerAMG (hypre)

- Ce qui tourne sur GPU (documentation de hypre) : grossissement PMIS,
  interpolations directe (3), BAMG directe (15), étendue (14), étendue+i
  (6) et étendue+e (18), grossissement agressif avec interpolation
  multipasse ou en deux étapes, lisseurs Jacobi (7), l1-Jacobi (18),
  Gauss-Seidel hybride (3, 4, 6), Gauss-Seidel en deux étapes (11, 12) et
  Chebyshev (16). L'environnement compile hypre sans mémoire unifiée : une
  option sans support GPU échoue au lieu de passer sur le CPU.
- MFEM choisit des valeurs compatibles avec le GPU quand hypre est compilé
  pour le GPU (PMIS, ext+i, l1-Jacobi, sans grossissement agressif ;
  vérifié dans `HypreBoomerAMG::SetDefaultOptions`, MFEM 4.9).
- **Réglages exposés dans mfem-mgis** (intérêt 8) : aujourd'hui seule la
  « stratégie » est réglable. Il faut exposer grossissement, interpolation,
  lisseur, nombre de balayages, seuil, nombre de niveaux, grossissement
  agressif, et les comparer sur GPU.
- **Stratégies de mfem-mgis sur GPU** (ajout, intérêt 8) :
  - « System » (`SetSystemsOptions`, approche par inconnue) fonctionne sur
    GPU : MFEM copie la table des composantes sur le GPU. C'est le défaut
    de mfem-mgis sans `Strategy` quand l'espace a plusieurs composantes.
  - « Elasticity » appelle `HypreBoomerAMG::SetElasticityOptions`, qui
    arrête le programme avec hypre GPU (« this method is not supported in
    hypre built with GPU support », vérifié dans MFEM 4.9). Cela ne
    concerne que la numérotation byVDIM : en byNODES, mfem-mgis passe déjà
    à « System » avec un avertissement (`LinearSolverFactory.cxx:84-112`).
    Sur GPU, mfem-mgis doit refuser cette stratégie avec un message clair
    ou passer à « System ». Les modes rigides (vecteurs d'interpolation)
    restent sur CPU avec BoomerAMG.
  - `SystemHypreBoomerAMG` rappelle `SetSystemsOptions` à chaque
    `SetOperator` (`LinearSolverFactory.cxx:46-50`), ce qui remet le
    nombre de niveaux de grossissement agressif à 0 et le seuil à 0,5
    (`linalg/hypre.cpp:5231-5232`) : des réglages exposés seraient écrasés
    à chaque itération de Newton. À corriger avec l'exposition des
    réglages.
  - `SetDefaultOptions` choisit les valeurs pour CPU ou GPU à la
    construction de l'AMG (`linalg/hypre.cpp:5029`) : l'AMG doit être
    créé après `mfem::Device`.
- **Filtrage des couplages entre composantes** (ajout, intérêt 4) :
  `HYPRE_BoomerAMGSetFilterFunctions` (hypre 2.32) retire les couplages
  entre composantes de la matrice utilisée pour l'AMG : complexité plus
  faible, convergence à mesurer en élasticité. Suppose des inconnues
  entrelacées (byVDIM). GPU à vérifier.
- **Petits systèmes : hypre sur l'hôte** (fait le 12/10). Sous 10 000
  inconnues environ, une résolution par GMRES et BoomerAMG sur le GPU est
  dominée par les lancements de noyaux : 100 ms pour les 9 516 inconnues
  de la thermique d'ex8, contre 42 ms sur deux cœurs. L'option
  `--hypre-on-host` de `mfem_mgis::initialize` garde les lois et les
  noyaux sur le GPU et fait calculer hypre sur l'hôte : ex7 complet
  (2 634 inconnues) 6,7 -> 4,4 s, ex8 12,5 -> 11,6 s
  (`results/hypre-hote-2026-10-12`). MUMPS sur l'hôte donne le même gain
  pour la thermique d'ex8 (`results/ex8-thermique-gpu-2026-10-11`).

### GAMG (PETSc)

Agrégation lissée avec espace quasi nul (modes rigides) : l'AMG de
référence pour l'élasticité. Avec PETSc compilé avec Kokkos ou CUDA, les
matrices sont sur GPU (`aijkokkos`, `aijcusparse`) ; la part de la
construction faite sur GPU dépend de la version, à mesurer. PETSc 3.26
utilise Kokkos 5.2.2 par défaut, la version de l'environnement. Le chemin
PETSc de mfem-mgis plantait lors d'essais antérieurs : à réparer d'abord
(voir [le prérequis des solveurs non linéaires](#prérequis--un-opérateur-qui-intègre-la-loi)).
MFEM 4.10 exige PETSc compilé avec hypre. Le chemin par blocs sur GPU de
l'article d'Adams (arXiv 2606.24748) n'est pas dans PETSc 3.26.0, qui n'a
que les matrices `aijkokkos` et `aijcusparse` sur GPU. Intérêt 6.

### AmgX (NVIDIA)

- `AmgXSolver` de MFEM : matrice scalaire, double précision seulement
  (`AMGX_mode_dDDI`, commentaire « currently only double is supported »
  dans `linalg/amgxsolver.hpp`), modes MPI, réutilisation complète par
  `UpdateOperator` (coefficients remplacés, pas de nouvelle
  construction). Ses ressources, partagées par tous ses solveurs, gardent
  le communicateur et la configuration du premier : détruire le premier
  solveur fait échouer les suivants ; le GPU 0 est imposé ; le statut de
  convergence d'AmgX n'est pas accessible. **mfem-mgis ne l'utilise pas.**
- **Enveloppe propre à mfem-mgis** (fait, option `enable-amgx`, AmgX donné
  par `AMGX_DIR`) : interface C d'AmgX, préconditionneur et solveur `AmgX`
  (cycles du multigrille), solveur `AmgXPCG` (gradient conjugué d'AmgX).
  - Ressources d'AmgX uniques pour le processus : elles portent des objets
    globaux (réserve de mémoire du GPU, poignées cuSPARSE et cuBLAS),
    libérés par la destruction de n'importe quelles ressources, et gardent
    l'adresse de leur communicateur. Elles sont créées pour le premier
    solveur sur `MPI_COMM_WORLD` (AmgX en fait lui-même une copie pour
    ses échanges) et détruites avec le dernier ; création et destruction
    ne sont pas collectives.
  - Communicateur d'un solveur parallèle : tous les processus dans le
    même ordre, ou un seul processus (la matrice est alors téléversée
    comme une matrice locale). Chaque processus doit posséder des lignes
    de la matrice, au plus 2^31 - 1 lignes en tout : contrôles collectifs
    dans `SetOperator`. Après un échec de `SetOperator`, le solveur refuse
    de résoudre jusqu'à la mise en place suivante.
  - Les solveurs doivent être détruits avant `MPI_Finalize` (AmgX libère
    alors ses communicateurs).
  - Planchers d'AmgX : la convergence est déclarée dès que la norme du
    résidu passe sous le maximum de 1e-20 et de `e` fois la norme du
    résidu initial (`e` vaut 1e-12 en double, 1e-6 en simple), et
    d'emblée si la norme initiale est sous 1e-20. Sans tolérance
    absolue, `AmgXPCG` multiplie le second membre par une puissance de 2
    qui ramène le maximum des valeurs absolues de ses composantes entre
    1 et 2 (calculé sur tous les processus ; la norme elle-même
    déborderait aux extrêmes), puis la solution par l'inverse
    (opérations exactes) : seul reste le plancher relatif, et la
    tolérance relative doit valoir au moins `e`. Les statistiques
    affichées par AmgX sont alors celles du système mis à l'échelle.
    Avec une tolérance absolue (double précision seulement), le second
    membre n'est pas modifié : la tolérance absolue doit valoir au moins
    1e-20, et avec une tolérance absolue seule la résolution s'arrête
    aussi quand la norme du résidu a été multipliée par `e`. En simple
    précision, le second membre est toujours mis à l'échelle : une
    tolérance relative est exigée et une tolérance absolue refusée. Le
    préconditionneur en simple précision (configuration par défaut) met
    aussi le second membre à l'échelle, pour rester dans la plage de la
    simple précision : une tolérance passée par `AmgXParameters`
    s'applique alors au second membre mis à l'échelle.
  - Double copie d'AmgX, en attendant la reconstruction de MFEM :
    libmfem.so contient l'AmgX statique non patché de MFEM. Un programme
    doit lier libamgxsh.so (patchée) avant MFEM, comme le fait
    `MFEMMGISConfig.cmake` ; sinon le multigrille classique échoue
    (« Thrust failure »). Un programme qui se termine avec un solveur
    encore vivant (`std::exit`, objet statique) peut s'arrêter ou se
    bloquer sur une double libération des objets globaux des deux
    copies.
  - GPU courant du processus, statut de convergence d'AmgX, erreurs
    rendues par les codes de retour (sans le paramètre
    `exception_handling`, qui fait au contraire abandonner AmgX).
  - Paramètres d'AmgX par leur nom (`AmgXParameters`, par exemple
    `strength_threshold=0.8, interpolator=MULTIPASS`), qui complètent la
    configuration par défaut (celle de MFEM).
  - Simple précision (`Precision` : `Single`, mode `dFFI`, vecteurs
    convertis sur le GPU). La précision mixte `dDFI` (matrice en simple,
    vecteurs en double) n'est plus prise en charge par AmgX 2.5
    (« Mixed precision modes not currently supported for CUDA 10.1 or
    later »).
  - Reste à faire : matrices par blocs 3×3 (ordre par composantes),
    réutilisation partielle (`AMGX_solver_resetup`), mise en place sans
    passer par l'hôte.
- AmgX 2.5.0 annonce CUDA 13 (testé jusqu'à 13.0). Deux corrections
  complémentaires ont suivi sur la branche principale en 2026 (C++17 avec
  CUDA 13.1 et plus, tuples de CUDA 13) : l'environnement utilise un
  commit de cette branche. Avec CUDA 13 sur sm_120, le multigrille
  classique échouait (« invalid device function ») : patch de la recette
  spack, reproducteur `reproducers/amgx-classical-cuda13`.

### Ginkgo

Multigrille algébrique par agrégation (`Pgm`), distribuée depuis 1.8. Pas
d'enveloppe MFEM pour la multigrille de Ginkgo (ni en 4.9 ni en 4.10) : à
configurer par l'API de Ginkgo. Les enveloppes de MFEM ne proposent que
les exécuteurs de référence, OpenMP, CUDA et HIP (pas SYCL). Voir aussi
[Plusieurs GPU](#plusieurs-gpu) et [Précision mixte](#précision-mixte).

### MueLu (Trilinos)

Agrégation lissée avec modes rigides sur GPU par Kokkos. Pas d'interface
MFEM : intégration lourde. Intérêt 3, à garder en tête comme point de
comparaison.

## Précision mixte

Sur la RTX 5080, le calcul en double précision est environ 64 fois plus
lent qu'en simple. Les produits matrice-vecteur sont limités par la
mémoire : en CSR, passer les valeurs en simple précision réduit les
octets lus d'environ un tiers (12 octets par coefficient au lieu de 8,
les indices restant sur 4 octets). Les calculs
denses (lois, petits blocs) sont limités par le calcul en double. Sur
A100 ou H100, le rapport n'est que de 2.

Idée commune : le Krylov extérieur et le résidu restent en double, le
préconditionneur travaille en simple (ou moins). La solution garde la
précision double.

- **Ginkgo** : Jacobi par blocs à précision adaptative (`storage_opt` et
  `accuracy` de `JacobiPreconditioner` dans l'enveloppe de MFEM 4.9,
  vérifié), CB-GMRES (base de Krylov compressée) et raffinement itératif
  (`CBGMRESSolver`, `IRSolver` dans l'enveloppe de MFEM 4.9), demi
  précision et bfloat16 (1.9, 1.10), matrices distribuées en précision
  mixte (1.10), multigrille en précision mixte.
- **hypre 3** : fonctions en précision mixte (3.0 : choix de la précision
  à l'exécution, Krylov préconditionné en précision plus basse), sur GPU
  depuis 3.2 (notes de version). Mais la recette spack de hypre interdit
  `precision=mixed` avec CUDA (`conflicts("precision=mixed")` sous
  `+cuda`), et MFEM exige un hypre en double (`hypre precision=double`
  dans sa recette) : il faudrait changer la recette, puis vérifier la
  compatibilité avec MFEM.
- **AmgX** : mode `dDFI` (voir plus haut).
- **MFEM** : pas de précision mixte. MFEM se compile entièrement en
  simple ou en double (`precision`).

Intérêt 6, 8 sur GPU grand public.

## Solveurs directs sur GPU

Intérêt pour les petits cas, les configurations mixtes et comme référence.

- **cuDSS** (intérêt 5, fait le 12/10) : solveur direct de NVIDIA, un ou
  plusieurs GPU, MPI. MFEM 4.10 l'enveloppe (`CuDSSSolver`, avec
  `SetReorderingReuse` : l'analyse symbolique est gardée quand seule la
  valeur des coefficients change, ce qui est le cas d'une itération de
  Newton à l'autre). La recette spack de MFEM 4.10 accepte cuDSS 0.6 à
  0.7.1 ; l'environnement installe cuDSS 0.7.1 pour CUDA 13 et MFEM avec
  `+cudss`. Solveur `CuDSSSolver` de la fabrique de mfem-mgis, avec
  `Symmetric`, `PositiveDefinite` et `ReuseAnalysis` (motif comparé sur
  tous les processus, solveur de MFEM recréé quand le motif ou la taille
  change : MFEM ne permet de choisir la réutilisation qu'avant la première
  matrice et s'arrête sur une taille différente). Sur la matrice
  d'élasticité d'ex8, une factorisation prend 0,018 s contre 0,203 s pour
  MUMPS sur l'hôte ; ex8 sur le GPU 9,8 -> 4,4 s, puis 2,6 s avec la
  thermique aussi par cuDSS (9 516 inconnues : factorisation 6 ms,
  résolution 1 ms, contre 52 ms par GMRES et BoomerAMG sur l'hôte)
  (`results/cudss-2026-10-12`). Chaque solution est contrôlée par son
  résidu (cuDSS ne signale pas les échecs de la factorisation). L'analyse
  distribuée de cuDSS 0.7.1 échoue pour certaines matrices que MUMPS
  résout : lignes réduites à leur coefficient diagonal, matrice
  tridiagonale couplant les processus, certaines répartitions inégales
  (`reproducers/cudss-distributed`, plusieurs processus sur un GPU) ; un
  processus n'est pas touché, et les matrices de diffusion passent.
- **MUMPS sur l'hôte, analyse réutilisée** (fait le 12/10). Dans ex8 sur
  le GPU, les 29 factorisations de MUMPS de la mécanique sont le premier
  poste. Le paramètre `ReuseAnalysis` du solveur `MUMPSSolver` de
  mfem-mgis garde l'analyse (ordonnancement et factorisation symbolique)
  tant que le motif de la matrice ne change sur aucun processus :
  `SetReorderingReuse` de MFEM ne vérifie pas le motif, mfem-mgis le
  compare avant chaque factorisation. Factorisations −19 % sur le GPU
  (un processus) et −26 % sur CPU (deux processus), ex8 sur le GPU
  11,1 -> 9,7 s, mêmes résultats (`results/mumps-analyse-2026-10-12`).
  Les autres réglages de MUMPS (ordonnancements, factorisation approchée
  BLR) ne gagnent rien sur cette matrice. La même comparaison du motif
  sert aussi au solveur cuDSS (ci-dessus) : avec `SetReorderingReuse`,
  MFEM n'y vérifie que le nombre de coefficients (et s'arrête s'il
  change), pas leurs colonnes.
- **SuperLU_DIST** (intérêt 4) : factorisation sur GPU (CUDA). MFEM 4.10
  ajoute `SuperLUSolver::SetDeviceOffload`.
- **STRUMPACK** (intérêt 4) : `STRUMPACKSolverBase::EnableGPU` dans MFEM
  4.9 (vérifié).

## Plusieurs GPU

Utile pour l'étude sur grappe (un processus par GPU).

- **Ginkgo distribué** (ajout, intérêt 7) : matrices et vecteurs
  distribués, multigrille distribuée avec grossissement PGM (1.8),
  Schwarz à deux niveaux (1.10), CG pipeliné (1.10), communications recouvertes par le calcul
  (1.11). Enveloppes MFEM pour les solveurs et préconditionneurs
  distribués dans MFEM 4.10 (exemple `examples/ginkgo/ex1p`, demande
  Ginkgo 1.11 compilé avec MPI, ce qu'a l'environnement). Avant 4.10,
  l'enveloppe ne prenait que des matrices séquentielles.
- **Krylov pipeliné** (ajout, intérêt 5) : une seule réduction globale par
  itération, recouverte par le calcul. CG pipeliné de Ginkgo, `KSPPIPECG`
  et voisins dans PETSc. Gain seulement à grand nombre de GPU.
- **MPI compatible GPU** (outil, intérêt 5) : échanges directement entre
  mémoires GPU. `Device::SetGPUAwareMPI` ou variable `MFEM_GPU_AWARE_MPI`
  (depuis MFEM 4.9) ; `GroupCommunicator` sur le périphérique (MFEM
  4.10) ; variante `+gpu-aware-mpi` de hypre (pas activée). L'OpenMPI de
  l'environnement est compilé avec CUDA (hwloc interne sans OpenCL), pour
  que le code reste générique ; sur une seule carte, l'effet est nul. Sur
  une grappe, le MPI de la grappe. Sans MPI compatible GPU, les échanges
  de la prolongation et ceux de hypre passent par l'hôte à chaque
  application : « aucune copie cachée » n'est alors
  atteignable qu'en un seul processus.

## Assemblage et opérateur tangent

Le plan prévoit trois options pour l'opérateur tangent
([plan-mfem-mgis.md](plan-mfem-mgis.md#choix-retenus)). Ajouts :

- **Règle d'intégration** (ajout, intérêt 8) : mfem-mgis intègre avec
  `IntRules.Get(geom, 2 * OrderGrad)` (code des intégrateurs générés),
  soit 27 points sur l'hexaèdre d'ordre 1. Une règle à 8 points divise
  environ par 3,4 le coût de la loi, du résidu et de la tangente, sur CPU
  comme sur GPU. C'est le levier de coût le plus direct des maillages
  hexaédriques. Demande une option de règle et la vérification de
  l'absence de modes parasites.
- **Assemblage de la matrice globale sur GPU** (ajout, intérêt 9) : pour
  l'option c sur GPU, les matrices élémentaires calculées sur GPU doivent
  donner la matrice hypre sans passer par l'hôte. Deux chemins :
  - le remplissage de matrice creuse de MFEM sur le périphérique
    (`ElementRestriction::FillJAndData`, niveau d'assemblage `FULL` des
    formes bilinéaires). Il travaille sur des degrés de liberté
    scalaires ; l'élasticité de MFEM l'utilise composante par composante.
    Pour un espace à 3 composantes, à vérifier, ou à faire par blocs de
    composantes ;
  - l'assemblage IJ de hypre sur GPU (pointeurs en mémoire GPU, depuis
    hypre 2.20).
  Le motif creux, la `HypreParMatrix` et son schéma de communication
  doivent être gardés d'une itération à l'autre (seules les valeurs
  changent), et les degrés de liberté imposés éliminés sur le
  périphérique (possible depuis MFEM 4.5). Sans cela, l'option c sur GPU
  recopie la tangente sur l'hôte, ou reconstruit sa structure, à chaque
  itération.
- **Assemblage par élément** (ajout, intérêt 6) : garder les matrices
  élémentaires et appliquer la tangente élément par élément (niveau
  `ELEMENT` de MFEM pour les formes bilinéaires). Intermédiaire entre
  l'option a (rien de stocké que `K` aux points) et l'option c (matrice
  globale) : plus de mémoire que a, pas d'assemblage global, diagonale
  et blocs nodaux immédiats. `BatchedLinAlg` de MFEM (4.8) pour les
  produits par lots ; graphes CUDA pour réduire le coût de lancement des
  petits noyaux (petits cas, plusieurs matériaux).
- La mini-application `hooke` de MFEM (élasticité non linéaire sans
  matrice sur GPU, Newton et CG, préconditionneur diagonal ou par blocs)
  est la référence la plus proche de l'option a dans MFEM.

## Problèmes couplés

- **MGR de hypre** (intérêt 4) : réduction multigrille pour les systèmes
  monolithiques (thermomécanique, micromorphe). Construction sur GPU
  depuis 2.24, sans mémoire unifiée depuis 2.29, nombreuses options depuis (3.2 : réutilisation de
  composants, complément de Schur fourni, plusieurs cycles).
- Préconditionneurs par blocs de MFEM (`BlockDiagonalPreconditioner`,
  complément de Schur) : génériques, tournent sur GPU si les blocs y
  sont.

## Outils de mesure

Indispensables aux études de performance, sans être des solveurs.

- **Profils** (intérêt 7) : Nsight Systems et Nsight Compute (installés
  sur la machine) ; marqueurs NVTX dans hypre (variante
  `+gpu-profiling`) ; Caliper (pas encore dans l'environnement). Sans
  marqueurs, un profil de hypre est difficile à lire.
- **Pools mémoire** (intérêt 6) : Umpire pour hypre (`+umpire`, dans
  l'environnement) et pour MFEM (`+umpire`). La construction d'un AMG sur
  GPU alloue beaucoup de petits tableaux. MFEM compilé avec Umpire prend
  les types de mémoire Umpire par défaut : à garder en tête dans les
  mesures.
- **Périphérique `debug` de MFEM** : il simule le GPU sur le CPU (noyaux
  exécutés sur l'hôte, mémoire « périphérique » distincte) et détecte les
  accès à une copie périmée. Il ne peut pas exécuter une loi compilée
  pour `cuda` : avec lui, utiliser une loi `openmp` ou `stlpar` et un
  emplacement accessible par l'hôte.
- GPU sans HMM et profil Nsight Systems : vérifier qu'aucune copie cachée
  n'a lieu sur le vrai GPU (déjà dans le plan).

## Écartés (CPU seul)

- Solveurs directs MUMPS et UMFPACK (restent dans l'environnement : ce sont
  les solveurs actuels de mfem-mgis, et la référence des configurations
  mixtes) ;
- `HypreEuclid`, `HypreParaSails` ; options élasticité de BoomerAMG
  (`SetElasticityOptions`, voir plus haut) ;
- BDDC et FETI (PETSc), HPDDM : support GPU partiel ou absent ;
- (n'est plus écarté) libCEED par MFEM : backend `ceed-cpu` et
  `ceed-cuda` de MFEM, installé pour comparaison avec l'assemblage partiel
  natif ; aussi utile à travers Ratel.

## MFEM 4.10

Tout le prototype utilise MFEM 4.10 (décision du 09/10/2026 ; mfem-mgis
était jusqu'ici utilisé avec 4.9). MFEM 4.10 (01/09/2026) apporte, pour
cette étude (notes de version et sources de 4.10) :

- enveloppe de cuDSS (`CuDSSSolver`) ;
- solveurs et préconditionneurs distribués de Ginkgo ;
- `SuperLUSolver::SetDeviceOffload` ;
- `GroupCommunicator` sur le périphérique, avec MPI compatible GPU ;
- noyaux d'assemblage partiel améliorés (diagonale de l'élasticité) ;
- assemblage LOR par lots sur maillages très connectés.

Les faits vérifiés sur les sources de 4.9 et cités dans ce document
(BoomerAMG sur GPU, L-BFGS, `mfem::forall`, AmgX, Ginkgo, assemblage
partiel des formes non linéaires) ont été revérifiés sur celles de 4.10.

## Correspondance avec l'environnement spack

| Technique | Paquet et variante de `mmm-gpu` |
|---|---|
| Krylov de MFEM, lisseurs sans matrice, p-MG, LOR, `BatchedLinAlg` | `mfem@4.10.0 +cuda cxxstd=20` |
| Mêmes noyaux sur CPU multicœur, backend libCEED | `mfem +openmp +threadsafe +libceed`, `libceed +cuda +openmp` |
| Newton inexact, recherche linéaire, point fixe et Anderson, L-BFGS préconditionné | à écrire dans mfem-mgis |
| BoomerAMG, MGR, FSAI, ILU, Krylov de hypre | `hypre@3.2.0 +cuda +umpire` |
| Ginkgo, distribué compris | `ginkgo@1.11.0 +cuda +openmp +mpi`, `mfem +ginkgo` |
| GAMG, SNES, Krylov pipeliné | `petsc@3.26.0 +cuda +kokkos +openmp +hypre`, `kokkos-kernels`, `mfem +petsc` |
| KINSOL | `sundials +cuda +mpi`, `mfem +sundials` |
| AmgX | `amgx` (commit compatible CUDA 13), `mfem +amgx` |
| cuDSS | `cudss@0.7.1-13 +mpi`, `mfem +cudss` |
| SuperLU_DIST, STRUMPACK sur GPU | `superlu-dist +cuda`, `strumpack +cuda`, `mfem +superlu-dist +strumpack` |
| MPI compatible GPU | `openmpi@5 +cuda` (hwloc interne sans OpenCL) |
| Pools mémoire | `umpire +cuda`, `mfem +umpire`, `hypre +umpire` |
| Solveurs directs CPU de mfem-mgis | `mfem +mumps +suite-sparse` |

Pas encore dans l'environnement : Caliper, hypre `+gpu-profiling`, hypre
en précision mixte, MueLu, Ratel (qui demande PETSc 3.24).

## Sources

- Sources de MFEM 4.9.0 et 4.10.0 (archives de la recette spack, sommes de
  contrôle vérifiées) : `linalg/hypre.cpp` (`SetDefaultOptions`,
  `SetSystemsOptions`, `SetElasticityOptions`), `linalg/ginkgo.hpp`,
  `linalg/amgxsolver.hpp` et `.cpp`, `linalg/sundials.hpp`,
  `linalg/strumpack.hpp`, `linalg/solvers.hpp`, `linalg/cudss.hpp` (4.10),
  `miniapps/hooke`, `miniapps/diag-smoothers`, `CHANGELOG`, `INSTALL`.
- hypre : `CHANGELOG` de la version 3.2.0, `src/parcsr_ls/HYPRE_parcsr_ls.h`
  (3.2.0), documentation de BoomerAMG (options GPU).
- Ginkgo : `CHANGELOG.md` de la version 1.11.0.
- PETSc 3.24.6 et 3.26.0 : `config/BuildSystem/config/packages/kokkos.py`
  et `kokkos-kernels.py`.
- AmgX : historique de la branche principale (corrections CUDA 13 de
  janvier et juillet 2026, après la version 2.5.0 de décembre 2025).
- mfem-mgis : `src/LinearSolverFactory.cxx` (stratégies de BoomerAMG),
  `src/NewtonSolver.cxx`, `src/NonLinearSolverBase.cxx`,
  `src/MultiMaterialNonLinearIntegrator.cxx`, intégrateurs générés.
- Relecture indépendante de ce document le 09/10/2026 (constats corrigés).
- Recettes spack du dépôt `spack-packages` au commit `3b7fdffc24`.
- Ratel : `https://gitlab.com/micromorph/ratel`, arXiv 2204.01722.
