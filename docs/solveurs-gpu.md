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
13. [MFEM 4.9 ou 4.10](#mfem-49-ou-410)
14. [Correspondance avec l'environnement spack](#correspondance-avec-lenvironnement-spack)
15. [Sources](#sources)

## Tableau d'ensemble

Les lignes sans marque viennent de la liste de départ. Les lignes marquées
« ajout » complètent cette liste. Ginkgo est déjà retenu : il apparaît dans
plusieurs lignes.

| Technique | CPU/GPU | Intérêt | Où |
|---|---|---|---|
| Réutiliser le préconditionneur (AMG construit une fois, sur l'opérateur élastique, ou gardé plusieurs itérations ou pas de temps) | les deux | 9 | mfem-mgis |
| Ajout : réutilisation partielle de la hiérarchie AMG (prolongateurs gardés, opérateurs grossiers recalculés) | les deux | 7 | GAMG de PETSc, AmgX (`UpdateOperator` de MFEM) |
| Point fixe accéléré par Anderson | les deux | 9 | mfem-mgis, KINSOL |
| Ajout : assemblage de la matrice tangente globale sur GPU (option c sur GPU) | les deux | 9 | mfem-mgis, hypre (IJ sur GPU), MFEM |
| Newton inexact (Eisenstat-Walker, `NewtonSolver::SetAdaptiveLinRtol`) | les deux | 8 | MFEM |
| Lisseurs sans matrice Jacobi et Chebyshev sur la diagonale (`OperatorJacobiSmoother`, `OperatorChebyshevSmoother`), indispensables à l'option a | les deux | 8 | MFEM |
| Ajout : Jacobi par blocs nodaux 3×3 sans matrice (diagonale par blocs de la tangente) | les deux | 7 | mfem-mgis, comme la mini-application `hooke` de MFEM |
| Ajout : l1-Jacobi « en valeur absolue » sans matrice | les deux | 5 | MFEM 4.9, `miniapps/diag-smoothers` |
| Réglages GPU de BoomerAMG exposés dans mfem-mgis (PMIS, ext+i, l1-Jacobi, Chebyshev, nombre de niveaux, seuils) | les deux | 8 | hypre, MFEM |
| Ajout : stratégie « System » de BoomerAMG sur GPU (la stratégie « Elasticity » est refusée par MFEM avec hypre GPU) | les deux | 8 | MFEM, mfem-mgis |
| Précision mixte : préconditionneur en simple précision, Krylov en double ; gain fort sur la RTX 5080 (FP64 lent), plus faible sur A100 ou H100 | surtout GPU | 6 (8 sur GPU grand public) | Ginkgo, hypre 3, AmgX |
| L-BFGS préconditionné par un cycle AMG ou p-MG (`LBFGSSolver` de MFEM), comme Ratel | les deux | 7 | MFEM |
| p-multigrille sans matrice avec AMG sur le Q1 (ordre 2 et plus), ou LOR (`LORSolver`) | les deux | 7 | MFEM |
| Ajout : Ginkgo distribué (matrices MPI, Schwarz à deux niveaux, multigrille PGM distribuée) | les deux | 7 | Ginkgo 1.10 et 1.11 ; enveloppes MFEM 4.10 |
| Ajout : AmgX en précision mixte (mode `dDFI`) et matrices par blocs 3×3 | GPU NVIDIA | 6 | AmgX ; pas dans l'enveloppe de MFEM |
| Recherche linéaire pour Newton (crochet `ComputeScalingFactor`) | les deux | 6 | MFEM |
| AMG élasticité avec modes rigides sur GPU : GAMG de PETSc avec espace quasi nul (Kokkos), nouveau chemin par blocs sur GPU (Adams, arXiv 2606.24748) ; demande de réparer le chemin PETSc de mfem-mgis | les deux | 6 | PETSc 3.26 |
| KINSOL de SUNDIALS par MFEM (`KINSolver` : point fixe avec Anderson, JFNK, recherche linéaire) | les deux | 6 | SUNDIALS 7.9, MFEM |
| Ajout : prédiction par extrapolation de l'incrément de déplacement | les deux | 5 | mfem-mgis |
| AmgX par `AmgXSolver` de MFEM (point de comparaison avec BoomerAMG, modes MPI) | GPU NVIDIA | 5 | AmgX, MFEM |
| cuDSS (solveur direct NVIDIA, un ou plusieurs GPU, MPI, refactorisation) | GPU NVIDIA | 5 | cuDSS ; enveloppe dans MFEM 4.10, pas dans 4.9 |
| Ajout : Krylov pipeliné (CG pipeliné de Ginkgo, `KSPPIPECG` de PETSc) | les deux | 5 | Ginkgo 1.10, PETSc |
| Ajout : Ratel comme calcul de référence sur le même cas | les deux | 5 | Ratel, PETSc, libCEED |
| SuperLU_DIST et STRUMPACK sur GPU (enveloppes MFEM existantes) | les deux | 4 | MFEM |
| MGR de hypre (GPU depuis 2.29) pour les problèmes couplés monolithiques | les deux | 4 | hypre |
| FSAI (GPU depuis 2.30) et ILU de hypre sur GPU | les deux | 4 | hypre |
| Ajout : hypre récent (filtrage des couplages entre composantes, cycle flexible, Schwarz avec recouvrement, solveur de Chebyshev) | les deux | 4 | hypre 2.32 à 3.2 |
| Ajout : pas de temps adaptatif guidé par la convergence de Newton | les deux | 4 | mfem-mgis |
| Ajout : solveurs non linéaires de PETSc sur vecteurs GPU (NGMRES, Anderson, quasi-Newton, NASM) | les deux | 4 | PETSc SNES |
| Ajout : MueLu de Trilinos (agrégation lissée avec modes rigides, Kokkos) | les deux | 3 | Trilinos ; pas d'interface MFEM |
| Ajout : recyclage de sous-espaces de Krylov entre itérations de Newton | les deux | 3 | PETSc (HPDDM, GCRO-DR) ; GPU à vérifier |
| Ajout : lois de comportement en simple précision sur GPU grand public | GPU | 3 | demande TFEL et MGIS |
| Ajout (outils) : profils NVTX de hypre, Caliper, Nsight Systems | les deux | 7 | hypre `+gpu-profiling`, MFEM |
| Ajout (outils) : pools mémoire Umpire pour MFEM et hypre | surtout GPU | 6 | Umpire |
| Ajout (outils) : MPI compatible GPU | GPU | 5 | MPI de la grappe, MFEM, hypre |

## Réutiliser les préconditionneurs

Sur GPU, la construction d'un AMG coûte souvent plus que plusieurs
résolutions. mfem-mgis la refait à chaque itération de Newton
(`prec->SetOperator(J)`).

- **Réutilisation complète** : construire l'AMG une fois et le garder.
  Variantes : sur l'opérateur élastique (garde toute la simulation), sur
  la tangente du début du pas (garde un pas), sur la tangente gardée
  `n` itérations. Rien à ajouter aux bibliothèques : c'est un choix de
  mfem-mgis (ne pas rappeler `SetOperator` du préconditionneur quand
  l'opérateur change). Le Krylov reste exact, seul le préconditionneur
  vieillit. Intérêt 9.
- **Réutilisation partielle** (ajout) : garder la structure de la
  hiérarchie (agrégats, prolongateurs) et recalculer seulement les
  opérateurs grossiers `Pᵀ A P` avec les nouvelles valeurs. Coût de
  construction bien plus faible, qualité proche d'une construction
  complète tant que la tangente ne change pas trop. Existe dans :
  - GAMG de PETSc (`PCGAMGSetReuseInterpolation`) ;
  - AmgX : l'enveloppe de MFEM le fait déjà (`AmgXSolver::UpdateOperator`,
    `AMGX_matrix_replace_coefficients` puis nouvelle construction,
    vérifié dans `linalg/amgxsolver.cpp` de MFEM 4.9) ;
  - MGR de hypre 3.2 (« component reuse », notes de version) ;
  - BoomerAMG : pas d'équivalent trouvé.
  Intérêt 7.

## Solveurs non linéaires

Tous ces solveurs n'utilisent que des vecteurs, des actions d'opérateurs
et des produits scalaires : ils tournent sur GPU dès que les vecteurs et
les opérateurs y sont.

- **Point fixe accéléré par Anderson** (intérêt 9). À écrire dans
  mfem-mgis (le solveur de point fixe du plan, option b), ou par KINSOL
  (`KINSolver::EnableAndersonAcc`, pour `KIN_FP` et `KIN_PICARD`, vérifié
  dans `linalg/sundials.hpp` de MFEM 4.9). Le point fixe seul convergeait
  mal sur le polycristal des essais préalables : Anderson est la première
  chose à essayer.
- **Newton inexact** (intérêt 8) : `NewtonSolver::SetAdaptiveLinRtol`
  (vérifié dans MFEM 4.9). Une option de mfem-mgis suffit.
- **L-BFGS** (intérêt 7) : `LBFGSSolver` de MFEM (vérifié dans 4.9), avec
  un cycle d'AMG ou de p-multigrille comme inverse initial, comme Ratel.
  Demande un préconditionneur réglable dans `LBFGSSolver`, à vérifier.
- **Recherche linéaire** (intérêt 6) : redéfinir
  `NewtonSolver::ComputeScalingFactor` (vérifié dans 4.9).
- **KINSOL** (intérêt 6) : point fixe avec Anderson, Picard, Newton-Krylov
  sans matrice (`SetJFNK`), recherche linéaire. L'enveloppe de MFEM
  utilise les vecteurs CUDA de SUNDIALS quand MFEM est compilé avec CUDA
  (vérifié dans `linalg/sundials.hpp`). À essayer sur la machine.
- **Prédiction par extrapolation** (ajout, intérêt 5) : début de chaque pas
  à partir de l'incrément du pas précédent (linéaire), ou des deux
  précédents (quadratique). Moins chère que la prédiction par l'opérateur
  tangent de mfem-mgis : aucune résolution. Réduit le nombre
  d'itérations de Newton quand le chargement est régulier.
- **Pas de temps adaptatif** (ajout, intérêt 4) : réduire le pas quand
  Newton converge mal, l'augmenter quand il converge vite. Influence
  directement le coût total, donc la comparaison des variantes.
- **Solveurs non linéaires de PETSc** (ajout, intérêt 4) : NGMRES,
  Anderson, quasi-Newton (L-BFGS), NASM, avec des vecteurs sur GPU.
  Utiles comme implémentations de référence pour les points ci-dessus,
  si le chemin PETSc de mfem-mgis est réparé.

## Préconditionneurs sans matrice

Nécessaires à l'option a (tangente sans matrice) : l'AMG demande une
matrice, les lisseurs ci-dessous seulement la diagonale ou des blocs.

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
  modèle direct pour mfem-mgis.
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
- **Stratégies de mfem-mgis sur GPU** (ajout, intérêt 8) : la stratégie
  « Elasticity » appelle `HypreBoomerAMG::SetElasticityOptions`, qui
  arrête le programme avec hypre GPU (« this method is not supported in
  hypre built with GPU support », vérifié dans MFEM 4.9). Sur GPU,
  mfem-mgis doit refuser cette stratégie avec un message clair ou passer à
  « System » (`SetSystemsOptions`, approche par inconnue, qui fonctionne
  sur GPU : MFEM copie la table des composantes sur le GPU). Les modes
  rigides (vecteurs d'interpolation) restent sur CPU avec BoomerAMG.
- **Filtrage des couplages entre composantes** (ajout, intérêt 4) :
  `HYPRE_BoomerAMGSetFilterFunctions` (hypre 2.32) retire les couplages
  entre composantes de la matrice utilisée pour l'AMG : complexité plus
  faible, convergence à mesurer en élasticité. GPU à vérifier.

### GAMG (PETSc)

Agrégation lissée avec espace quasi nul (modes rigides) : l'AMG de
référence pour l'élasticité. Avec PETSc compilé avec Kokkos ou CUDA, les
matrices sont sur GPU (`aijkokkos`, `aijcusparse`) ; la part de la
construction faite sur GPU dépend de la version, à mesurer. PETSc 3.26
utilise Kokkos 5.2.2 par défaut, la version de l'environnement. Le chemin
PETSc de mfem-mgis plantait lors d'essais antérieurs : à réparer d'abord.
MFEM 4.10 exige PETSc compilé avec hypre. Intérêt 6.

### AmgX (NVIDIA)

- `AmgXSolver` de MFEM : matrice scalaire, double précision seulement
  (`AMGX_mode_dDDI`, commentaire « currently only double is supported »
  dans `linalg/amgxsolver.hpp`), modes MPI, réutilisation partielle par
  `UpdateOperator`. Intérêt 5.
- **Au-delà de l'enveloppe** (ajout, intérêt 6) : AmgX accepte des
  matrices par blocs (blocs 3×3 nodaux en élasticité) et le mode
  `dDFI` (vecteurs en double, matrice et préconditionneur en simple). Il
  faudrait une enveloppe propre à mfem-mgis.
- La version 2.5.0 d'AmgX précède les corrections pour CUDA 13 (faites sur
  la branche principale en 2026) : l'environnement utilise un commit de
  cette branche.

### Ginkgo

Multigrille algébrique par agrégation (`Pgm`), distribuée depuis 1.10. Pas
d'enveloppe MFEM pour la multigrille de Ginkgo (ni en 4.9 ni en 4.10) : à
configurer par l'API de Ginkgo. Voir aussi
[Plusieurs GPU](#plusieurs-gpu) et [Précision mixte](#précision-mixte).

### MueLu (Trilinos)

Agrégation lissée avec modes rigides sur GPU par Kokkos. Pas d'interface
MFEM : intégration lourde. Intérêt 3, à garder en tête comme point de
comparaison.

## Précision mixte

Sur la RTX 5080, le calcul en double précision est environ 64 fois plus
lent qu'en simple. Les produits matrice-vecteur sont limités par la
mémoire : la simple précision divise par deux les octets lus. Les calculs
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
  depuis 3.2 (notes de version). Spack propose `precision=mixed`, mais
  MFEM exige un hypre en double (`hypre precision=double` dans sa
  recette) : compatibilité avec MFEM à vérifier, probablement une
  installation séparée.
- **AmgX** : mode `dDFI` (voir plus haut).
- **MFEM** : pas de précision mixte. MFEM se compile entièrement en
  simple ou en double (`precision`).

Intérêt 6, 8 sur GPU grand public.

## Solveurs directs sur GPU

Intérêt pour les petits cas, les configurations mixtes et comme référence.

- **cuDSS** (intérêt 5) : solveur direct de NVIDIA, un ou plusieurs GPU,
  MPI. MFEM 4.10 l'enveloppe (`CuDSSSolver`, avec
  `SetReorderingReuse` : l'analyse symbolique est gardée quand seule la
  valeur des coefficients change, ce qui est le cas d'une itération de
  Newton à l'autre). MFEM 4.9 n'a pas d'enveloppe. La recette spack de
  MFEM 4.10 accepte cuDSS 0.6 à 0.7.1 ; l'environnement installe cuDSS
  0.7.1 pour CUDA 13.
- **SuperLU_DIST** (intérêt 4) : factorisation sur GPU (CUDA). MFEM 4.10
  ajoute `SuperLUSolver::SetDeviceOffload`.
- **STRUMPACK** (intérêt 4) : `STRUMPACKSolverBase::EnableGPU` dans MFEM
  4.9 (vérifié).

## Plusieurs GPU

Utile pour l'étude sur grappe (un processus par GPU).

- **Ginkgo distribué** (ajout, intérêt 7) : matrices et vecteurs
  distribués, Schwarz à deux niveaux (1.10), multigrille distribuée
  (1.10), CG pipeliné (1.10), communications recouvertes par le calcul
  (1.11). Enveloppes MFEM pour les solveurs et préconditionneurs
  distribués dans MFEM 4.10 seulement (exemple `examples/ginkgo/ex1p`,
  demande Ginkgo 1.11 compilé avec MPI). En MFEM 4.9, l'enveloppe ne
  prend que des matrices séquentielles : utilisable comme préconditionneur
  local à chaque processus (Jacobi par blocs de processus).
- **Krylov pipeliné** (ajout, intérêt 5) : une seule réduction globale par
  itération, recouverte par le calcul. CG pipeliné de Ginkgo, `KSPPIPECG`
  et voisins dans PETSc. Gain seulement à grand nombre de GPU.
- **MPI compatible GPU** (outil, intérêt 5) : échanges directement entre
  mémoires GPU. `Device::SetGPUAwareMPI` ou variable `MFEM_GPU_AWARE_MPI`
  (MFEM 4.9) ; variante `+gpu-aware-mpi` de hypre. Sur une seule carte,
  sans effet : l'OpenMPI de l'environnement est compilé sans CUDA (avec
  CUDA, il demande un hwloc avec CUDA, dont le changement reconstruirait
  LLVM et TBB, utilisés par TFEL). PETSc doit alors être lancé avec
  `-use_gpu_aware_mpi 0`. Sur une grappe, le MPI de la grappe.

## Assemblage et opérateur tangent

Le plan prévoit trois options pour l'opérateur tangent
([plan-mfem-mgis.md](plan-mfem-mgis.md#choix-retenus)). Deux ajouts :

- **Assemblage de la matrice globale sur GPU** (ajout, intérêt 9) : pour
  l'option c sur GPU, les matrices élémentaires calculées sur GPU doivent
  donner la matrice hypre sans passer par l'hôte. Deux chemins : le
  remplissage de matrice creuse de MFEM sur le périphérique (utilisé par
  le niveau d'assemblage `FULL` des formes bilinéaires), ou l'assemblage
  IJ de hypre sur GPU (pointeurs en mémoire GPU, depuis hypre 2.20). Sans
  cela, l'option c sur GPU recopie la tangente sur l'hôte à chaque
  itération.
- **Assemblage par élément** (ajout, intérêt 6) : garder les matrices
  élémentaires et appliquer la tangente élément par élément (niveau
  `ELEMENT` de MFEM pour les formes bilinéaires). Intermédiaire entre
  l'option a (rien de stocké que `K` aux points) et l'option c (matrice
  globale) : plus de mémoire que a, pas d'assemblage global, diagonale
  et blocs nodaux immédiats.
- La mini-application `hooke` de MFEM (élasticité non linéaire sans
  matrice sur GPU, Newton et CG, préconditionneur diagonal ou par blocs)
  est la référence la plus proche de l'option a dans MFEM.

## Problèmes couplés

- **MGR de hypre** (intérêt 4) : réduction multigrille pour les systèmes
  monolithiques (thermomécanique, micromorphe). Construction sur GPU
  depuis 2.29, nombreuses options depuis (3.2 : réutilisation de
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
  GPU alloue beaucoup de petits tableaux.
- **Périphérique `debug` de MFEM**, GPU sans HMM : vérifier qu'aucune
  copie cachée n'a lieu (déjà dans le plan).

## Écartés (CPU seul)

- Solveurs directs MUMPS et UMFPACK (restent dans l'environnement : ce sont
  les solveurs actuels de mfem-mgis, et la référence des configurations
  mixtes) ;
- `HypreEuclid`, `HypreParaSails` ; options élasticité de BoomerAMG
  (`SetElasticityOptions`, voir plus haut) ;
- BDDC et FETI (PETSc), HPDDM : support GPU partiel ou absent ;
- libCEED par MFEM : l'assemblage partiel de MFEM couvre le besoin ;
  libCEED reste utile à travers Ratel.

## MFEM 4.9 ou 4.10

Le plan retient MFEM 4.9, la version de mfem-mgis. MFEM 4.10 (01/09/2026)
apporte, pour cette étude (notes de version et sources de 4.10) :

- enveloppe de cuDSS (`CuDSSSolver`) ;
- solveurs et préconditionneurs distribués de Ginkgo ;
- `SuperLUSolver::SetDeviceOffload` ;
- `GroupCommunicator` sur le périphérique, avec MPI compatible GPU ;
- noyaux d'assemblage partiel améliorés (diagonale de l'élasticité) ;
- assemblage LOR par lots sur maillages très connectés.

Proposition : garder 4.9 dans `mmm-gpu` (référence de mfem-mgis), et
ajouter plus tard un second environnement avec MFEM 4.10 (mêmes
dépendances, cuDSS et Ginkgo distribué), pour ces comparaisons. À
décider.

## Correspondance avec l'environnement spack

| Technique | Paquet et variante de `mmm-gpu` |
|---|---|
| Krylov de MFEM, lisseurs sans matrice, Newton, L-BFGS, p-MG, LOR | `mfem@4.9.0 +cuda` |
| BoomerAMG, MGR, FSAI, ILU, Krylov de hypre | `hypre@3.2.0 +cuda +umpire` |
| Ginkgo (local en MFEM 4.9, distribué ensuite) | `ginkgo@1.11.0 +cuda +openmp +mpi`, `mfem +ginkgo` |
| GAMG, SNES, Krylov pipeliné | `petsc@3.26.0 +cuda +kokkos +hypre`, `kokkos-kernels`, `mfem +petsc` |
| KINSOL | `sundials +cuda +mpi`, `mfem +sundials` |
| AmgX | `amgx` (commit compatible CUDA 13), `mfem +amgx` |
| cuDSS | `cudss@0.7.1-13 +mpi` (pas d'enveloppe en MFEM 4.9) |
| SuperLU_DIST, STRUMPACK sur GPU | `superlu-dist +cuda`, `strumpack +cuda`, `mfem +superlu-dist +strumpack` |
| MPI | `openmpi@5` sans CUDA (voir [Plusieurs GPU](#plusieurs-gpu)) |
| Pools mémoire | `umpire +cuda`, `mfem +umpire`, `hypre +umpire` |
| Solveurs directs CPU de mfem-mgis | `mfem +mumps +suite-sparse` |

Pas encore dans l'environnement : Caliper, hypre `+gpu-profiling`, hypre
en précision mixte, MueLu, Ratel (qui demande PETSc 3.24), MFEM 4.10, MPI
compatible GPU.

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
- mfem-mgis : `src/LinearSolverFactory.cxx` (stratégies de BoomerAMG).
- Recettes spack du dépôt `spack-packages` au commit `3b7fdffc24`.
- Ratel : `https://gitlab.com/micromorph/ratel`, arXiv 2204.01722.
