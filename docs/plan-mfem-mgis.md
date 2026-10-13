# Plan pour mfem-mgis

Statut : choix retenus le 08/10/2026, première étape (API de vues de MGIS)
faite le 09/10/2026, plan détaillé le 09/10/2026
([plan-mfem-mgis-detaille.md](plan-mfem-mgis-detaille.md)).
Les lots L0 à L10 sont faits (base master `a011d97d6f`, branche rebasée
le 13/10/2026 sur master `8fac67336e`, dépôt `mfem-mgis/` ; état dans
[vue-ensemble.md](vue-ensemble.md#état-davancement)). Ce plan s'appuie sur
TFEL et MGIS tels que décrits dans [tfel.md](tfel.md) et
[mgis.md](mgis.md).

## Sommaire

1. [Esprit](#esprit)
2. [Point de départ](#point-de-départ)
3. [Choix retenus](#choix-retenus)
4. [Architecture cible](#architecture-cible)
5. [Organisation du code](#organisation-du-code)
6. [Portabilité](#portabilité)
7. [Solveurs](#solveurs)
8. [Configurations mixtes](#configurations-mixtes)
9. [Performance](#performance)
10. [Validation](#validation)
11. [Étapes](#étapes)
12. [Questions ouvertes](#questions-ouvertes)
13. [Risques](#risques)

## Esprit

Le travail sur mfem-mgis a deux buts, compatibles :

- **explorer** : mesurer, comparer, essayer beaucoup de variantes ;
- **proposer** aux développeurs de mfem-mgis un code propre, complet,
  flexible et modulaire, dans l'esprit de TFEL et de MGIS.

Donc tout est choisissable et rien n'est exclu. Le chemin actuel de
mfem-mgis reste le chemin par défaut. Les nouveaux chemins s'activent par
des options. Quand plusieurs façons de faire existent, elles sont toutes
implémentées, choisies par une option, et comparées sur les mêmes cas.
L'API doit être prête pour les expériences : mieux vaut prévoir trop que
pas assez, tant que le code reste propre.

## Point de départ

### mfem-mgis aujourd'hui

- Assemblage élément par élément sur l'hôte. Un intégrateur MFEM par
  matériau (`MultiMaterialNonLinearIntegrator`, `AssembleElementVector`,
  `AssembleElementGrad`), puis une classe par hypothèse, cinématique et
  symétrie (`include/MFEMMGIS/*BehaviourIntegrator.hxx` : petites et
  grandes déformations, contraintes et déformations planes, 3D, thermique
  non linéaire, micromorphe, FBar, régularisations...).
- Dans chaque élément : gradients aux points de Gauss, intégration de la
  loi point par point
  (`BehaviourIntegratorBase::performsLocalBehaviourIntegration`), puis
  contributions au résidu et à la matrice tangente. Ce chemin demande une
  bibliothèque `generic` : MGIS refuse une bibliothèque de lots dessus,
  avec un message clair.
- Données aux points de Gauss : `PartialQuadratureFunction` (stockage
  `std::vector` ou vue externe). Le matériau (`Material`) dérive du
  `MaterialDataManager` de MGIS.
- Le paramètre `FiniteElementSpaceOrdering` (`byNODES` par défaut, ou
  `byVDIM`) choisit l'ordre des **degrés de liberté** de l'espace éléments
  finis. C'est une autre question que le rangement des **données aux points
  de Gauss** (voir [Choix retenus](#choix-retenus), point 2).
- Solveur non linéaire : `Newton` seulement, propre à mfem-mgis (il ne
  dérive pas de `mfem::NewtonSolver`). Il donne la tangente au solveur
  linéaire à chaque itération (`prec->SetOperator(J)`, où `prec` est le
  solveur linéaire), ce qui reconstruit le préconditionneur à chaque
  itération. Pas de Newton inexact, pas de quasi-Newton, pas de point
  fixe ; la correction n'est divisée par 2 que si la loi échoue (pas de
  vraie recherche linéaire). La loi est intégrée par une action du
  solveur, pas par le `Mult` du formulaire : un solveur non linéaire
  extérieur ne l'intégrerait pas.
- Solveurs linéaires : voir [Solveurs](#solveurs). Ginkgo n'est pas
  branché.
- Aucun usage de l'assemblage partiel de MFEM (`AssemblyLevel`), de
  `mfem::Device`, `ElementRestriction`, `QuadratureInterpolator` ni
  `mfem::forall`.

### Ce que MGIS offre désormais

- Lots de lois, de modèles et de propriétés matériau, sur les tableaux de
  l'appelant, sans copie.
- Basculement automatique des fonctions de MGIS qui travaillent sur un
  `MaterialDataManager`.
- Politique d'exécution (file), emplacement des tableaux vérifié,
  ressources mémoire.
- Calcul synchrone par défaut, ou asynchrone : plusieurs matériaux peuvent
  être soumis avant une seule synchronisation.
- Pas encore : rangement des données composante par composante, mise à
  jour des états sur GPU (`update`, `revert`). Les deux sont la première
  étape de ce plan.

### Essais préalables

Des essais faits hors de ce dépôt (dossier `mfem-mgis-gpu`) ont comparé, sur
une poutre en plasticité et 8 processus :

- Newton avec PCG et AMG : 14,4 s sur CPU ;
- point fixe avec l'opérateur élastique et un cycle AMG : 8,1 s sur CPU ;
- sur GPU, avec une loi J2 écrite à la main en CUDA : 0,69 s, 0,34 s et
  0,15 s pour trois variantes.

Sur un polycristal, le point fixe converge mal : Newton avec GMRES et AMG
prend 22,6 s, le point fixe 67 s. Aucune méthode ne gagne partout.

## Choix retenus

1. **Point d'entrée dans MFEM : les deux.** L'assemblage partiel de
   `NonlinearForm` (`AssemblyLevel::PARTIAL` : `AddMultPA`,
   `AssembleGradPA`, `AddMultGradPA`, `AssembleGradDiagonalPA`) et un
   opérateur propre à mfem-mgis. Choix par option, sans dupliquer le cœur :
   les noyaux par lot (gradients, intégration par MGIS, résidu, action de
   la tangente) sont écrits une fois et servent aux deux.
2. **Rangement des données aux points de Gauss : les deux.** Point par
   point (composantes d'un point contiguës) et composante par composante
   (valeurs d'une composante contiguës, accès coalescés sur GPU). TFEL
   accepte déjà les deux. MGIS n'accepte aujourd'hui que le point par
   point : il lui faut une API sur des vues, conçue en premier puisque le
   reste en dépend. L'ordre des degrés de liberté
   (`FiniteElementSpaceOrdering`) est un réglage indépendant, déjà
   disponible.
3. **Opérateur tangent : trois options, choisies par une option.**
   - a. sans matrice : action `Σ Bᵀ K B` aux points, solveur de Krylov
     préconditionné par l'AMG de l'opérateur élastique ;
   - b. point fixe avec l'opérateur élastique : résidu seul à chaque
     itération ;
   - c. tangente assemblée : matrices élémentaires, puis matrice globale
     hypre.

   Il faut aussi un solveur non linéaire de point fixe. Le plan détaillé
   ajoute une variante de l'option a, l'assemblage par élément (matrices
   élémentaires gardées), et précise les combinaisons de paramètres de
   chaque option.
4. **Hypothèses, cinématiques, intégrateurs** : au final, tout mfem-mgis.
   On commence par les petites déformations 3D. Chaque étape reste testable
   sur CPU et sur GPU.
5. **Propriétés matériau et variables externes** : toutes les façons, au
   final (uniformes, tableaux calculés par MFEM, propriétés
   `generic-parallel` évaluées par MGIS).
6. **Plusieurs matériaux** : traités, avec le calcul asynchrone de MGIS.
7. **Sorties, post-traitements, reprises** : fonctionnement actuel par
   défaut (copie vers l'hôte).
8. **Version de MFEM : 4.10** pour tout le prototype (décision du
   09/10/2026 ; mfem-mgis était utilisé avec 4.9). MFEM 4.10 apporte
   cuDSS, Ginkgo distribué et `GroupCommunicator` sur le périphérique
   (voir [solveurs-gpu.md](solveurs-gpu.md#mfem-410)). Compilé avec CUDA,
   Ginkgo, PETSc, SUNDIALS, AmgX, cuDSS, SuperLU_DIST, STRUMPACK, MUMPS et
   SuiteSparse dans `deps/spack/mmm-gpu`, avec un OpenMPI compatible GPU.
9. **Mesures fines plus tard.** Priorité à un cadre propre qui marche.
   Les premières mesures servent à comparer les variantes.
10. **Solveurs** : choisir les combinaisons qui marchent sur GPU (Krylov de
    MFEM ou de hypre, BoomerAMG sur GPU, Jacobi), brancher Ginkgo dans la
    fabrique de solveurs (CPU et GPU, choisis par l'exécuteur), garder les
    solveurs directs sur CPU pour les configurations mixtes.
11. **mfem-mgis-examples et mm-opera-hpc** : suivis dans ce dépôt depuis
    leur base. Chaque exemple recevra des paramètres pour choisir CPU ou
    GPU, et les autres options (niveau d'assemblage, opérateur tangent,
    solveur, rangement des données).

## Architecture cible

À chaque itération, pour chaque matériau, tout reste sur le GPU :

1. **Gradients aux points de Gauss.** Les valeurs élémentaires sont
   extraites de la solution (`ElementRestriction` de MFEM), dérivées aux
   points (`QuadratureInterpolator`), puis converties en gradients de la loi
   (déformation, gradient de la transformation, gradient de température),
   écrits directement dans les gradients de fin de pas du matériau, dans
   son rangement.
2. **Intégration du lot** par MGIS, sur ces mêmes tableaux.
3. **Résidu** : `Bᵀ σ w det J` par point, somme par élément, puis
   assemblage parallèle par MFEM sur le GPU.
4. **Opérateur tangent** : une des trois options (a, b ou c).
5. **Solveurs linéaires** : hypre sur GPU (PCG, GMRES, BoomerAMG), Ginkgo
   en option (CUDA, HIP, OpenMP).
6. **Fin de pas** : `update` et `revert` de MGIS sur le GPU. mfem-mgis
   possède la mémoire, MGIS sait quels tableaux copier. Les transferts vers
   l'hôte se limitent aux sorties.
7. **Post-traitements** : sur l'hôte d'abord, sur le GPU ensuite.

## Organisation du code

- **Le chemin actuel ne change pas** (`AssemblyLevel::LEGACY`, élément par
  élément sur l'hôte). Il reste le chemin par défaut.
- **Noyaux par lot communs.** Gradients, appel de MGIS, résidu, action de
  la tangente et diagonale sont écrits une fois, génériques sur le type de
  gradient, avec `mfem::forall` (portable CUDA, HIP et CPU).
- **Deux points d'entrée** sur ces noyaux :
  - assemblage partiel : les intégrateurs de mfem-mgis implémentent les
    méthodes de `NonlinearFormIntegrator`. `AddMultPA(x, y)` calcule les
    gradients, intègre le lot et ajoute `Bᵀ σ`. `AssembleGradPA` prépare
    l'opérateur tangent (déjà calculé par `AddMultPA` dans le cas
    courant). `AddMultGradPA(v, y)` ajoute `Σ Bᵀ K B v`.
    `AssembleGradDiagonalPA` donne la diagonale, pour un lisseur de
    Jacobi ;
  - opérateur propre à mfem-mgis : un `mfem::Operator` pour le résidu et un
    pour la tangente, sans passer par `NonlinearForm`. Il laisse plus de
    liberté (plusieurs matériaux soumis en asynchrone, ordre des noyaux,
    fusion de noyaux).
- Les classes actuelles (une par hypothèse, symétrie et cinématique)
  gardent leur rôle.
- **Données aux points de Gauss** dans la mémoire de MFEM
  (`mfem::Memory`), point par point ou composante par composante. MGIS les
  voit sans copie, par l'API de vues à concevoir. Le chemin hôte lit la
  même mémoire.
- Rotations des lois orthotropes : fonctions de rotation de lots de MGIS,
  avec des matrices données en tableau.
- Choix à l'exécution : périphérique de MFEM (`mfem::Device`),
  bibliothèque de la loi, point d'entrée, rangement des données, option de
  l'opérateur tangent, solveur. MGIS vérifie que la mémoire de MFEM
  convient à la bibliothèque (emplacement contre espace mémoire).

### Piste pour plus tard : la loi dans le noyau élément

Le chemin prévu sépare le noyau des gradients, le lot MGIS et le noyau du
résidu. Les données passent par la mémoire globale entre les trois. C'est
modulaire, mais demande plus de bande passante. Ratel, libCEED et dFEM
font au contraire tout le travail d'un point dans un seul noyau.

Pour essayer cette variante, TFEL devrait offrir une nouvelle sortie : une
fonction d'un point, générée en en-tête, compilée dans le noyau du code
appelant. Avec CUDA, appeler depuis un noyau une fonction d'une autre
bibliothèque demande du code relogeable. Cette sortie permettrait aussi
d'utiliser des lois MFront dans libCEED ou Ratel. Gros chantier dans TFEL,
à comparer au chemin séparé.

## Portabilité

| Périphérique de MFEM | Backend TFEL de la loi |
|---|---|
| CUDA | `cuda`, ou `kokkos` avec CUDA |
| HIP | `hip` |
| CPU (séquentiel, OpenMP) | `openmp`, `stlpar`, ou bibliothèque `generic` |

MFEM n'a pas de backend SYCL : les backends SYCL de TFEL servent d'autres
codes. Les lois HIP de TFEL s'exécutent ici sur la plateforme NVIDIA ;
mfem-mgis avec un MFEM compilé avec HIP pour un GPU AMD se compile et se
lie (`scripts/build-hip-amd.sh`), sans exécution sur un GPU AMD.

## Solveurs

### Ce qui existe dans mfem-mgis

Fabrique `src/LinearSolverFactory.cxx` :

- MFEM : `CGSolver`, `GMRESSolver`, `BiCGSTABSolver`, `MINRESSolver`,
  `SLISolver` (séquentiel et parallèle) ;
- parallèle : `HyprePCG`, `HypreGMRES`, `HypreFGMRES`, `MUMPSSolver` (si
  MFEM a MUMPS) ;
- séquentiel : `UMFPackSolver` (si SuiteSparse) ;
- préconditionneurs : `HypreBoomerAMG`, `HypreEuclid`, `HypreILU`,
  `HypreParaSails`, `HypreDiagScale`, ou aucun ;
- PETSc : un chemin existe (`Config.cxx`), il plantait lors d'essais
  antérieurs.

### Portabilité sur GPU (à confirmer sur la machine)

| Solveur ou préconditionneur | GPU |
|---|---|
| Krylov de MFEM (CG, GMRES, BiCGSTAB, MINRES, SLI) | oui, si l'opérateur et le préconditionneur sont sur le GPU |
| Krylov de hypre (PCG, GMRES, FGMRES) | oui, avec hypre compilé pour CUDA ou HIP |
| `HypreBoomerAMG` | oui avec hypre GPU, pour une partie des options (grossissement PMIS, interpolations directe et étendues, lisseurs Jacobi, l1-Jacobi, Gauss-Seidel hybride, Chebyshev) ; MFEM choisit des valeurs compatibles sur GPU |
| `HypreDiagScale` | oui |
| `HypreILU` | oui avec hypre 2.29 et plus (sans mémoire unifiée) |
| `HypreEuclid`, `HypreParaSails` | non |
| `MUMPSSolver`, `UMFPackSolver` | non : solveurs directs sur CPU, matrice copiée vers l'hôte |
| PETSc | seulement avec PETSc compilé pour GPU ; chemin à réparer |
| Ginkgo (enveloppes `mfem::Ginkgo`) | CPU (OpenMP) et GPU (CUDA, HIP) avec le même code, choisi par l'exécuteur (Ginkgo lui-même a aussi SYCL, pas les enveloppes de MFEM) ; Krylov, Jacobi, ILU, IC, ISAI dans les enveloppes, pas l'AMG de Ginkgo ; matrices séquentielles ou opérateurs sans matrice ; matrices séquentielles seulement dans les enveloppes de MFEM 4.9 (distribuées en MFEM 4.10) ; demande MFEM compilé avec Ginkgo |

Les options « élasticité » de BoomerAMG (grossissement nodal, vecteurs de
modes rigides) ne marchent pas sur GPU : `SetElasticityOptions` arrête le
programme avec hypre GPU (vérifié dans MFEM 4.9 et 4.10). La stratégie
« Elasticity » de mfem-mgis devra être refusée sur GPU, ou remplacée par
« System » (`SetSystemsOptions`, qui marche sur GPU).

### Ajouts envisagés

La liste complète des solveurs, préconditionneurs et techniques retenus
pour l'étude (CPU et GPU, ou GPU seul), avec leur intérêt, ce qui a été
vérifié et les paquets de l'environnement correspondants, est dans
[solveurs-gpu.md](solveurs-gpu.md). L'ordre dans lequel ils entrent dans
mfem-mgis est donné par le [plan détaillé](plan-mfem-mgis-detaille.md).

Référence principale : Ratel (`https://gitlab.com/micromorph/ratel`,
arXiv 2204.01722). p-multigrille sans matrice avec lisseurs Chebyshev ou
Jacobi et AMG sur le Q1 assemblé, Newton-CG avec des pas résolus à 1e-3,
L-BFGS avec un cycle de multigrille comme hessien inverse initial,
recherche linéaire, données aux points gardées pour la tangente.

## Configurations mixtes

Les choix se font à l'exécution, indépendamment : périphérique de MFEM,
bibliothèque de la loi, point d'entrée, rangement des données, option de
l'opérateur tangent, solveur linéaire.

| Configuration | Données | Copies hôte et GPU par itération |
|---|---|---|
| Tout CPU | hôte | aucune |
| Tout GPU | GPU | aucune, sauf les sorties |
| Loi sur GPU, assemblage et solveur sur CPU (lois lourdes) | partagées | gradients vers le GPU, contraintes et tangentes vers l'hôte |
| Assemblage sur GPU, solveur direct sur CPU | GPU, matrice copiée | matrice et second membre vers l'hôte |
| Une loi non portable dans un calcul GPU | ce matériau en mémoire partagée, loi `generic` | données de ce seul matériau |

Les copies des mélanges sont explicites : emplacement `SHARED` déclaré à
MGIS, ou copies de MFEM. Sur le chemin tout GPU, l'absence de copie cachée
est vérifiée de trois façons : périphérique `debug` de MFEM, GPU sans HMM,
profil `nsys`.

## Performance

- Les données restent sur le GPU d'un pas à l'autre. Pas de mémoire
  managée sur le chemin tout GPU.
- Les noyaux `generic-parallel` n'ont pas encore été mesurés sérieusement.
  Premières mesures : loi J2 par lot contre le noyau J2 écrit à la main des
  essais préalables, rangement point par point contre composante par
  composante, puis des lois plus lourdes (Norton implicite, monocristal).
- Plusieurs matériaux : soumission asynchrone et une seule
  synchronisation.
- MPI : un processus par GPU (pile des polycristaux), MPI compatible GPU si
  disponible.

## Validation

- Chemins par lot sur CPU : mêmes résultats que le chemin actuel sur la
  suite de tests de mfem-mgis, pour les deux points d'entrée et les deux
  rangements.
- GPU contre CPU sur les exemples (poutre en plasticité, grandes
  déformations, thermique, polycristal), avec une tolérance d'arrondi.
- Comparaison des variantes : temps par pas et par phase, CPU multicœur
  contre GPU, les trois options de l'opérateur tangent sur les mêmes cas.

## Étapes

Le découpage détaillé de ces étapes en lots (fichiers touchés, tests,
critères de fin, paramètres proposés) est dans
[plan-mfem-mgis-detaille.md](plan-mfem-mgis-detaille.md).

1. **MGIS : API sur des vues** pour les deux rangements des données aux
   points (lois, modèles, propriétés matériau ; synchrone et asynchrone),
   et mise à jour des états sur GPU. Tests sur CPU et CUDA, relecture
   indépendante. Fait le 2026-10-09 (voir
   [mgis.md](mgis.md#vues-et-rangements)).
2. **MFEM 4.10 avec CUDA** dans l'environnement `mmm-gpu`, avec toutes
   les bibliothèques de solveurs (peu de tâches à la fois : `ptxas` prend
   22 Go pour un fichier).
3. **mfem-mgis sur CPU** : données aux points dans la mémoire de MFEM, puis
   noyaux par lot (gradients, MGIS, résidu, action de la tangente), petites
   déformations 3D, résultats identiques au chemin actuel. Les deux points
   d'entrée et les deux rangements.
4. **Opérateur tangent** : les trois options et le point fixe ; solveurs
   compatibles GPU ; Ginkgo dans la fabrique.
5. **GPU** : périphérique de MFEM, mêmes résultats que sur CPU à l'arrondi
   près, sans copie cachée (périphérique `debug`, GPU sans HMM, `nsys`).
6. **Fin de pas sans copie vers l'hôte**, plusieurs matériaux en
   asynchrone.
7. **Tout mfem-mgis** : hypothèses, cinématiques, rotations orthotropes,
   propriétés matériau, variables externes.
8. **Exemples** de mfem-mgis-examples et de mm-opera-hpc avec les options
   CPU/GPU et les autres options ; mesures et comparaison des variantes.
9. **Documentation** du prototype et préparation des propositions amont
   ([integration-amont.md](integration-amont.md)).

À chaque étape : tests verts (passe ciblée), relecture indépendante avant
de déclarer l'étape finie, journal à jour, commit, synchronisation, push,
bundle.

## Questions ouvertes

- Propriétaire de la mémoire aux points sur GPU : `mfem::Memory` dans les
  `PartialQuadratureFunction`, vues par MGIS. La fin d'un pas de temps
  sans copie vers l'hôte est réglée côté MGIS : copie sur le GPU par
  `update` avec des opérations fournies par l'appelant. L'échange des
  tableaux est écarté : il n'équivaut pas à `update` (une prédiction
  n'écrit pas l'état de fin de pas).
- Proposition amont : un chemin parallèle activé par option, acceptable
  par les mainteneurs de mfem-mgis.

## Risques

- Compilation de MFEM 4.10 avec CUDA 13 : durée et mémoire.
- mfem-mgis n'a été utilisé qu'avec MFEM 4.8 et 4.9 : compatibilité avec
  4.10 à vérifier en premier.
- Ampleur : mfem-mgis a beaucoup d'intégrateurs. Le chemin par lot est
  générique et commence par les cas les plus utilisés.
- Coût de construction de l'AMG sur GPU, mémoire.
