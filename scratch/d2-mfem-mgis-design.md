# D2, mfem-mgis : calcul sur GPU (plan)

Statut : plan du 2026-10-08, avant tout code. Bases : mfem-mgis master
`4c8a91d815` (`mfem-mgis/`), MGIS selon `scratch/d1-mgis-design.md` (étapes 1 à 4a
faites le 2026-10-09 : chargement, intégration sans copie, fonctions
d'initialisation, post-traitements, rotations, emplacement des données), TFEL du
proto (R1 à R8). Essais préalables : point fixe à opérateur constant et
solveur GPU écrit à la main (loi J2 en CUDA), dossier autonome
`mfem-mgis-gpu` hors de ce dépôt.

## 1. Point de départ

- mfem-mgis assemble élément par élément sur l'hôte : un intégrateur MFEM
  par matériau (`MultiMaterialNonLinearIntegrator`, `AssembleElementVector`,
  `AssembleElementGrad`) ; dans chaque élément, gradients au point de
  Gauss, intégration de la loi point par point
  (`BehaviourIntegratorBase::performsLocalBehaviourIntegration`), puis
  contributions au résidu et à la matrice tangente.
- Les données aux points de Gauss sont dans des
  `PartialQuadratureFunction` (stockage `std::vector` ou vue externe) ; le
  matériau est un `MaterialDataManager` de MGIS.
- Dépendances de l'environnement mmm-gpu : hypre et Ginkgo 1.11 compilés
  avec CUDA. MFEM 4.10 avec CUDA reste à compiler (plusieurs heures, un
  fichier demande jusqu'à 22 Go : compilation seule, peu de tâches).
- Essais préalables (poutre, 8 processus) : Newton PCG+AMG 14,4 s, point
  fixe avec l'opérateur élastique et un cycle AMG 8,1 s sur CPU ; sur GPU,
  0,69 s, 0,34 s et 0,15 s pour les trois variantes. Polycristal : le
  point fixe converge mal (Newton GMRES+AMG 22,6 s contre 67 s).

## 2. Architecture cible

À chaque itération, pour chaque matériau, tout reste sur le GPU :

1. **Gradients aux points de Gauss** : valeurs élémentaires extraites de
   la solution (`ElementRestriction` de MFEM), dérivées aux points
   (`QuadratureInterpolator`), puis gradient de la loi (déformation,
   gradient de la transformation, gradient de température) écrit
   directement dans les gradients de fin de pas du matériau.
2. **Intégration du lot** : MGIS (D1) appelle la bibliothèque
   `generic-parallel` sur ces mêmes tableaux.
3. **Résidu** : `Bᵀ σ w det J` par point, somme par élément, puis
   assemblage parallèle (opérations de MFEM sur le GPU).
4. **Opérateur tangent** (méthodes d'assemblage partiel de MFEM, section
   3), au choix :
   - a. **sans matrice** : action `Σ Bᵀ K B` aux points de Gauss, solveur
     de Krylov préconditionné par un AMG de l'opérateur élastique,
     assemblé une seule fois ;
   - b. **point fixe** avec l'opérateur élastique (résidu seul à chaque
     itération, efficace pour les lois peu non linéaires) ;
   - c. **tangente assemblée sur le GPU** (matrices élémentaires, puis
     matrice globale hypre) : Newton exact, plus coûteux.
5. **Solveurs linéaires** : hypre sur GPU (PCG, GMRES, BoomerAMG) ;
   Ginkgo en option (CUDA, HIP, OpenMP).
6. **Fin de pas** : `update`/`revert` de MGIS avec la copie de MFEM
   (mfem-mgis possède la mémoire, MGIS connaît les blocs à copier) ;
   transferts vers l'hôte seulement pour les sorties.
7. **Post-traitements** : sur l'hôte d'abord ; plus tard, sur le GPU avec
   le module `mgis::function` de MGIS (évaluateurs, algorithmes CUDA).

Les essais préalables montrent qu'aucune option ne gagne partout : b gagne
sur les lois peu non linéaires, a ou c sur les lois très non linéaires. Les
trois options sont donc implémentées, choisies par une option du problème,
et comparées sur les mêmes cas (CPU et GPU).

## 3. Organisation du code

Le chemin GPU est l'**assemblage partiel de MFEM** (`AssemblyLevel::PARTIAL`
de `NonlinearForm`), pas un opérateur à côté de MFEM :

- **Le chemin actuel (`AssemblyLevel::LEGACY`, élément par élément sur
  l'hôte) ne change pas** et reste le chemin par défaut.
- Les intégrateurs de mfem-mgis implémentent en plus les méthodes
  d'assemblage partiel de `NonlinearFormIntegrator`, matériau par matériau
  (`MultiMaterialNonLinearIntegrator` boucle sur les matériaux comme
  aujourd'hui) :
  - `AddMultPA(x, y)` : gradients aux points de Gauss, intégration du lot
    par MGIS (avec l'opérateur tangent si le problème le demande, comme
    `AssembleElementVector` aujourd'hui), `y += Bᵀ σ` ;
  - `AssembleGradPA(x, fes)` : opérateur tangent aux points pour l'état
    `x` (déjà calculé par `AddMultPA` dans le cas courant) ;
  - `AddMultGradPA(v, y)` : `y += Σ Bᵀ K B v` (option a) ;
  - `AssembleGradDiagonalPA` : diagonale, pour un préconditionneur de
    Jacobi.
- Le problème non linéaire choisit le niveau d'assemblage par une option ;
  le point fixe (option b) et l'opérateur élastique de l'AMG passent par
  les fabriques de solveurs existantes ; l'option c remplit une matrice
  globale à partir des matrices élémentaires calculées sur le GPU.
- Les classes actuelles (une par hypothèse, symétrie et cinématique)
  gardent leur rôle (choix de la loi, rotations, options) ; leur code par
  lot est écrit une fois, générique sur le type de gradient (petites
  déformations, grandes déformations, thermique), avec `mfem::forall`
  (portable CUDA, HIP et CPU).
- Stockage des `PartialQuadratureFunction` dans la mémoire de MFEM
  (`mfem::Memory`), visible par MGIS par les tableaux externes de ses
  initialiseurs (aucune copie) ; le chemin hôte lit la même mémoire.
- Autres éléments :
  - rotations des lois orthotropes par lot (fonctions de rotation de TFEL,
    par MGIS) ;
  - propriétés matériau et variables externes : uniformes ou tableaux sur
    le GPU (fonctions de MFEM aux points de Gauss, ou propriétés matériau
    `generic-parallel` par MGIS, étape 6 de D1).
- Choix à l'exécution : périphérique de MFEM (`mfem::Device`) et
  bibliothèque de la loi ; mfem-mgis vérifie que le backend de la loi
  (donné par la description de la loi dans MGIS, jamais lu directement
  dans la bibliothèque) correspond au périphérique.

## 4. Portabilité

| MFEM | Backend TFEL de la loi |
|---|---|
| CUDA | `cuda` (ou `kokkos` avec CUDA) |
| HIP | `hip` |
| CPU (séquentiel, OpenMP) | `openmp`, `stlpar` |

MFEM n'a pas de backend SYCL : les backends SYCL de TFEL servent d'autres
codes. HIP se compile ici mais ne s'exécute que sur la plateforme NVIDIA
(pas de GPU AMD).

## 5. Configurations : tout CPU, tout GPU, mixtes

Les choix se font à l'exécution, indépendamment : périphérique de MFEM,
bibliothèque de la loi (`generic` par point ou `generic-parallel` et son
backend), niveau d'assemblage, option de l'opérateur tangent, solveur
linéaire. Seule contrainte vérifiée : le backend de la loi doit pouvoir
lire la mémoire du périphérique de MFEM.

| Configuration | Données | Copies hôte/GPU par itération |
|---|---|---|
| Tout CPU (chemin actuel ou assemblage partiel, loi `generic`, `openmp` ou `stlpar`) | hôte | aucune |
| Tout GPU (MFEM CUDA ou HIP, loi `cuda`, `hip` ou `kokkos`, hypre ou Ginkgo sur GPU) | GPU | aucune (sorties seulement) |
| Loi sur GPU, assemblage et solveur sur CPU (lois lourdes : polycristaux) | partagées (`SHARED`) | gradients vers le GPU, contraintes et tangentes vers l'hôte : inhérentes au mélange |
| Assemblage sur GPU, solveur direct sur CPU (MUMPS) | GPU, matrice copiée | matrice et second membre vers l'hôte |
| Matériaux mélangés : une loi non portable sur GPU dans un calcul GPU | ce matériau en `SHARED`, loi `generic` | données de ce seul matériau |

Les copies des mélanges sont le prix du mélange, pas des copies inutiles :
elles sont explicites (emplacement `SHARED` déclaré à MGIS, ou copies de
MFEM), et mesurées par la campagne de performance.

**Vérification de l'absence de copies cachées** sur le chemin tout GPU :

- périphérique `debug` de MFEM, qui protège la mémoire du GPU simulé et
  arrête le calcul à tout accès de l'hôte ;
- GPU sans HMM (comme pour la validation de TFEL du 08/10) ;
- profil `nsys` : nombre et volume des transferts par itération (seulement
  les sorties attendues).

## 6. Performance

- Données résidentes sur le GPU d'un pas à l'autre ; aucune migration de
  pages (pas de mémoire managée dans le chemin tout GPU).
- Intégration des lois : **les noyaux `generic-parallel` n'ont pas encore
  été mesurés**. Première mesure avant le travail GPU de mfem-mgis : loi J2
  par lot contre le noyau J2 écrit à la main des essais préalables, et
  quelques lois plus lourdes (Norton implicite, monocristal).
- Plusieurs matériaux : soumission asynchrone, une synchronisation (D1,
  étape 5).
- Lois simples : le solveur domine ; options a et b, AMG de l'opérateur
  élastique construit une fois.
- MPI : un processus par GPU (pile des polycristaux), MPI compatible GPU
  si disponible.

## 7. Validation

- Chemin par lot sur CPU : mêmes résultats que le chemin actuel sur la
  suite de tests de mfem-mgis.
- GPU contre CPU sur les exemples (poutre en plasticité, grandes
  déformations, thermique, polycristal), avec une tolérance d'arrondi.
- Campagne de performance : temps par pas et par phase, CPU multicœur
  contre GPU, les trois options de l'opérateur tangent sur les mêmes cas
  (lois peu et très non linéaires).

## 8. Étapes

0. Compiler MFEM 4.10 avec CUDA dans mmm-gpu (une nuit, peu de tâches).
   Mesurer les noyaux `generic-parallel` (section 5).
1. Stockage dans la mémoire de MFEM, CPU seulement, tests inchangés.
2. Assemblage partiel sur CPU (gradients de tous les points, MGIS,
   résidu, action de la tangente), petites déformations 3D d'abord ;
   résultats identiques au chemin actuel.
3. Noyaux GPU (gradients, résidu) ; résidu GPU égal au résidu CPU.
4. Opérateur tangent, options a et b ; hypre sur GPU.
5. Option c : tangente assemblée sur le GPU (matrices élémentaires à
   partir des tangentes locales, matrice globale hypre, AMG reconstruit).
6. Autres hypothèses et types de gradients, rotations orthotropes.
7. Plusieurs matériaux en asynchrone, `update`/`revert` sur le GPU,
   sorties.
8. Ginkgo, compilation HIP.
9. Campagne de performance, documentation.

## 9. Risques

- Compilation de MFEM avec CUDA 13 (durée, compatibilité).
- Ampleur : nombreux intégrateurs dans mfem-mgis ; le chemin par lot est
  générique et commence par les cas les plus utilisés.
- Coût de construction de l'AMG sur GPU, mémoire.
- Proposition amont : chemin parallèle à l'existant, activé par option,
  pour rester acceptable.
