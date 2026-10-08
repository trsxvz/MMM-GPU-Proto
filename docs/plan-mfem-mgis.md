# Plan pour mfem-mgis

Statut : proposition, à discuter avant tout code. mfem-mgis n'est pas
encore modifié (base master `4c8a91d815`, dépôt `mfem-mgis/`). Ce plan
s'appuie sur TFEL et MGIS tels que décrits dans [tfel.md](tfel.md) et
[mgis.md](mgis.md).

## Point de départ

### mfem-mgis aujourd'hui

- mfem-mgis assemble élément par élément sur l'hôte. Il y a un intégrateur
  MFEM par matériau (`MultiMaterialNonLinearIntegrator`,
  `AssembleElementVector`, `AssembleElementGrad`).
- Dans chaque élément, il calcule les gradients aux points de Gauss,
  intègre la loi point par point
  (`BehaviourIntegratorBase::performsLocalBehaviourIntegration`), puis
  ajoute les contributions au résidu et à la matrice tangente.
- Les données aux points de Gauss sont dans des
  `PartialQuadratureFunction` (stockage `std::vector` ou vue externe). Le
  matériau est un `MaterialDataManager` de MGIS.

### Ce que MGIS offre désormais

- Lots de lois, de modèles et de propriétés matériau, sur les tableaux de
  l'appelant, sans copie.
- Basculement automatique des fonctions de MGIS qui travaillent sur un
  `MaterialDataManager`. Le chemin actuel de mfem-mgis intègre point par
  point, sur une vue d'un point : il demande une bibliothèque `generic`.
  Une bibliothèque de lots y est refusée avec un message clair.
- Politique d'exécution (file), emplacement des tableaux vérifié,
  ressources mémoire.
- Calcul synchrone par défaut, ou asynchrone avec une politique
  asynchrone : plusieurs matériaux peuvent être soumis avant une seule
  synchronisation.
- Pas encore : mise à jour des états sur GPU (`update`, `revert`). Elle
  sera ajoutée pendant ce travail (voir [Étapes](#étapes)).

### Essais préalables

Des essais faits hors de ce dépôt (dossier `mfem-mgis-gpu`) ont comparé, sur
une poutre en plasticité et 8 processus :

- Newton avec PCG et AMG : 14,4 s sur CPU ;
- point fixe avec l'opérateur élastique et un cycle AMG : 8,1 s sur CPU ;
- sur GPU, avec une loi J2 écrite à la main en CUDA : 0,69 s, 0,34 s et
  0,15 s pour trois variantes.

Sur un polycristal, le point fixe converge mal : Newton avec GMRES et AMG
prend 22,6 s, le point fixe 67 s. Aucune méthode ne gagne partout.

## Architecture cible

À chaque itération, pour chaque matériau, tout reste sur le GPU :

1. **Gradients aux points de Gauss.** Les valeurs élémentaires sont
   extraites de la solution (`ElementRestriction` de MFEM), dérivées aux
   points (`QuadratureInterpolator`), puis converties en gradients de la loi
   (déformation, gradient de la transformation, gradient de température),
   écrits directement dans les gradients de fin de pas du matériau.
2. **Intégration du lot** par MGIS, sur ces mêmes tableaux.
3. **Résidu** : `Bᵀ σ w det J` par point, somme par élément, puis
   assemblage parallèle par MFEM sur le GPU.
4. **Opérateur tangent**, au choix :
   - a. **sans matrice** : action `Σ Bᵀ K B` aux points de Gauss, solveur de
     Krylov préconditionné par un AMG de l'opérateur élastique, construit
     une fois ;
   - b. **point fixe** avec l'opérateur élastique : seul le résidu est
     calculé à chaque itération, efficace pour les lois peu non linéaires ;
   - c. **tangente assemblée sur le GPU** : matrices élémentaires, puis
     matrice globale hypre. Newton exact, plus coûteux.
5. **Solveurs linéaires** : hypre sur GPU (PCG, GMRES, BoomerAMG), Ginkgo
   en option (CUDA, HIP, OpenMP).
6. **Fin de pas** : `update` et `revert` de MGIS sur le GPU. mfem-mgis
   possède la mémoire, MGIS sait quels blocs copier. Les transferts vers
   l'hôte se limitent aux sorties.
7. **Post-traitements** : sur l'hôte d'abord, sur le GPU ensuite.

Les trois options de l'opérateur tangent sont implémentées, choisies par
une option du problème, et comparées sur les mêmes cas, sur CPU et sur
GPU.

## Organisation du code

Le chemin GPU est l'**assemblage partiel de MFEM**
(`AssemblyLevel::PARTIAL` de `NonlinearForm`), pas un opérateur écrit à
côté de MFEM.

- **Le chemin actuel ne change pas** (`AssemblyLevel::LEGACY`, élément par
  élément sur l'hôte). Il reste le chemin par défaut.
- Les intégrateurs de mfem-mgis implémentent en plus les méthodes
  d'assemblage partiel de `NonlinearFormIntegrator`, matériau par
  matériau :
  - `AddMultPA(x, y)` : gradients, intégration du lot par MGIS, `y += Bᵀ σ` ;
  - `AssembleGradPA(x, fes)` : opérateur tangent aux points pour l'état
    `x`, déjà calculé par `AddMultPA` dans le cas courant ;
  - `AddMultGradPA(v, y)` : `y += Σ Bᵀ K B v` (option a) ;
  - `AssembleGradDiagonalPA` : diagonale, pour un préconditionneur de
    Jacobi.
- Les classes actuelles (une par hypothèse, symétrie et cinématique)
  gardent leur rôle. Leur code par lot est écrit une fois, générique sur le
  type de gradient, avec `mfem::forall` (portable CUDA, HIP et CPU).
- Les `PartialQuadratureFunction` sont stockées dans la mémoire de MFEM
  (`mfem::Memory`). MGIS les voit par les tableaux externes de ses
  initialiseurs, sans copie. Le chemin hôte lit la même mémoire.
- Rotations des lois orthotropes : fonctions de rotation de lots de MGIS,
  avec des matrices données en tableau.
- Propriétés matériau et variables externes : uniformes, tableaux sur le
  GPU calculés par MFEM, ou propriétés matériau `generic-parallel`
  évaluées par MGIS.
- Choix à l'exécution : périphérique de MFEM (`mfem::Device`) et
  bibliothèque de la loi. MGIS vérifie que la mémoire de MFEM convient à la
  bibliothèque (emplacement contre espace mémoire).

## Portabilité

| Périphérique de MFEM | Backend TFEL de la loi |
|---|---|
| CUDA | `cuda`, ou `kokkos` avec CUDA |
| HIP | `hip` |
| CPU (séquentiel, OpenMP) | `openmp`, `stlpar`, ou bibliothèque `generic` |

MFEM n'a pas de backend SYCL : les backends SYCL de TFEL servent d'autres
codes. HIP se compile ici mais ne s'exécute que sur la plateforme NVIDIA.

## Configurations mixtes

Les choix se font à l'exécution, indépendamment : périphérique de MFEM,
bibliothèque de la loi, niveau d'assemblage, option de l'opérateur tangent,
solveur linéaire.

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
- **Les noyaux `generic-parallel` n'ont pas encore été mesurés
  sérieusement.** Première mesure avant le travail GPU de mfem-mgis : loi
  J2 par lot contre le noyau J2 écrit à la main des essais préalables, puis
  des lois plus lourdes (Norton implicite, monocristal).
- Plusieurs matériaux : soumission asynchrone et une seule
  synchronisation.
- MPI : un processus par GPU (pile des polycristaux), MPI compatible GPU si
  disponible.

## Validation

- Chemin par lot sur CPU : mêmes résultats que le chemin actuel sur la
  suite de tests de mfem-mgis.
- GPU contre CPU sur les exemples (poutre en plasticité, grandes
  déformations, thermique, polycristal), avec une tolérance d'arrondi.
- Campagne de performance : temps par pas et par phase, CPU multicœur
  contre GPU, les trois options de l'opérateur tangent sur les mêmes cas.

## Étapes

0. Compiler MFEM 4.10 avec CUDA dans l'environnement `mmm-gpu` (une nuit,
   peu de tâches : un fichier demande jusqu'à 22 Go). Mesurer les noyaux
   `generic-parallel`.
1. Stockage dans la mémoire de MFEM, CPU seulement, tests inchangés.
2. Assemblage partiel sur CPU (gradients de tous les points, MGIS, résidu,
   action de la tangente), petites déformations 3D d'abord. Résultats
   identiques au chemin actuel.
3. Noyaux GPU des gradients et du résidu. Résidu GPU égal au résidu CPU.
4. Opérateur tangent, options a et b. hypre sur GPU.
5. Option c : tangente assemblée sur le GPU.
6. Autres hypothèses et types de gradients, rotations orthotropes.
7. Mise à jour des états sur GPU dans MGIS, plusieurs matériaux soumis en
   asynchrone, sorties.
8. Ginkgo, compilation HIP.
9. Campagne de performance, documentation.

## Questions à trancher

1. **Point d'entrée dans MFEM** : assemblage partiel de `NonlinearForm`
   (proposé), ou opérateur propre à mfem-mgis ?
2. **Options de l'opérateur tangent** : les trois sont prévues et comparées.
   Faut-il en prioriser une pour le premier cas complet ?
3. **Mémoire des données aux points de Gauss** : `mfem::Memory` dans les
   `PartialQuadratureFunction` convient-il à tous les usages actuels
   (sorties, post-traitements, reprises) ?
4. **Premier cas cible** : poutre en plasticité en petites déformations 3D
   (proposé), ou un cas plus proche des études (polycristal, grandes
   déformations) ?
5. **Calcul asynchrone dans mfem-mgis** : MGIS l'offre. Faut-il s'en servir
   dès les premiers matériaux multiples, ou à l'étape 7 ?
6. **Propriétés matériau aux points de Gauss** : fonctions de MFEM, ou
   propriétés matériau `generic-parallel` évaluées par MGIS ?
7. **Proposition amont** : un chemin parallèle activé par option, pour
   rester acceptable par les mainteneurs de mfem-mgis.

## Risques

- Compilation de MFEM avec CUDA 13 : durée et compatibilité.
- Ampleur : mfem-mgis a beaucoup d'intégrateurs. Le chemin par lot est
  générique et commence par les cas les plus utilisés.
- Coût de construction de l'AMG sur GPU, mémoire.
