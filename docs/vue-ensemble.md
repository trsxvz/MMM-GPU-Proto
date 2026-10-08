# Vue d'ensemble

Ce document présente le prototype : son but, le rôle de chaque
bibliothèque, les principes suivis et l'état d'avancement. Les détails
techniques sont dans [tfel.md](tfel.md) et [mgis.md](mgis.md).

## But

Une simulation mécanique avec mfem-mgis passe beaucoup de temps dans les
lois de comportement. Elles sont évaluées en chaque point d'intégration, à
chaque itération. Ce calcul est très parallèle : les points sont
indépendants.

Le prototype rend toute la chaîne MFront, MGIS et mfem-mgis capable de
calculer sur GPU comme sur CPU. Le même fichier MFront doit servir partout.
Le choix entre CPU et GPU se fait à la compilation de la loi ou à
l'exécution, pas dans le code de la loi.

## La chaîne et le rôle de chaque bibliothèque

Chaque bibliothèque ne connaît que l'interface de la précédente.

| | TFEL/MFront | MGIS | mfem-mgis |
|---|---|---|---|
| Rôle | code d'un point matériel, traitement d'un lot de points | données matériau, chargement des bibliothèques, aiguillage | éléments finis, solveurs, orchestration |
| Code GPU | noyaux des fonctions de lot (un backend par bibliothèque générée) | aucun dans le cœur | assemblage et solveurs (MFEM, hypre, Ginkgo) |
| Mémoire | n'alloue rien (sauf les copies SYCL de valeurs inaccessibles) | utilise les tableaux de l'appelant, ou une ressource mémoire donnée par l'appelant | possède les données, décide où elles vivent |
| Interface exposée | structures C versionnées, symboles de description | `MaterialDataManager`, fonctions de lot | formes non linéaires de MFEM |

Un lot est un ensemble de points d'intégration d'un même matériau. Une
fonction de lot traite tous ces points en un appel.

## Principes

1. **Même sens que l'interface `generic`.** Un appel sur un lot équivaut à
   l'appel de la fonction `generic` sur chaque point. Les résultats de
   l'interface `generic` sur CPU servent de référence à tous les tests.
2. **Pas de physique dupliquée.** La classe générée par MFront pour une loi
   sert telle quelle sur CPU et sur GPU. L'interface `generic-parallel`
   n'ajoute qu'un constructeur et un peu de code d'adaptation.
3. **Logique dans des en-têtes, code généré mince.** Le travail en un point
   est écrit une fois, en gabarits utilisables sur hôte et sur GPU. Chaque
   backend n'écrit que sa boucle parallèle.
4. **Portabilité par backend.** Une bibliothèque générée vise un backend :
   CUDA, HIP, algorithmes parallèles du C++, Kokkos, SYCL ou OpenMP.
   Plusieurs bibliothèques d'une même loi peuvent être chargées ensemble.
   Le choix entre CPU et GPU se fait alors à l'exécution.
5. **L'appelant choisit la mémoire.** TFEL ne choisit ni le GPU ni la
   mémoire. Chaque bibliothèque déclare seulement où les données doivent
   être accessibles. MGIS vérifie que les tableaux de l'appelant
   conviennent.
6. **Erreurs sans coût en cas de succès.** Les noyaux n'écrivent rien
   tant qu'aucun point n'échoue. Une erreur d'argument ou du backend lève
   une exception. Un échec du calcul en un point est rapporté comme tel,
   jamais confondu avec une erreur.
7. **Symétrie.** Lois, modèles et propriétés matériau ont la même
   interface d'exécution : file d'exécution, calcul synchrone ou
   asynchrone, espace mémoire, version du contrat. Les différences
   restantes sont justifiées dans [tfel.md](tfel.md#comparaison-entre-lois-et-propriétés-matériau).
8. **Aucun changement des résultats CPU des interfaces existantes**, sauf
   les corrections de défauts. Chaque correction est isolée et justifiée
   dans [tfel.md](tfel.md#corrections-du-code-existant).

## Vocabulaire

- **Lot** : ensemble de points d'intégration traités par un appel.
- **Backend** : modèle de programmation visé par une bibliothèque générée
  (`cuda`, `hip`, `stlpar`, `kokkos`, `sycl`, `openmp`).
- **Champ** : valeurs d'une variable sur tous les points d'un lot, données
  par un pointeur et des pas. Un pas nul entre points désigne une valeur
  uniforme.
- **Espace mémoire** d'une bibliothèque : où ses calculs lisent les
  données (`host`, `device` ou `any`). Il est fixé à la compilation.
- **Emplacement** des tableaux : où l'appelant les a placés (hôte, GPU ou
  mémoire partagée). MGIS le compare à l'espace mémoire.
- **File d'exécution** : flux CUDA ou HIP, file SYCL, instance Kokkos.
- **Synchrone** : la fonction rend la main une fois le calcul fini.
  **Asynchrone** : elle rend la main dès le calcul soumis à la file.
  L'appelant synchronise la file avant de lire les résultats.

## État d'avancement

| Partie | État | Détails |
|---|---|---|
| TFEL/MFront : interface `generic-parallel` des lois et des modèles | fait, validé sur 14 configurations | [tfel.md](tfel.md) |
| TFEL/MFront : interface `generic-parallel` des propriétés matériau | étendue (options d'exécution, asynchrone, nouveaux backends) | [tfel.md](tfel.md#interface-generic-parallel-des-propriétés-matériau) |
| TFEL/MFront : code GPU et corrections | fait | [tfel.md](tfel.md#code-gpu-dans-tfel-et-mfront) |
| MGIS : lots de lois, de modèles et de propriétés matériau | fait | [mgis.md](mgis.md) |
| MGIS : calcul asynchrone, mise à jour des états sur GPU | à faire, après les premiers essais de mfem-mgis | [mgis.md](mgis.md#limites-actuelles) |
| mfem-mgis | non modifié, plan à discuter | [plan-mfem-mgis.md](plan-mfem-mgis.md) |
| Performance des noyaux | non mesurée sérieusement | [plan-mfem-mgis.md](plan-mfem-mgis.md#performance) |

La validation est décrite dans [validation.md](validation.md). Le
découpage en contributions pour les dépôts officiels est dans
[integration-amont.md](integration-amont.md). Le cheminement, jour par
jour, est dans le [journal](../JOURNAL.md).
