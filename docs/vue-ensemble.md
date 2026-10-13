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
| Code GPU | noyaux des fonctions de lot (un backend par bibliothèque générée) | aucun dans la bibliothèque compilée ; copies et fin de pas sur GPU par des noyaux en en-têtes (`FieldOperations` : CUDA, HIP, Kokkos, SYCL, OpenMP déporté), compilés par l'appelant | noyaux de mfem-mgis (gradients, forces intérieures, tangente, prédiction, post-traitements), assemblage et solveurs de MFEM, hypre, AmgX et cuDSS |
| Mémoire | n'alloue pas les données des points ; alloue une fois par périphérique les copies des données partagées et les blocs d'état et de statistiques, réutilisés d'un appel à l'autre (et, avec SYCL, les copies des valeurs inaccessibles) | utilise les tableaux de l'appelant, ou une ressource mémoire donnée par l'appelant | possède les données, décide où elles vivent |
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
| MGIS : calcul asynchrone | fait, synchrone par défaut | [mgis.md](mgis.md#calcul-asynchrone) |
| MGIS : vues, rangement composante par composante, fin de pas sur GPU | fait | [mgis.md](mgis.md#vues-et-rangements) |
| Documentation de TFEL et de MGIS (pages web, Doxygen, commentaires) | à jour avec les changements | [tfel.md](tfel.md#documentation-dans-larbre-de-tfel), [mgis.md](mgis.md#documentation-dans-larbre-de-mgis) |
| mfem-mgis | lots L0 à L10 faits : données des matériaux dans la mémoire de MFEM, intégration par lots, mécanique par noyaux, point fixe et NewtonK0, thermique, micromorphe, FBar, Faltus, conditions aux limites et autres intégrateurs par noyaux, modèles point par point sur le GPU, post-traitements par noyaux, toute la suite de tests sur le GPU ; Kokkos et AmgX facultatifs ; mémoire de l'hôte standard sur le GPU. Depuis : solveur direct cuDSS et réutilisation de l'analyse de MUMPS, assemblage ordonné de la tangente (`MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY=Ordered`), réglages d'exécution (`--execution cpu|gpu|mixed`, placement de chaque phase, de la loi de chaque matériau par `BehaviourExecutionLocation`), noyaux compilés pour HIP (GPU AMD : compilation vérifiée, sans essai). Tous les exemples de mfem-mgis-examples et les cas de mm-opera-hpc tournent par lots et sur le GPU | [plan-mfem-mgis-detaille.md](plan-mfem-mgis-detaille.md), [plan-execution.md](plan-execution.md), [solveurs-gpu.md](solveurs-gpu.md), [validation.md](validation.md) |
| Performance | poutre r4 sur une RTX 5080 : 1,96 s avec NewtonK0 contre 9,5 s sur 8 cœurs avec la même méthode (13,6 s avec Newton) ; polycristal 9,6 s sur le GPU avec NewtonK0 (8,8 s le 12/10) contre 23,5 s sur 8 cœurs avec la même méthode (29,6 s avec Newton) ; ex8 (5 pas) 2,6 s sur le GPU (solveur direct cuDSS pour les deux problèmes, analyse réutilisée) contre 19,2 s sur 2 cœurs par noyaux (MUMPS symétrique, analyse réutilisée) ; ex7 (40 pas) 2,0 s sur le GPU (cuDSS) contre 5,4 s sur 2 cœurs (MUMPS) ; bulle de mm-opera-hpc 5,3 s sur le GPU contre 25,2 s sur un cœur. Sur le GPU, les petits systèmes linéaires (moins de 10 000 inconnues) sont mieux traités par cuDSS que par GMRES et BoomerAMG | `results/`, [journal](../JOURNAL.md) |

La validation est décrite dans [validation.md](validation.md). Le
découpage en contributions pour les dépôts officiels est dans
[integration-amont.md](integration-amont.md). Le cheminement, jour par
jour, est dans le [journal](../JOURNAL.md).
