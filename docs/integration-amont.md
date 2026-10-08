# Intégration amont

Ce document propose un découpage des changements en contributions pour les
dépôts officiels de TFEL et de MGIS. Chaque contribution doit pouvoir être
relue et testée seule.

Certains commits du prototype mêlent plusieurs sujets. Les contributions
seront donc reconstruites à partir de l'état final des branches, fichier par
fichier, plutôt que par sélection de commits.

## TFEL

L'ordre suit les dépendances : chaque contribution s'appuie sur les
précédentes. Les détails techniques sont dans [tfel.md](tfel.md).

1. **Corrections du code existant**, sans code GPU. Ce sont les corrections
   listées dans [tfel.md](tfel.md#corrections-du-code-existant). Elles
   peuvent changer les résultats sur CPU, puisque ce sont des corrections.
   Chacune est couverte par un test existant ou nouveau.
2. **Code GPU dans TFEL** (en-têtes mathématiques et matériau) :
   annotations hôte et GPU, `TFEL_DEVICE_CODE` pour les compilateurs à une
   passe, lois sans fonctions virtuelles, `TrivialGenTypeBase`,
   `getConstantObject`, valeurs propres sans récursion, contournements de
   compilateurs (par exemple les itérateurs des splines cubiques). Les
   résultats sur CPU ne changent pas. Les incompatibilités sont listées dans
   les notes de version.
3. **Code généré compatible GPU dans MFront** : annotations des classes
   générées, classe de loi compatible GPU, pilotes de l'interface
   `generic` coupés en un cœur hôte et GPU et les fonctions historiques,
   fonctions de rotation et appels de méthodes utilisables sur GPU, type de
   l'opérateur tangent vérifié avant sa conversion, paramètres transmis à la
   classe de loi et à ses variables de comportement, erreurs portables
   (`tfel::fail`, `abortResolution`, `abortIntegration`), données partagées
   (`@SharedData`). Les résultats des interfaces existantes ne changent pas.
4. **Interface `generic-parallel` des lois et des modèles**, avec un
   premier ensemble de backends : CUDA, HIP, algorithmes parallèles du C++
   et OpenMP sur l'hôte. Contrat de données versionné, options d'exécution,
   calcul synchrone ou asynchrone, erreurs sans coût en cas de succès,
   symbole d'espace mémoire, données partagées copiées sur chaque GPU, tests
   et documentation.
5. **Backends et compilateurs supplémentaires** : Kokkos, SYCL, déport de
   `nvc++` (`stdpar` et OpenMP), clang CUDA et HIP. La plupart des
   contournements de compilateurs vont dans cette contribution.
6. **Extensions de l'interface des propriétés matériau** : backends Kokkos
   et OpenMP, propriétés définies par des données ou appelant d'autres lois
   sur CUDA et HIP, options d'exécution en dernier argument avec le même sens
   pour tous les backends (version 2 des fonctions), évaluation asynchrone
   et fonction de compte rendu, symboles de backend, de version et d'espace
   mémoire, code des noyaux partagé par les backends, statuts dans des blocs
   réutilisés, statut 1 avec la politique `Warning`.

Le contrat de données des lois
(`MFront/GenericParallelBehaviour/BatchData.h`) est versionné dès la
quatrième contribution, pour que MGIS puisse le vérifier.

### Vérifications pour chaque contribution

- la suite de tests complète de TFEL (`ctest`) ;
- les comparaisons des tests `generic-parallel` pour les configurations de
  la contribution, à `-O3` ;
- aucun changement des résultats des interfaces existantes, sauf les
  corrections de la première contribution ;
- formatage par `clang-format` limité aux lignes modifiées.

## MGIS

Les contributions de MGIS suivent la quatrième contribution de TFEL, dont
elles lisent le contrat. Les détails sont dans [mgis.md](mgis.md).

1. **Corrections** : axes des matrices de rotation, propriété matériau
   donnée par un `std::span` d'une valeur, statut d'échec des fonctions
   d'initialisation et des post-traitements point par point.
2. **Lots des lois et des modèles** : chargement et vérifications, copies
   des structures de TFEL et test de leur disposition, intégration sans
   copie, fonctions d'initialisation, post-traitements, rotations,
   basculement automatique des fonctions existantes, politique d'exécution,
   refus des fonctions d'un point.
3. **Mémoire** : emplacement des tableaux vérifié contre l'espace mémoire
   des bibliothèques, refus des accès de l'hôte aux tableaux placés sur un
   GPU, ressources mémoire et leurs en-têtes.
4. **Propriétés matériau** : chargement des propriétés `generic-parallel`
   et fonction `evaluate`.
5. **Documentation** : page web des lots et notes de version.

Le test CUDA de MGIS demande un GPU NVIDIA. Il est désactivé par défaut
(option `enable-generic-parallel-cuda-tests`), pour que l'intégration
continue de MGIS reste sur CPU.
