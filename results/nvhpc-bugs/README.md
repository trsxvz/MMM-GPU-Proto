# Défauts de compilateurs GPU rencontrés (reproducteurs minimaux)

| Fichier | Compilateur | Défaut | Contournement dans le proto |
|---|---|---|---|
| `signed-loop-trip-count.cu` | nvc++ 25.11 (`-cuda`, `-stdpar=gpu`) | boucle à compteur signé et test `!=` : N+1 itérations sur GPU à -O2/-O3 | boucles non signées dans le code concerné |
| `empty-kernel-stdpar.cxx` | nvc++ 25.11 | mode mémoire unifiée (`-gpu=mem:unified`, défaut) : écriture à une adresse de pile | `-gpu=mem:managed` |
| `aggregate-zero-init.cxx` | nvc++ 25.11 `-stdpar=gpu` | initialiseur par défaut `= {}` d'un membre tableau : constante LLVM invalide (« parse integer constant must have integer type ») | membre initialisé dans le constructeur (`PointView`) |
| `member-function-pointer.cxx` | nvc++ 25.11 `-stdpar=gpu` | appel par pointeur sur fonction membre ignoré **silencieusement** sur GPU (-O0 à -O3) | objets fonctions nommés générés par MFront |
| (dans le journal) | nvc++ 25.11 `-stdpar=gpu` | `std::variant` (libstdc++) : « Unexpected branch type », alias LLVM invalides | `TrivialGenTypeBase` pour l'opérateur tangent |
| `nvcc-static-constexpr-object.cu` | nvcc 13.0 | copie device d'un objet `static constexpr` local remplie de **zéros, sans diagnostic** si un constructeur est appelé | variables `__device__ constexpr` d'espace de noms (`tfel::utilities::getConstantObject`) |
| `nvcc-decltype-functional-cast.cu` | nvcc 13.0 (cudafe++) | `auto i = decltype(s){};` avec `const auto s = c.size();` dans un patron utilisé en code device : réécrit pour le compilateur hôte en une conversion fonctionnelle invalide (g++ : « expected primary-expression before ‘auto’ ») | `size_type{}` (`LinearInterpolation.ixx`) |
| (dans le journal, 05/10 soir, C5) | nvc++ 25.11 `-mp=gpu` (OpenMP, -O2 et -O3) | lois des DSL implicites : `computeConsistentTangentOperator` renvoie faux sur GPU (décomposition LU du jacobien) alors que le Newton a convergé ; juste à -O0 et -O1 ; un `printf` dans la branche d'échec corrige le résultat ; non reproduit hors contexte | configuration `openmp-offload-nvhpc-gpu` à -O1 |
| (dans le journal, 05/10 soir, C5) | nvc++ 25.11 `-mp=gpu` (OpenMP) | le code device de `HosfordPerfectPlasticity` déclare un bloc statique interne (`_st_<hash>__T1_0`, 317 octets, jamais utilisé) que l'optimiseur supprime ; le support d'exécution le cherche au chargement (« Could not find symbol ... on the device ») et aucune loi de la bibliothèque ne peut être exécutée | configuration `openmp-offload-nvhpc-gpu` non activée |
| `nvcc-cicc-initialize-function-reproducer.cu` (non minimal : inclut le début de la source générée de `SaintVenantKirchhoffElasticity` pour l'interface `generic-parallel`, compiler avec `-DUSE_BACKEND -I include` dans un dossier de génération de la loi) | nvcc 13.0, optimiseur NVVM (`cicc`) | la contrainte écrite par une fonction d'initialisation de loi en grandes déformations est perdue quand l'état de fin de pas est copié après la construction de la vue du point ; juste avec `-G`, `-Xcicc -O0` et clang CUDA ; le résultat change avec des `printf` | copie de l'état de fin de pas dans le constructeur de `PointView` |

Contraintes (non des défauts) : nvc++ `-stdpar=gpu` (mémoire gérée)
interdit l'accès aux variables globales et statiques depuis le code GPU
(erreur explicite) : `getConstantObject` y renvoie une copie.

nvcc accepte une définition sans `__host__ __device__` quand la
déclaration en porte (et inversement), clang CUDA la refuse. Ce n'est pas
un défaut de nvcc mais une incohérence de TFEL, corrigée dans 85
définitions (balayage clang de tous les en-têtes, voir le journal).
