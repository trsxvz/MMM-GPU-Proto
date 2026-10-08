# Défauts de compilateurs GPU rencontrés (reproducteurs minimaux)

| Fichier | Compilateur | Défaut | Contournement dans le proto |
|---|---|---|---|
| `signed-loop-trip-count.cu` | nvc++ 25.11 (`-cuda`, `-stdpar=gpu`) | boucle à compteur signé et test `!=` : N+1 itérations sur GPU à -O2/-O3 | boucles non signées dans le code concerné |
| `empty-kernel-stdpar.cxx` | nvc++ 25.11 | mode mémoire unifiée (`-gpu=mem:unified`, défaut) : écriture à une adresse de pile | `-gpu=mem:managed` |
| `aggregate-zero-init.cxx` | nvc++ 25.11 `-stdpar=gpu` | initialiseur par défaut `= {}` d'un membre tableau : constante LLVM invalide (« parse integer constant must have integer type ») | membre initialisé dans le constructeur (`PointView`) |
| `member-function-pointer.cxx` | nvc++ 25.11 `-stdpar=gpu` | appel par pointeur sur fonction membre ignoré **silencieusement** sur GPU (-O0 à -O3) | objets fonctions nommés générés par MFront |
| `std-variant-visit.cxx` | nvc++ 25.11 `-stdpar=gpu` | `std::variant` (libstdc++) : « Unexpected branch type », alias LLVM invalides | `TrivialGenTypeBase` pour l'opérateur tangent |
| `nvcc-static-constexpr-object.cu` | nvcc 13.0 | copie device d'un objet `static constexpr` local remplie de **zéros, sans diagnostic** si un constructeur est appelé | variables `__device__ constexpr` d'espace de noms (`tfel::utilities::getConstantObject`) |
| `nvcc-decltype-functional-cast.cu` | nvcc 13.0 (cudafe++) | `auto i = decltype(s){};` avec `const auto s = c.size();` dans un patron utilisé en code device : réécrit pour le compilateur hôte en une conversion fonctionnelle invalide (g++ : « expected primary-expression before ‘auto’ ») | `size_type{}` (`LinearInterpolation.ixx`) |
| `nvcc-builtin-constant-p.cu` | nvcc 13.0 | fonction `constexpr` utilisant `__builtin_constant_p` (fonction hôte pour nvcc) compilée **sans diagnostic** en une fonction vide dans le code device (-O0 à -O3) : `std::advance`, `std::next` et `std::prev` de libstdc++ 14 ne font rien sur GPU (interpolation par spline cubique à NaN) | `advance` et `previous` internes (`CubicSpline.ixx`) |
| (journal, 06/10 nuit) | nvc++ 25.11 à 26.9 `-mp=gpu` (OpenMP, -O2 et -O3) | lois des DSL implicites : `computeConsistentTangentOperator` renvoie faux sans passer par aucun de ses `return false`. Déclencheur : objet local généré par MFront (`GetPartialJacobianInvert`) tenant une référence `bool&` vers un drapeau local | objets de MFront propriétaires de leur drapeau (commit f008859df) |
| (journal, 06/10 nuit) | nvc++ 25.11 et 26.3 `-mp=gpu` (OpenMP) | bloc statique `_st_<hash>__T1_0` (317 octets) introuvable au chargement : texte `__PRETTY_FUNCTION__` d'un `assert` du calcul des vecteurs propres 3D, supprimé du code device mais enregistré | `assert` redondants retirés (commit a22a3cff6) |
| (journal, 06/10 nuit) | nvc++ 25.11 `-mp=gpu` (-O2) et nvcc 13.0 | contrainte écrite par une fonction d'initialisation de loi en grandes déformations perdue quand les champs sont rangés point par point. Déclencheur : sauvegarde puis restauration des pointeurs des données de la loi autour de l'appel, dans les pilotes de l'interface `generic` | pilotes travaillant sur une copie locale des données (branche `batchdata-by-variable`) |
| `omp-array-value-init.cxx` | nvc++ 25.11 `-mp=gpu` (OpenMP, -O2) | initialisation à zéro d'un membre tableau (`v()`, `v{}` ou `= {}`) d'un objet du noyau : copie d'un tableau global `.common` non initialisé, puis l'optimiseur **supprime le code du noyau** (aucun résultat écrit) | mise à zéro par une boucle dans le constructeur par défaut de `GenericFixedSizeArray` |
| `omp-lambda-pointer-array.cxx` | nvc++ 25.11 `-mp=gpu` (OpenMP) | lambda capturant par valeur un `std::array` de pointeurs, transmise par `map(to : kernel)` : pointeurs faux sur GPU à -O2 et -O3, **sans diagnostic** ; erreur NVVM à -O0 et -O1. Un pointeur capturé seul, une structure tenant un pointeur ou un tableau d'entiers fonctionnent | pointeurs des arguments capturés un par un (propriétés matériau) |
| `stdpar-staticinit-link.cxx` | nvc++ 25.11 `-stdpar=gpu` (-O0 à -O3) | données constantes du code device exportées sous un nom qui n'est unique que dans un fichier (`__staticinit_<compteur>_<fonction>`) : deux fichiers qui donnent le même compteur à la même fonction inline ne peuvent pas être liés ensemble (« nvlink error : Multiple definition »). Survient au hasard dans une bibliothèque de plusieurs lois | code device de chaque fichier lié seul (`-gpu=nordc`) |
| `omp-tensor-view-switch.cxx` (non minimal : en-têtes de TFEL où la vue `tensor` d'un `stensor` utilise un `switch`) | nvc++ 25.11 `-mp=gpu` (OpenMP, -O2 et -O3) | erreur interne « Computed goto control expression must be a 32-bit integer » sur le produit d'un `tensor` par un `stensor` : `switch` de la vue inliné dans celui du produit | vue `tensor` d'un `stensor` écrite sans `switch` (table d'indices) |
| `nvcc-cicc-initialize-function-reproducer.cu` (non minimal : inclut le début de la source générée de `SaintVenantKirchhoffElasticity` pour l'interface `generic-parallel`, compiler avec `-DUSE_BACKEND -I include` dans un dossier de génération de la loi) | nvcc 13.0, optimiseur NVVM (`cicc`) | la contrainte écrite par une fonction d'initialisation de loi en grandes déformations est perdue quand l'état de fin de pas est copié après la construction de la vue du point ; juste avec `-G`, `-Xcicc -O0` et clang CUDA ; le résultat change avec des `printf` | copie de l'état de fin de pas dans le constructeur de `PointView` |

`run-nvcxx-reproducers.sh <nvc++>` compile et exécute les reproducteurs de
nvc++ de -O0 à -O3. Résultats (06/10/2026, RTX 5080) :

| Reproducteur | 25.11 | 26.3 | 26.9 |
|---|---|---|---|
| `signed-loop-trip-count.cu` (-O2, -O3) | échec | échec | échec |
| `empty-kernel-stdpar.cxx` (-O2, -O3) | échec | échec | échec |
| `aggregate-zero-init.cxx` | échec | échec | **corrigé** |
| `member-function-pointer.cxx` | échec | échec | échec |
| `std-variant-visit.cxx` | échec | échec | échec |

Les deux défauts de nvcc (`nvcc-static-constexpr-object.cu`,
`nvcc-decltype-functional-cast.cu`) sont toujours présents avec nvcc 13.0.88,
13.1.115, 13.2.86, 13.3.73 et 13.4.92 (vérifié le 06/10/2026 avec les
compilateurs distribués en paquets pip `nvidia-cuda-nvcc`,
`nvidia-cuda-crt`, `nvidia-nvvm`). Les notes de version de NVHPC 26.1 à 26.9
ne donnent aucune liste de défauts corrigés : seuls les reproducteurs disent
si les défauts de nvc++ le sont.

Contraintes (non des défauts) : nvc++ `-stdpar=gpu` (mémoire gérée)
interdit l'accès aux variables globales et statiques depuis le code GPU
(erreur explicite) : `getConstantObject` y renvoie une copie.

nvcc accepte une définition sans `__host__ __device__` quand la
déclaration en porte (et inversement), clang CUDA la refuse. Ce n'est pas
un défaut de nvcc mais une incohérence de TFEL, corrigée dans 85
définitions (balayage clang de tous les en-têtes, voir le journal).

Défauts d'édition de liens et d'environnement (NVHPC 25.11) :

- `libnvcpumath.so` utilise `sincos` de `libm.so.6` sans en dépendre : un
  programme qui ne charge pas la bibliothèque mathématique ne peut pas
  charger une bibliothèque compilée par nvc++
  (`dlopen-without-libm.cxx`). Sans effet sur les programmes C++ qui
  utilisent `<cmath>`, comme MGIS.
- `nvc++ -mp=gpu -gpu=mem:managed` ne lie pas `libnvhpcmanaux`, dont la
  fonction `free` reconnaît la mémoire gérée, alors que `-stdpar=gpu` la
  lie. Un programme compilé ainsi qui charge des bibliothèques compilées par
  g++ peut s'arrêter sur `free(): invalid pointer` (lecture d'un fichier par
  `CxxTokenizer` dans une loi de polycristal). Contournement : lier le
  programme avec `-lnvhpcmanaux`.
