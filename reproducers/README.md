# Reproducteurs

Petits programmes autonomes qui vérifient un point précis, en dehors des
grandes suites de tests. Chaque dossier a un script `run.sh` qui travaille
dans `scratch/reproducers` et utilise MFront de l'arbre de développement
(`scratch/build-dev2`, sources de `tfel-dev`).

| Dossier | Ce qui est vérifié | Résultat attendu |
|---|---|---|
| `rotations` | Rotation des blocs de l'opérateur tangent d'une loi orthotrope (bibliothèque `generic`), indépendamment du code qui les tourne : pour tout opérateur `K` et tous gradients `x`, `rotateTangentOperatorBlocks(K) x = rotateThermodynamicForces(K rotateGradients(x))`, et rotation sur place identique. La loi utilisée par défaut, `OrthotropicTangentOperatorBlocksTest` (tests `generic-parallel` de `tfel-dev`), a des blocs de tous les types (tableaux, vecteurs, tenseurs symétriques, tenseurs). Une autre loi peut être donnée en argument. | `OK` pour chaque hypothèse (écart de l'ordre de 1e-15) |
| `device-strings` | Constante chaîne utilisée dans du code GPU (comparaisons, taille), avec six compilateurs. `std::string_view` passe par `std::char_traits` (`strlen`, `memcmp`) : résultats faux sans message avec nvcc et hipcc, erreurs avec nvc++. Un petit type sans `char_traits` (`constant-string.hxx`) fonctionne partout. | variantes A et B fausses ou refusées, C et D justes (`1 7 13 19`) |
| `hypre-ilu-gpu` | Factorisation incomplète de hypre (ILU) préconditionnant son GMRES flexible, sur CPU et sur GPU (hypre 3.2.0 construit avec CUDA), trois résolutions, l'opérateur étant redonné avant chacune. `mfem::HypreILU` recrée le préconditionneur de hypre avec les options par défaut de MFEM (niveau de remplissage 1) à chaque nouvel opérateur : un niveau donné une fois (mode 2) est perdu après la première résolution. Un préconditionneur gardé d'un opérateur à l'autre (mode 1) garde son niveau. Programme autonome (MFEM et hypre seuls), `run.sh` sans MFront. | CPU, 2 processus : 17 itérations avec le niveau 1, niveau 0 donné une fois 18 puis 17, 17 (niveau perdu). GPU : ILU(0) gardé 18 itérations à chaque résolution ; niveau 1 gardé ou par défaut : erreur de segmentation dans `hypre_ILUSolveLUDevice` dès la première résolution ; niveau 0 donné une fois : 18 puis erreur de segmentation à la deuxième |
| `amgx-classical-cuda13` | Exemple `amgx_capi` d'AmgX (gradient conjugué préconditionné par son multigrille, configurations livrées avec AmgX) sur un Laplacien 20³, avec l'AmgX de l'environnement (commit 91a8413e de la branche principale, CUDA 13.0.2, sm_120) et, s'il est construit, celui de `build/amgx-stub` (même commit, `-static-global-template-stub=false`). Le multigrille classique lance un noyau absent du module de l'objet qui le lance : `amg.cu.o` garde une copie faible de `cub::DeviceFor::Bulk` (appelé par `thrust::sequence` depuis `Matrix::colorMatrix`) dont le stub du noyau est interne (défaut de CUDA 13) alors que son cubin ne contient pas le noyau. | AmgX de l'environnement : classique « invalid device function » dans `parallel_for`, agrégation 22 itérations |
| `sync-call-overhead` | Durée d'un appel synchrone du backend CUDA sur un lot d'une loi élastique, données en mémoire managée, avec les en-têtes de `tfel-dev` à un commit de référence (avant R8 par défaut) et ceux de l'arbre de travail. | 08/10, RTX 5080 : 1000 points 135 → 29 µs, 100000 points 546 → 442 µs |

Les comparaisons entre les interfaces `generic` et `generic-parallel` des
lois touchées se lancent seules, sans construire toutes les lois de test,
avec `scratch/targeted-test.sh configuration motif loi...`.
