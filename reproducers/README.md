# Reproducteurs

Petits programmes autonomes qui vérifient un point précis, en dehors des
grandes suites de tests. Chaque dossier a un script `run.sh` qui travaille
dans `scratch/reproducers` et utilise MFront de l'arbre de développement
(`scratch/build-dev2`, sources de `tfel-dev`).

| Dossier | Ce qui est vérifié | Résultat attendu |
|---|---|---|
| `rotations` | Rotation des blocs de l'opérateur tangent d'une loi orthotrope (bibliothèque `generic`), indépendamment du code qui les tourne : pour tout opérateur `K` et tous gradients `x`, `rotateTangentOperatorBlocks(K) x = rotateThermodynamicForces(K rotateGradients(x))`, et rotation sur place identique. La loi utilisée par défaut, `OrthotropicTangentOperatorBlocksTest` (tests `generic-parallel` de `tfel-dev`), a des blocs de tous les types (tableaux, vecteurs, tenseurs symétriques, tenseurs). Une autre loi peut être donnée en argument. | `OK` pour chaque hypothèse (écart de l'ordre de 1e-15) |
| `device-strings` | Constante chaîne utilisée dans du code GPU (comparaisons, taille), avec six compilateurs. `std::string_view` passe par `std::char_traits` (`strlen`, `memcmp`) : résultats faux sans message avec nvcc et hipcc, erreurs avec nvc++. Un petit type sans `char_traits` (`constant-string.hxx`) fonctionne partout. | variantes A et B fausses ou refusées, C et D justes (`1 7 13 19`) |
| `sync-call-overhead` | Durée d'un appel synchrone du backend CUDA sur un lot d'une loi élastique, données en mémoire managée, avec les en-têtes de `tfel-dev` à un commit de référence (avant R8 par défaut) et ceux de l'arbre de travail. | 09/10, RTX 5080 : 1000 points 135 → 29 µs, 100000 points 546 → 442 µs |

Les comparaisons entre les interfaces `generic` et `generic-parallel` des
lois touchées se lancent seules, sans construire toutes les lois de test,
avec `scratch/targeted-test.sh configuration motif loi...`.
