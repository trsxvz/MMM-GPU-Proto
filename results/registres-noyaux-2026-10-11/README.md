# Registres des noyaux GPU, matrices de rigidité sans former B (11/10/2026)

## Constat

`cuobjdump -res-usage` sur les objets des noyaux de mfem-mgis (CUDA 13.0,
sm_120). Colonnes : fichier, registres par fil, pile (octets), mémoire
locale (octets), noyau. `kind=N` : `GradientsKind` dans l'ordre de sa
déclaration (0 petites déformations 3D, 1 petites déformations 2D, 2
gradient de la transformation 3D, 3 gradient de la transformation 2D, 4 et
5 gradient de température, 6 et 7 température, 8 et 9 micromorphe).

- `registres.txt` : branche `positions` (b64d5dda96). Les matrices de
  rigidité en grandes transformations 3D (kind 2, cas ssna303 et
  polycristal) prenaient 254 registres par fil, le maximum, et 576 octets
  de pile : au plus 8 warps par multiprocesseur. En petites déformations
  3D : 122 à 126 registres, 576 octets de pile. La pile venait des
  tableaux `Ba[9][3]`, `Bb[9][3]`, alors que `B` est creux (en grandes
  transformations, `F_ij` ne dépend que de `u_i`, par `dN/dX_j`).
- `registres-apres.txt` : branche `registres` (3570cd54dd, premier
  commit), produits
  `K B_b` et `B_aᵀ (K B_b)` écrits avec une table constexpr des entrées non
  nulles de `B` (`getBOperatorEntries`, vérifiée à la compilation contre
  `computeB`), dans le même ordre d'opérations. Grandes transformations
  3D : 64 registres, 72 octets de pile ; petites déformations 3D : 70 à 80
  registres, 72 octets ; les cas 2D n'ont plus de pile. Le gradient de
  température 3D et le micromorphe 3D passent de 64 et 72 à 84 registres.

## Mesures

Banc `cas/assemblage` (assemblage par noyaux de toutes les matrices de
rigidité dans la matrice creuse locale, blocs de la tangente arbitraires,
un processus), construit contre les installations séparées des branches
`positions` (d04acab66c, avant) et `registres` (28842e74dd, après), mêmes
commandes. Machine calme, un cœur (`taskset -c 4`) ou la RTX 5080. Temps
médians de 10 assemblages. Sorties : `hote-2.txt`, `gpu-2.txt` (les
fichiers `hote.txt` et `gpu.txt` sont une première série où le temps
mesuré sur le GPU comptait la copie de toute la matrice vers l'hôte).

Un cœur (avec les positions des entrées, chemin par défaut sans
périphérique) :

| Cas | Avant | Après | Gain | Gain total |
|---|---|---|---|---|
| grandes transformations, hexaèdres Q1 27 points (r4) | 164 ms | 52 ms | ×3,1 | ×3,4 |
| petites déformations, hexaèdres Q1 27 points (r4) | 116 ms | 48 ms | ×2,4 | ×2,7 |
| poutre, tétraèdres P1 (r4) | 116 ms | 75 ms | ×1,5 | ×3,1 |
| grandes transformations, tétraèdres P2 (r3) | 277 ms | 133 ms | ×2,1 | ×3,0 |
| thermique, tétraèdres P1 (r5) | 137 ms | 121 ms | ×1,1 | ×3,2 |

« Gain total » : par rapport à l'assemblage sans positions de la branche
`positions`, c'est-à-dire à l'état de gpu-proto avant ces deux lots.
Valeurs des matrices identiques au bit près entre les deux branches
(fichiers écrits par `--output-values`).

GPU (sans positions, chemin par défaut avec un GPU) :

| Cas | Avant | Après | Gain |
|---|---|---|---|
| grandes transformations, hexaèdres Q1 27 points (r5) | 59,8 ms | 33,0 ms | ×1,81 |
| petites déformations, hexaèdres Q1 27 points (r5) | 42,5 ms | 30,5 ms | ×1,40 |
| poutre, tétraèdres P1 (r5) | 25,0 ms | 18,1 ms | ×1,38 |
| grandes transformations, tétraèdres P2 (r4) | 102,2 ms | 52,9 ms | ×1,93 |
| thermique, tétraèdres P1 (r5) | 13,0 ms | 13,8 ms | ×0,94 |

Sur le GPU, les valeurs diffèrent d'un calcul à l'autre, avec le même
programme, à cause de l'ordre des additions atomiques : écart maximal
1,4e-17 pour des valeurs de 5,5e-2, le même entre les deux branches. La
thermique perd 6 % sur le GPU (84 registres au lieu de 64 pour le
gradient de température 3D).

Le noyau de l'action de la tangente sans assemblage (option `MatrixFree`)
utilise la même table (second commit) : 80 puis 72 registres en grandes
transformations 3D, résultats identiques ; non mesuré par ce banc.
