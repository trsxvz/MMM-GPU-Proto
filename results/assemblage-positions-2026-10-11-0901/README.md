# Positions des entrées dans la matrice creuse (11/10/2026)

Banc `cas/assemblage` : assemblage par noyaux des matrices de rigidité de
tous les éléments dans la matrice creuse locale de la tangente, avec les
positions des entrées calculées une fois
(`batch::computeStiffnessMatricesPositions`) et sans elles (recherche de
chaque entrée dans sa ligne, par dichotomie), dans le même programme, en
alternant les deux, 10 fois chacun. Blocs de la tangente arbitraires, un
processus. mfem-mgis branche `positions` (d04acab66c), construction
séparée installée dans `install/mfem-mgis-positions`. Machine calme, un
cœur (`taskset -c 4`) ou la RTX 5080.

Sorties : `hote-2.txt` (un cœur) et `gpu-2.txt` (`--device cuda`). Les
fichiers `hote.txt` et `gpu.txt` sont la première série (b64d5dda96) :
sur le GPU, chaque assemblage chronométré y comptait la copie de toute la
matrice vers l'hôte (temps gonflés de 10 à 40 %), et le temps des
positions la copie des positions ; les rapports sont les mêmes.

## Un cœur

| Cas | Éléments | Entrées de la matrice | Positions | Sans positions | Avec positions | Gain | Rentabilisé à |
|---|---|---|---|---|---|---|---|
| poutre, tétraèdres P1, mécanique (r4) | 196 608 | 4,7 M | 57 Mo, 419 ms | 235 ms | 117 ms | ×2,0 | 4 assemblages |
| thermique, tétraèdres P1 (r5) | 1 572 864 | 4,1 M | 50 Mo, 272 ms | 382 ms | 133 ms | ×2,9 | 2 |
| hexaèdres Q1, 27 points, petites déformations (r4) | 4 096 | 1,1 M | 4,7 Mo, 31 ms | 130 ms | 116 ms | ×1,1 | 3 |
| hexaèdres Q1, 27 points, grandes transformations (r4) | 4 096 | 1,1 M | 4,7 Mo, 31 ms | 178 ms | 165 ms | ×1,1 | 3 |
| tétraèdres P2, mécanique (r3) | 24 576 | 8,9 M | 44 Mo, 357 ms | 337 ms | 216 ms | ×1,6 | 3 |
| hexaèdres Q2, mécanique (r2) | 64 | 0,3 M | 0,8 Mo, 5 ms | 44 ms | 41 ms | ×1,1 | 2 |

Temps médians ; le temps des positions est celui d'un second calcul.
Matrices identiques au bit près, aucune entrée cherchée avec les
positions. Le calcul des positions coûte de 0,1 à 1,8 assemblage ; il est
rentabilisé entre la deuxième et la quatrième tangente. Le gain est grand
quand le calcul des blocs est léger (éléments linéaires, peu de points,
thermique) ; avec 27 points par hexaèdre, le calcul des blocs domine
(voir `results/registres-noyaux-2026-10-11`, qui le réduit). Mémoire : 2
n² octets par élément, du même ordre que la matrice locale pour les
éléments linéaires.

## GPU

| Cas | Éléments | Positions | Sans positions | Avec positions | Gain |
|---|---|---|---|---|---|
| poutre, tétraèdres P1 (r5) | 1 572 864 | 453 Mo | 24,6 ms | 24,7 ms | ×0,99 |
| thermique, tétraèdres P1 (r5) | 1 572 864 | 50 Mo | 13,0 ms | 13,0 ms | ×1,00 |
| hexaèdres Q1, 27 points (r5) | 32 768 | 38 Mo | 42,9 ms | 42,5 ms | ×1,01 |
| tétraèdres P2 (r4) | 196 608 | 354 Mo | 67,5 ms | 67,1 ms | ×1,01 |

Aucun gain sur le GPU : la recherche y est masquée par le calcul des blocs
et les additions atomiques, et les positions prendraient jusqu'à 450 Mo.
Décision : `MultiMaterialNonLinearIntegrator` ne les utilise pas avec un
GPU (ni avec le périphérique `debug` de MFEM). Il les garde pour chaque
disposition des lignes (colonnes triées par un solveur direct, ou
diagonale en tête), les recalcule après un assemblage où des entrées ont
été cherchées, au plus trois fois, puis les abandonne.

Le banc ne mesure que le noyau standard ; les noyaux du FBar (rang un) et
de Faltus lisent les mêmes positions.
