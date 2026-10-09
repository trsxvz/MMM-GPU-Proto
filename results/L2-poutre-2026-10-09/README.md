# Mesures du lot L2 sur la poutre (09/10/2026)

Poutre r4 (111 843 inconnues), plasticité de von Mises à écrouissage
linéaire, PCG + BoomerAMG (tolérance 1e-4), 10 pas. mfem-mgis 4a68166ca9,
loi `generic-parallel` compilée avec les options corrigées (patch 0093 de
TFEL), machine calme. Remplace les mesures provisoires de
`results/L2-2026-10-09/` (machine chargée, fils mal placés).

| Cas | Total | Intégration | Résidu | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|---|
| 8 processus, point par point | 15,62 | 0,49 | 0,43 | 2,22 | 9,55 | 2,87 |
| 8 processus, lots, `generic` | 15,56 | 0,50 | 0,43 | 2,21 | 9,52 | 2,84 |
| 8 processus, lots, `openmp` 1 fil | 15,53 | 0,49 | 0,43 | 2,21 | 9,49 | 2,83 |
| 1 processus, point par point | 65,36 | 3,39 | 2,92 | 11,63 | 33,73 | 13,50 |
| 1 processus, lots, `generic` | 65,06 | 3,19 | 2,91 | 11,55 | 33,63 | 13,56 |
| 1 processus, lots, `openmp` 1 fil | 65,10 | 3,22 | 2,95 | 11,61 | 33,65 | 13,46 |
| 1 processus, lots, `openmp` 8 fils | 64,82 | 2,63 | 2,94 | 11,55 | 33,96 | 13,49 |
| 1 processus, lots, `stlpar` 8 cœurs | 65,11 | 2,68 | 2,95 | 11,57 | 34,15 | 13,54 |
| 2 processus, point par point | 36,62 | 1,67 | 1,53 | 6,53 | 19,17 | 7,57 |
| 2 processus, lots, `openmp` 4 fils | 36,34 | 1,35 | 1,56 | 6,54 | 19,26 | 7,49 |

Même solution partout (norme 48,83704586).

Constat : la loi pèse peu (moins de 5 % du temps) ; plusieurs fils ne
réduisent que sa part, pas le calcul des gradients élément par élément
qui reste dans la phase d'intégration. Les gains de la poutre viendront de
la résolution (plus de 50 %), de la tangente et de la prédiction (lots L3
à L5).
