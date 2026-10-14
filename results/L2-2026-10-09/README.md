# Mesures du lot L2 sur la poutre (09/10/2026, provisoires)

**Remplacées** par `results/L2-beam-2026-10-09/` (machine calme, fils
placés correctement, options de compilation de la loi corrigées).

Intégration des lois par lots (option `BehaviourIntegration`, mfem-mgis
ef52adb445), poutre r4 (111 843 inconnues), PCG + BoomerAMG, tolérance
1e-4, 10 pas. Script : `cases/measure-L2.sh` ; résumé par phase :
`cases/summary-profiling.py *.txt`.

**Mesures provisoires** : une autre compilation tournait sur la machine
(charge de 5 à 9). Le temps de référence du lot L0 pour le même cas, sur
machine calme, est de 15,45 s à 8 processus ; on mesure ici 20,02 s pour
le même chemin (point par point). Les écarts de moins de 10 à 20 % ne sont
pas significatifs.

## Résultats

Toutes les variantes donnent la même solution (norme 48,83704586) et les
mêmes itérations de Newton (31 à 8 processus, 28 à 1 et 2 processus).

Temps par phase (secondes, maximum sur les processus ; entre parenthèses,
nombre d'appels). Intégration : `NS::processNewUnknownsEstimate`, qui
comprend le calcul des gradients élément par élément et l'appel de la loi.

| Cas | Total | Intégration | Résidu | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|---|
| 8 processus, point par point | 20,02 | 0,86 | 0,74 | 2,89 | 11,84 | 3,61 |
| 8 processus, lots, `generic` | 19,59 | 0,73 | 0,76 | 2,97 | 11,55 | 3,49 |
| 8 processus, lots, `openmp` 1 fil | 17,19 | 0,59 | 0,53 | 2,43 | 10,40 | 3,16 |
| 1 processus, point par point | 67,06 | 3,36 | 2,97 | 11,72 | 34,92 | 13,87 |
| 1 processus, lots, `generic` | 67,62 | 3,25 | 2,96 | 11,73 | 35,25 | 14,20 |
| 1 processus, lots, `openmp` 1 fil | 67,55 | 3,24 | 3,01 | 11,76 | 35,16 | 14,16 |
| 1 processus, lots, `openmp` 8 fils | 70,35 | 2,93 | 3,20 | 12,19 | 36,95 | 14,85 |
| 1 processus, lots, `stlpar` 8 cœurs | 69,02 | 3,37 | 3,03 | 11,96 | 36,07 | 14,38 |
| 2 processus, lots, `openmp` 4 fils | 37,48 | 1,40 | 1,63 | 6,71 | 19,96 | 7,66 |

## Constats

- **Exactitude** : le mode par lots ne change aucun résultat, quelle que
  soit la bibliothèque (`generic`, `generic-parallel` openmp ou stlpar) et
  le nombre de fils.
- **La poutre n'est pas le bon cas pour L2** : l'intégration de la loi
  (plasticité de von Mises à écrouissage linéaire) pèse 5 % du temps. La
  résolution linéaire (52 %), la prédiction (21 %) et l'assemblage de la
  tangente (17 %) dominent. Plusieurs fils ne réduisent que l'appel de la
  loi, pas le calcul des gradients, qui reste séquentiel élément par
  élément dans la même phase.
- **Plusieurs fils par processus** : sans effet visible ici (1 processus
  à 8 fils : 2,93 s d'intégration contre 3,24 s à 1 fil). À 8 processus
  d'un fil, la bibliothèque `generic-parallel` gagne 0,27 s sur 0,86 s
  d'intégration, dans le bruit de la machine chargée.
- **Placement** : avec `--bind-to none`, OpenMPI 5 attache tous les
  processus au premier cœur (premier essai abandonné : 8 processus sur un
  cœur). Le script réserve maintenant des cœurs par processus
  (`--map-by slot:PE=<fils> --bind-to core`).
- **Cas à mesurer ensuite** : le polycristal (55 % du temps dans
  l'intégration au lot L0, loi cristalline orthotrope). Il demande le
  support des lois orthotropes `generic-parallel` (rotations par champs),
  en cours.
