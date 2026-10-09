# Mesures du lot L3 sur le polycristal (09/10/2026)

Polycristal de mm-opera-hpc (5 grains, loi de monocristal de l'UO2,
grandes transformations, orthotrope, problème périodique), 10 pas, sur
CPU, bibliothèque `generic`. Gradients, forces intérieures et matrices
élémentaires de la tangente calculés élément par élément ou par noyaux
(mfem-mgis ed91564643). Script : `cas/mesures-L3.sh polycristal`.

| Cas | Total | Intégration | Résidu | Tangente | Résolution |
|---|---|---|---|---|---|
| 1 processus, point par point | 149,73 | 104,52 | 4,16 | 8,53 | 32,20 |
| 1 processus, lots | 148,03 | 102,74 | 4,16 | 8,34 | 32,47 |
| 1 processus, gradients par noyaux | 144,97 | 99,52 | 4,19 | 8,43 | 32,51 |
| 1 processus, résidu par noyaux | 145,46 | 102,55 | 1,49 | 8,61 | 32,49 |
| 1 processus, tangente par noyaux | 151,04 | 102,70 | 4,14 | 11,50 | 32,38 |
| 1 processus, tout par noyaux | 145,37 | 99,55 | 1,49 | 11,53 | 32,47 |
| 8 processus, point par point | 30,35 | 14,45 | 0,59 | 1,69 | 13,38 |
| 8 processus, lots | 30,21 | 14,31 | 0,58 | 1,68 | 13,40 |
| 8 processus, tout par noyaux | 29,89 | 13,85 | 0,21 | 2,16 | 13,43 |

Contraintes macroscopiques : identiques au calcul point par point sur le
même nombre de processus à 5e-9 près (rapporté au maximum de chaque
colonne), au niveau de la tolérance du solveur itératif.

Constats : le résidu par noyaux est 2,8 fois plus rapide ; les gradients
par noyaux font gagner 3 s sur 103 (la loi domine : voir L2 pour ses fils
d'exécution) ; les matrices élémentaires par noyaux restent plus lentes
sur CPU (+38 %), comme sur la poutre.
