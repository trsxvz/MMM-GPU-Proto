# Tangente par noyaux de bout en bout, 8 processus (11/10/2026)

Script : `cases/measure-positions.sh`. Thermique stationnaire r5 (279 873
inconnues, 5 pas) et poutre r4 (10 pas), 8 processus d'un fil, GMRES ou
gradient conjugué de hypre avec BoomerAMG. Trois calculs par cas : point
par point, par noyaux avec mfem-mgis installé dans `install/mfem-mgis`
(8a6325335e, avant), par noyaux avec mfem-mgis installé dans
`install/mfem-mgis-registres` (4a007e85f2 : positions des entrées dans la
matrice creuse et matrices sans former B, après). La ligne « mfem-mgis »
de `summary.txt` donne la version des sources, pas celle des installations.
Machine calme, GPU inutilisé.

Temps en secondes (`cases/summary-profiling.py`, nombre d'appels entre
parenthèses) :

| Cas | Total | Intégration | Résidu | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|---|
| poutre, point par point | 15,70 | 0,51 | 0,45 | 2,25 (31) | 9,56 | 2,87 |
| poutre, noyaux avant | 13,05 | 0,18 | 0,08 | 1,40 (29) | 8,76 | 2,57 |
| poutre, noyaux après | 12,25 | 0,18 | 0,08 | 0,81 (29) | 8,71 | 2,42 |
| thermique, point par point | 9,44 | 2,15 | 2,30 | 3,12 (18) | 1,59 | |
| thermique, noyaux avant | 4,21 | 0,51 | 0,23 | 1,64 (18) | 1,56 | |
| thermique, noyaux après | 3,65 | 0,51 | 0,22 | 1,06 (18) | 1,59 | |

- La tangente par noyaux gagne ×1,73 (poutre) et ×1,55 (thermique) ; elle
  est maintenant 2,8 à 2,9 fois plus rapide que la tangente point par
  point. Le reste de la phase est surtout le produit PtAP de hypre.
- Totaux : ×1,07 (poutre, la résolution domine) et ×1,15 (thermique).
- Mêmes itérations de Newton et mêmes normes de la solution que le calcul
  avant ces lots.
