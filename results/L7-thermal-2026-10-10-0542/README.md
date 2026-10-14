# Thermique stationnaire par noyaux, première mesure (10/10/2026)

Thermique stationnaire non linéaire de la poutre (`cases/thermal`, r5,
279 873 inconnues, 5 pas), conductivité fonction de la température,
température imposée aux deux bouts. mfem-mgis 6cae30af99. Script :
`cases/measure-L7-thermal.sh`. Une RTX 5080, 16 cœurs, charge 1,9.

| Cas | Temps | Itérations de Newton |
|---|---|---|
| CPU 8 processus, point par point | 9,17 s | 18 |
| CPU 8 processus, noyaux | 4,89 s | 18 |
| GPU, noyaux | 3,01 s | 18 |
| GPU, noyaux, Eisenstat-Walker 1 | 4,03 s | 25 |
| GPU, noyaux, sans matrice (GMRES + Jacobi) | 183,40 s | 18 |

Sur le GPU, la première intégration coûtait 1,76 s : calcul unique de la
géométrie des éléments (`computeElementsGeometry`), fait sur l'hôte par un
espace d'éléments finis d'ordre un et une projection des coordonnées
(1,0 s et 0,83 s, mesurés à part). Ce coût est supprimé par
le calcul des jacobiens par un noyau : voir
`results/L7-thermal-2026-10-10-0606/`.
