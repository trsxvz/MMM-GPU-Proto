# Thermique stationnaire par noyaux (10/10/2026)

Thermique stationnaire non linéaire de la poutre (`cases/thermal`, r5,
279 873 inconnues, 1,57 million de tétraèdres, 5 pas), conductivité
fonction de la température, température imposée aux deux bouts. Script :
`cases/measure-L7-thermal.sh`. Une RTX 5080, 16 cœurs, charge 2,4 au
départ. mfem-mgis 6cae30af99 plus le calcul des jacobiens par un noyau
(8b0b47acfb, `summary.txt` indique le commit de base).

La tangente thermique est définie négative et la diagonale des
températures imposées vaut un : GMRES remplace le gradient conjugué
(GMRES de hypre et BoomerAMG, tolérance 1e-8, ou Eisenstat-Walker 1 ;
GMRES de MFEM et Jacobi sans matrice).

| Cas | Avant (0542) | Après | Itérations de Newton |
|---|---|---|---|
| CPU 8 processus, point par point | 9,17 s | 9,03 s | 18 |
| CPU 8 processus, noyaux | 4,89 s | 4,76 s | 18 |
| GPU, noyaux | 3,01 s | **1,25 s** | 18 |
| GPU, noyaux, Eisenstat-Walker 1 | 4,03 s | 1,36 s | 25 |
| GPU, noyaux, sans matrice | 183,40 s | 179,09 s | 18 |

Temps de la boucle sur les pas. Même norme de la solution dans tous les
cas (306 232,626).

Détail (`cases/summary-profiling.py`, secondes et nombre d'appels) :

| Cas | Intégration | Résidu | Tangente | Résolution |
|---|---|---|---|---|
| CPU 8 processus, point par point | 2,11 (23) | 2,22 (23) | 3,05 (18) | 1,60 (18) |
| CPU 8 processus, noyaux | 0,50 (23) | 0,23 (23) | 2,45 (18) | 1,52 (18) |
| GPU, noyaux, avant | 2,03 (23) | 0,03 (23) | 0,34 (18) | 0,54 (18) |
| GPU, noyaux | 0,27 (23) | 0,03 (23) | 0,35 (18) | 0,53 (18) |

## Constats

- **Le GPU fait ×3,8 contre 8 processus avec les noyaux, ×7,2 contre le
  point par point.** Avant, la géométrie des éléments coûtait 1,76 s sur
  l'hôte à la première intégration (espace d'éléments finis d'ordre un et
  projection des coordonnées) : 60 % du temps sur le GPU. Le noyau des
  jacobiens la ramène à quelques millisecondes.
- **Sur CPU, la tangente domine** (2,45 s sur 4,76 s avec les noyaux,
  3,05 s point par point) : le profilage ne sépare pas encore le noyau des
  matrices élémentaires, l'assemblage creux et la matrice parallèle de
  hypre. C'est la mesure suivante sur CPU.
- **Eisenstat-Walker ne paie pas ici** : le Newton converge en 3 ou 4
  itérations, quadratiquement ; le Newton inexact demande 25 itérations au
  lieu de 18 et chaque résolution est déjà courte.
- **Sans matrice, Jacobi ne suffit pas** : 178,7 s de résolution, 10 s
  par résolution, contre 0,03 s avec BoomerAMG. Il faudrait un
  préconditionneur construit une fois (AMG sur une tangente de référence,
  comme NewtonK0 en mécanique).
