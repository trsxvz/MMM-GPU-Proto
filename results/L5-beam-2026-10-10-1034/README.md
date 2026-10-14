# Options de BoomerAMG pour NewtonK0 (10/10/2026)

Poutre r4 du prototype (111 843 inconnues, 10 pas), tous les noyaux,
mfem-mgis a066c82baa. Meilleure méthode des mesures précédentes : NewtonK0
(gradient conjugué préconditionné par BoomerAMG construit une fois sur K0,
prédiction par ce gradient conjugué) avec Eisenstat-Walker 1, tolérance
initiale 0,1. Options de hypre passées par `-amgo`. Script :
`cases/measure-L5-amg.sh`. Une RTX 5080 (loi `cuda`), 16 cœurs (loi
`generic`, 8 processus). Charge 0,2 au départ.

| Options de BoomerAMG | GPU | Itérations | CPU, 8 processus | Itérations |
|---|---|---|---|---|
| par défaut | 1,85 s | 33 | 8,78 s | 38 |
| `StrengthThreshold=0.5` | 1,87 s | 35 | 8,80 s | 38 |
| `StrengthThreshold=0.7` | 2,55 s | 38 | 10,05 s | 38 |
| `AggressiveCoarseningLevels=1` | **1,71 s** | 34 | 7,41 s | 34 |
| `StrengthThreshold=0.5`, `AggressiveCoarseningLevels=1` | 1,72 s | 37 | **7,22 s** | 34 |
| `CoarseningType=10` (HMIS) | 1,87 s | 38 | 9,74 s | 38 |
| `InterpolationType=14` | 1,79 s | 36 | 9,27 s | 35 |
| `RelaxationType=16` (Chebyshev) | 2,31 s | 37 | 9,15 s | 35 |
| `RelaxationSweeps=2` | 2,25 s | 33 | 10,45 s | 38 |
| `RelaxationType=7` (Jacobi) | échec | | échec | |

Temps de la boucle sur les pas, une exécution par cas ; itérations de
Newton au total.

Trois exécutions de plus sur le GPU (machine moins calme, charge 4) :
par défaut 1,91, 2,05, 1,91 s ; grossissement agressif 1,87, 1,86, 1,73 s ;
seuil 0,5 et grossissement agressif 1,73, 1,81, 1,72 s ; interpolation 14
1,93, 1,88, 1,82 s.

## Constats

- **Le grossissement agressif sur un niveau est le meilleur réglage** :
  -16 à -18 % sur 8 processus CPU (7,22 s avec le seuil 0,5, contre
  8,78 s), -5 à -8 % sur le GPU. Il réduit le coût de chaque cycle de
  BoomerAMG (hiérarchie plus petite) sans augmenter le nombre de
  résolutions.
- Les lisseurs plus coûteux (Chebyshev, deux balayages) et un seuil de
  connexion de 0,7 ralentissent ; le Jacobi simple de hypre (7) fait
  échouer le gradient conjugué dès le premier pas (préconditionneur qui
  n'est plus symétrique défini positif).
- Par rapport à la mesure de 05h (`results/L5-beam-2026-10-10-0429`,
  NewtonK0 + EW1 tolérance initiale 0,1 : 1,96 s sur GPU, 9,46 s sur 8
  cœurs), les noyaux de la géométrie et des matrices sur l'hôte donnent
  1,85 s et 8,78 s avec les options par défaut.
- Sur le GPU, la résolution linéaire domine toujours (0,94 à 1,07 s sur
  1,7 à 2,1 s, `cases/summary-profiling.py`), puis la prédiction (0,42 s).
