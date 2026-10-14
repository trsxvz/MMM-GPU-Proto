# Cermet en float : seuil du point fixe sur les contraintes latérales (14/10/2026)

Cermet de mm-opera-hpc (test `CermetNullStrainXXYY` : maillage
`5grains.msh`, deux pas de 0,2 s, deux processus, sans post-traitement).
Constructions de la CI simulée en local (`build/ci-sim-clean` en double,
`build/ci-sim-float` en float, lois compilées par gcc 14 avec les options
de `tfel-config --oflags`) et construction float du prototype
(`build/mm-opera-hpc-float`, lois compilées par gcc 13). Fichiers `.res` :
temps, F11, F22, F33, puis σ11, σ22, σ33 (Pa) à la fin de chaque pas.
`point-fixe.txt` : résidu du point fixe sur les contraintes latérales,
|res| = (σ11² + σ22²)^½, à chaque itération.

| Calcul | σ33 à 0,2 s (Pa) | σ33 à 0,4 s (Pa) | Écart au double (pas de 0,2 s) |
|---|---|---|---|
| Double, pas de 0,2 s | 17 219 167 | 29 377 289 | |
| Double, pas de 0,1 s | 17 642 604 | 29 886 489 | +2,46 % et +1,73 % |
| Float, gcc 13, seuil 1e4 Pa | 17 245 144 | 29 393 436 | +0,15 % et +0,05 % |
| Float, gcc 14, seuil 1e4 Pa | 17 690 582 (pas coupé) | 29 498 396 | +2,74 % et +0,41 % |
| Float, gcc 13, seuil 5e4 Pa | 17 260 680 | 29 393 780 | +0,24 % et +0,06 % |
| Float, gcc 14, seuil 5e4 Pa | 17 262 862 | 29 397 298 | +0,25 % et +0,07 % |
| Cast3M (référence du test) | 17 246 000 | 29 462 000 | |

- Avec le seuil de 1e4 Pa et la loi compilée par gcc 14, le point fixe du
  premier pas oscille entre 17 234,6 et 15 454,4 Pa pendant ses 20
  itérations : le pas est rejeté puis coupé en deux pas de 0,1 s. La
  contrainte à 0,2 s s'écarte alors du double de 2,74 %, mais de 0,27 %
  seulement du double calculé avec les mêmes pas de 0,1 s : l'écart vient
  du pas de temps, pas de l'arrondi.
- Cause de l'oscillation : les composantes latérales de F, juste sous 1,
  ne bougent que d'un ulp (6e-8) en float, soit environ 1e4 Pa de
  contrainte latérale pour ce module. Le point fixe ne peut pas
  descendre de façon sûre sous ce pas.
- Avec un seuil de 5e4 Pa en float (mm-opera-hpc 2e4881c), les pas de
  0,2 s sont acceptés avec les deux compilateurs, qui ne diffèrent plus
  que de 0,013 % ; l'écart à Cast3M reste sous 0,25 %. Le test vérifie de
  nouveau les deux instants en float. Le double ne change pas (mêmes
  fichiers au bit près).
