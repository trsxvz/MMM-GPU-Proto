# Cermet en float : seuil du point fixe sur les contraintes latérales (14/10/2026)

Cermet de mm-opera-hpc (test `CermetNullStrainXXYY` : maillage
`5grains.msh`, deux pas de 0,2 s, sans post-traitement). La loi est
passée par `--library` : bibliothèque float compilée par gcc 13
(construction float du prototype) ou par gcc 14 (CI simulée), dans les
deux cas avec les options de `tfel-config --oflags`. Exécutables :

- « CI » : cermet de la CI simulée en local (`build/ci-sim-float`,
  environnement `ci-cpu-float` ; `build/ci-sim-clean` en double), dans un
  environnement vide ;
- « prototype » : cermet de la construction float du prototype
  (`build/mm-opera-hpc-float`, environnement `mmm-float-cpu`).

Le seuil du point fixe sur les contraintes latérales, constante de
`UniaxialMacroscopicStressPeriodicSimulation.hxx`, est changé pour chaque
ligne, puis l'exécutable est reconstruit. Fichiers `<exécutable>-<n>p-<loi>-<seuil>.res`
(n processus) et `double-pas-*.res` : temps, F11, F22, F33, puis σ11, σ22,
σ33 (Pa) à la fin de chaque pas. `point-fixe.txt` : résidu du point fixe,
|res| = (σ11² + σ22²)^½, à chaque itération de chaque calcul.

Double, exécutable de la CI, un processus (le test) : σ33 = 17 219 167 Pa
à 0,2 s et 29 377 289 Pa à 0,4 s par pas de 0,2 s ; 17 642 604 Pa
(+2,46 %) et 29 886 489 Pa (+1,73 %) par pas de 0,1 s.

Float : σ33 à 0,2 s puis à 0,4 s, et écart au double par pas de 0,2 s.

| Exécutable | Processus | Seuil (Pa) | Loi compilée par gcc 13 | Loi compilée par gcc 14 |
|---|---|---|---|---|
| CI | 1 | 1e4 | 17 684 464 (pas coupé), +2,70 % ; 29 501 400, +0,42 % | 17 684 506 (pas coupé), +2,70 % ; 29 501 338, +0,42 % |
| CI | 1 | 2e4 | 17 247 006, +0,16 % ; 29 387 932, +0,04 % | 17 247 020, +0,16 % ; 29 387 976, +0,04 % |
| CI | 1 | 3e4 | 17 247 006, +0,16 % ; 29 396 894, +0,07 % | 17 247 020, +0,16 % ; 29 396 888, +0,07 % |
| CI | 1 | 5e4 | 17 262 850, +0,25 % ; 29 397 188, +0,07 % | 17 262 862, +0,25 % ; 29 397 298, +0,07 % |
| CI | 2 | 1e4 | 17 684 538 (pas coupé), +2,70 % ; 29 501 340, +0,42 % | 17 684 532 (pas coupé), +2,70 % ; 29 501 422, +0,42 % |
| CI | 2 | 2e4 | 17 246 848, +0,16 % ; 29 387 992, +0,04 % | 17 246 870, +0,16 % ; 29 387 964, +0,04 % |
| CI | 2 | 3e4 | 17 246 848, +0,16 % ; 29 396 982, +0,07 % | 17 246 870, +0,16 % ; 29 396 966, +0,07 % |
| CI | 2 | 5e4 | 17 262 702, +0,25 % ; 29 397 332, +0,07 % | 17 262 718, +0,25 % ; 29 397 140, +0,07 % |
| prototype | 2 | 1e4 | 17 245 144, +0,15 % ; 29 393 436, +0,05 % | 17 690 582 (pas coupé), +2,74 % ; 29 498 396, +0,41 % |
| prototype | 2 | 3e4 | 17 237 086, +0,10 % ; 29 393 648, +0,06 % | 17 237 132, +0,10 % ; 29 393 532, +0,06 % |

Référence de Cast3M du test : 17 246 000 Pa à 0,2 s et 29 462 000 Pa à
0,4 s, à 1 % près.

- Avec le seuil de 1e4 Pa, le point fixe du premier pas ne converge pas
  dans cinq calculs sur six : il stagne vers 10 700 à 10 830 Pa
  avec l'exécutable de la CI, quelle que soit la loi, et oscille entre
  17 234,6 et 15 454,4 Pa avec l'exécutable du prototype et la loi de
  gcc 14. Le pas est rejeté puis coupé en deux pas de 0,1 s. L'écart de
  2,7 % au double vient de ce pas de temps : par pas de 0,1 s, le double
  donne +2,46 %, et le float s'en écarte de 0,27 %. Le compilateur de la
  loi n'est pas en cause : l'arrondi de tout le calcul (exécutable, MFEM,
  nombre de processus) décide du point où le point fixe s'arrête.
- Cause : les composantes latérales de F, juste sous 1, ne sont connues
  qu'à un ulp près (6e-8). Quand les deux contraintes latérales sont
  proches d'une même valeur S, la correction du point fixe vaut environ
  (1 + ν) (1 - 2 ν) S / E (0,47 S / E pour E = 173 GPa et ν = 0,321, modules
  effectifs du cermet), et elle s'arrondit à zéro sous un demi-ulp, soit
  |res| inférieur à environ 1,5e4 Pa. Le point fixe ne peut pas descendre
  de façon sûre sous ce plancher.
- Seuil retenu en float : 3e4 Pa, environ deux fois le plancher
  (mm-opera-hpc 19d64ba, après 2e4881c qui avait pris 5e4 Pa). Tous les calculs
  passent, avec un écart au double de +0,10 à +0,16 % à 0,2 s et de
  +0,06 à +0,07 % à 0,4 s. L'essentiel de l'écart à 0,2 s vient de
  l'arrondi de la déformation imposée : 1 + 1e-4 vaut 1,000100017 en
  float, soit +0,17 % sur la déformation. Les contraintes latérales
  laissées par le seuil décalent σ33 de ν (σ11 + σ22), au plus 0,08 % avec
  3e4 Pa (0,13 % avec 5e4 Pa, d'où les +0,25 %). 2e4 Pa marche aussi ici,
  mais à 1,16 fois seulement de l'oscillation la plus haute mesurée
  (17 235 Pa).
- Le polycristal utilise la même classe et les mêmes modules : ses
  mesures en float (results/precision-polycristal-2026-10-14) sont
  identiques au bit près avec le seuil de 3e4 Pa.
- Le double ne change pas (mêmes fichiers au bit près).
