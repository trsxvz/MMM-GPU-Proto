# Cermet en float : tolérance du point fixe sur les contraintes latérales (14/10/2026)

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

Fichiers `<exécutable>-<n>p-<loi>-<tolérance>.res` (n processus) et
`double-*.res` : temps, F11, F22, F33, puis σ11, σ22, σ33 (Pa) à la fin de
chaque pas. `fixed-point.txt` : résidu du point fixe,
|res| = (σ11² + σ22²)^½, à chaque itération de chaque calcul.

## Règle retenue

La tolérance sur les contraintes latérales vaut 1e4 Pa, et le point fixe
prend au moins deux fois sa résolution (mm-opera-hpc 05358bd et
ae5768b). Quand les deux contraintes latérales sont proches d'une même
valeur S, la correction des composantes latérales de F vaut environ
(1 + ν) (1 - 2 ν) S / E. Elle s'arrondit à zéro quand elle est plus petite
que la moitié de l'écart entre deux flottants voisins de ces composantes.
Avec les modules effectifs de Kröner des grains d'UO2 (E = 173 GPa,
ν = 0,321), utilisés par le cermet et le polycristal, et des composantes
latérales juste sous 1 (écart de 6e-8), la résolution vaut 1,54e4 Pa, et
la tolérance 30 804 Pa. En compression, les composantes latérales
dépassent 1 et la résolution double. En double, elle vaut environ
3e-5 Pa : rien ne change (mêmes fichiers au bit près).

Les fichiers `*-resolution.res` (règle retenue) sont identiques au bit
près aux fichiers `*-3e4.res` de la même configuration.

## Balayage des tolérances constantes

Double, exécutable de la CI, un processus (le test) : σ33 = 17 219 167 Pa
à 0,2 s et 29 377 289 Pa à 0,4 s par pas de 0,2 s ; 17 642 604 Pa
(+2,46 %) et 29 886 489 Pa (+1,73 %) par pas de 0,1 s.

Float : σ33 à 0,2 s puis à 0,4 s, et écart au double par pas de 0,2 s.

| Exécutable | Processus | Tolérance (Pa) | Loi compilée par gcc 13 | Loi compilée par gcc 14 |
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

- Avec 1e4 Pa, le point fixe du premier pas ne converge pas dans cinq
  calculs sur six : il stagne vers 10 684 à 10 830 Pa avec l'exécutable de
  la CI, et oscille entre 17 234,6 et 15 454,4 Pa avec l'exécutable du
  prototype et la loi de gcc 14 ; avec ce même exécutable, la loi de
  gcc 13 converge. L'arrondi de tout le calcul, loi comprise, décide du
  point où le point fixe s'arrête. Le pas est rejeté puis coupé en deux
  pas de 0,1 s : l'écart de 2,7 % au double vient de ce pas de temps (par
  pas de 0,1 s, le float s'écarte du double de 0,24 à 0,27 %).
- 2e4 Pa marche ici, mais à 1,16 fois seulement de l'oscillation la plus
  haute mesurée (17 235 Pa) ; 5e4 Pa laisse des contraintes latérales
  plus grandes (+0,25 % à 0,2 s).

## Écart au double avec la règle retenue

+0,10 à +0,16 % à 0,2 s et +0,06 à +0,07 % à 0,4 s. À 0,2 s, sur les
calculs de la CI (+0,16 %), à peu près :

- +0,07 % : tolérance absolue de l'équilibre microscopique en float
  (6e4), mesurée en double avec cette tolérance
  (`double-microscopic-tolerance-6e4.res` : +0,072 % à 0,2 s, +0,023 %
  à 0,4 s) ;
- +0,03 % : contraintes latérales laissées par le point fixe, qui
  décalent σ33 de ν (σ11 + σ22) (-0,045 % sur les calculs du prototype) ;
- +0,017 % : arrondi de la déformation imposée (1 + 1e-4 vaut
  1,000100017 en float) ;
- le reste, environ +0,04 %, est l'arrondi du calcul en float.

## Polycristal

Il utilise la même classe et les mêmes modules. Ses mesures en float
(results/precision-polycrystal-2026-10-14), refaites avec la règle
retenue : mêmes fichiers au bit près pour des pas de 0,1 et 0,4 s ; pour
les autres, voir le README de ce dossier.
