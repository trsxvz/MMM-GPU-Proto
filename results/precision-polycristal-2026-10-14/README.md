# Polycristal en float : écart au double et limites (14/10/2026)

Polycristal de mm-opera-hpc (`uniaxial-polycrystal`, maillage
`5crystals.msh`, loi `Mono_UO2_Cosh_Jaco3` en grandes transformations),
un seul pas de durée d (0,1, 0,4 ou 1), intégration point par point sur le
CPU, deux processus. Float : build/mm-opera-hpc-float-gpu (environnement
mmm-float-gpu). Double : build/mm-opera-hpc-9b. La loi est générée par le
mfront du prototype avec l'interface `generic`, en float
(`--dsl-option=real_type:float`) ou en double, et passée par `--library`.
Variantes de la loi : telle qu'écrite (`ref`), avec `cosh(x) - 1` et
`sinh(x)` calculés depuis `expm1(x)` sans cancellation (`expm1`), et avec
`@Epsilon 1e-6` au lieu de 1e-11 (`eps6`, float seulement). Fichiers `.res` :
temps, F11, F22, F33, puis σ11, σ22, σ33 (Pa) à la fin du calcul.

| d | Double σ33 (Pa) | Float σ33 (Pa), loi telle qu'écrite | Écart | Remarque |
|---|---|---|---|---|
| 0,1 | 8 997 894 | 9 053 098 | +6,1e-3 | F33 - 1 vaut 4,9949e-5 en float au lieu de 5e-5 (arrondi de F) |
| 0,4 | 35 338 761 | 35 373 936 | +1,0e-3 | |
| 1 | 60 410 661 | arrêt à 0,86 s | | échecs répétés de l'intégration de la loi (51 messages sur les deux processus), pas coupé jusqu'à 0,002 s, nombre maximal de sous-pas atteint ; même comportement dans un premier essai (arrêt à 0,85 s) |

- L'écart en float décroît quand la déformation croît (6,1e-3 à 5e-5 de
  déformation axiale, 1,0e-3 à 2e-4).
- La forme en `expm1` donne le même résultat à d = 0,1 et presque le
  même à d = 0,4 (écart 1,0e-3) : la cancellation de `cosh(x) - 1` n'est
  pas en cause. En double, elle ne change rien. À d = 1, elle échoue aussi
  (arrêt à 0,53 s).
- Avec `@Epsilon 1e-6`, la loi ne s'arrête plus, mais l'écart à d = 0,1
  devient -1,5e-2, et à d = 1 le point fixe sur la contrainte
  macroscopique stagne à 78 000 Pa, au-dessus de son seuil de 1e4 Pa
  (arrêt à 0,5 s).

Interprétation : F est stocké en float près de l'identité. Il ne bouge
que par pas d'un ulp (6e-8 sous 1, 1,2e-7 au-dessus), soit 1e-3 d'une
déformation de 5e-5, et un quantum de contrainte d'environ 1,5e4 Pa pour
un module de 2,5e11 Pa, du même ordre que la tolérance du point fixe de
ces mesures (1e4 Pa).
La déformation imposée elle-même est arrondie. Ces mesures plaident pour
une interface en H = F - I (gradient du déplacement, petit, gardant toute
la précision relative du float) pour les lois en grandes transformations
en float ; décision à prendre (voir docs/plan-precision.md).

Ces mesures ont été faites avec une tolérance de 1e4 Pa pour le point
fixe sur les contraintes latérales. Depuis mm-opera-hpc 05358bd, il prend
au moins deux fois sa résolution, soit 30 804 Pa en float ici
(results/cermet-float-2026-10-14). Refaites avec cette règle (fichiers
`*-resolution.res`) :

- d = 0,1 et 0,4 : mêmes fichiers au bit près (le point fixe s'arrête aux
  mêmes itérations) ;
- d = 1, loi telle qu'écrite : arrêt à 0,82 s, par les seuls échecs de
  l'intégration de la loi (50 messages), sans échec du point fixe (avant,
  deux échecs du point fixe s'y ajoutaient) ; forme en `expm1` : arrêt à
  0,57 s, même cause ;
- `@Epsilon 1e-6` : à d = 0,1, le pas n'est plus coupé et l'écart devient
  -1,3e-2 ; à d = 1, le point fixe stagne encore vers 69 660 Pa, au-dessus
  de la tolérance (arrêt à 0,66 s) : vraie non-convergence, due au
  critère lâche de la loi.

La conclusion ne change pas : avec un pas de 1 s, la loi échoue en float.
