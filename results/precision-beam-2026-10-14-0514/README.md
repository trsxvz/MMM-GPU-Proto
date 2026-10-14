# Précision float contre double sur le CPU, poutre (14/10/2026)

Poutre r4 du prototype (111 843 inconnues), 10 pas, Newton, PCG de hypre
préconditionné par BoomerAMG (tolérance 1e-4), loi de von Mises à
écrouissage linéaire. Deux processus, un par cœur physique (cœurs 0 et 1,
avec leurs jumeaux 8 et 9), machine calme (la session voisine ne compile
pas, la compilation de MFEM en float pour le GPU est suspendue). Trois
répétitions, les précisions alternant. Script :
`cases/measure-precision-beam.sh`, synthèse : `cases/summary-precision.py`,
temps par phase : `cases/summary-profiling.py`.

Les deux constructions ne diffèrent que par la précision : environnements
spack jumeaux `deps/spack/mmm-double-cpu` et `deps/spack/mmm-float-cpu`
(MFEM 4.10, hypre 3.2.0 et MUMPS en double ou en simple, sans CUDA), MGIS
en double ou en float, même mfem-mgis (gpu-proto 283e0166e9), loi générée
par mfront dans le type de la construction.

| Calcul | Critère de Newton | Temps min. (s) | Médiane (s) | Itérations | Norme des déplacements |
|---|---|---|---|---|---|
| double | 1e-8 (celui du cas) | 34,27 | 34,37 | 28 | 48,83704586 |
| double | 1e-4 | 23,12 | 23,40 | 16 | 48,83715795 |
| float | 1e-4 | 20,74 | 20,83 | 16 | 48,83312988 |

Temps par phase du premier passage à 1e-4 (secondes, nombre d'appels) :

| Calcul | Total | Intégration | Résidu | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|---|
| double | 23,49 | 0,76 (26) | 0,49 (26) | 3,53 (16) | 11,33 (16) | 7,25 (10) |
| float | 20,93 | 0,69 (26) | 0,48 (26) | 3,42 (16) | 9,52 (16) | 6,71 (10) |

- À critère égal, le float est 1,11 fois plus rapide, avec les mêmes
  itérations de Newton. Le gain vient surtout de la résolution linéaire
  (1,19 fois) et de la prédiction (1,08 fois) ; l'intégration des lois,
  le résidu et la tangente ne changent presque pas (calculs élément par
  élément sans vectorisation qui profiterait de la simple précision).
- La norme des déplacements en float s'écarte de 8,0e-5 de celle du double
  convergé à 1e-8 (contre 2,3e-6 pour le double à 1e-4) : l'écart vient de
  la précision des calculs, accumulée sur les dix pas de plasticité.
- Le critère du cas (1e-8), borné à 100 ε = 1,2e-5 par mfem-mgis en
  float, n'est pas atteint à r4 : le résidu relatif plafonne vers 1,3e-5
  et Newton s'arrête au bout de 50 itérations (premier passage,
  `results/precision-beam-2026-10-14-0506`). À r2 (2 475 inconnues), le
  plancher était sous la borne. Le plancher croît avec la taille du
  problème : en float, la tolérance de Newton doit être choisie au-dessus
  (1e-4 ici).
