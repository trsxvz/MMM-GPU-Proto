# Barreau entaillé ssna303 (13/10/2026)

Campagne de mesures, cas 2, reprise de la série du 11/10
(`results/ssna303-ssna303-2026-10-11-0625`, arrêtée par manque de
mémoire). Exemple ex4 de mfem-mgis-examples, variante résolue par hypre
(Newton, FGMRES de hypre préconditionné par BoomerAMG), plasticité en
grandes transformations (déformation logarithmique de Hencky), chargement
complet en 10 pas, maillage raffiné zéro, une et deux fois. Script :
`cases/measure-ssna303.sh`, avec `examples_build=build/examples-f57072220f`
(ex4 construit contre une installation séparée de mfem-mgis f57072220f,
`install/mfem-mgis-f57072220f`). Une RTX 5080 (16 Go), 8 cœurs (16 fils),
30 Go de mémoire ; une autre session de travail sur la machine, sans
calcul pendant la série.

| Raffinement | Inconnues | 8 cœurs, point par point | 8 cœurs, noyaux | GPU, noyaux | Gain du GPU |
|---|---|---|---|---|---|
| 0 | 27 945 | 7,13 s | 5,18 s | 6,89 s | ×0,75 |
| 1 | 196 047 | 71,0 s | 59,3 s | 26,0 s | ×2,3 |
| 2 | 1 461 915 | 693 s | 601 s | ne tient pas | |

FBar (`--use-fbar`) :

| Raffinement | 8 cœurs, point par point | 8 cœurs, noyaux | GPU, noyaux | Gain du GPU |
|---|---|---|---|---|
| 0 | 7,64 s | 5,86 s | 7,28 s | ×0,80 |
| 1 | 90,9 s | 74,0 s | 27,9 s | ×2,7 |
| 2 | non lancé | non lancé | ne tient pas | |

Temps de la boucle sur les pas ; le gain compare le GPU aux noyaux sur 8
cœurs.

- Par rapport au 11/10 (mfem-mgis 8a6325335e), les noyaux sur 8 cœurs
  gagnent 10 % au raffinement 1 (65,5 → 59,3 s ; FBar 78,7 → 74,0 s) et
  au raffinement 2 (670 → 601 s) ; le GPU gagne 4 à 5 % (27,1 → 26,0 s ;
  FBar 29,5 → 27,9 s). Ces gains viennent des lots du 11/10 (positions
  des entrées, matrices sans former B, inverses des jacobiennes).
- Raffinement 2 sur 8 processus : la mémoire de la machine est presque
  pleine pendant le calcul par noyaux (27,5 Go sur 30, 4,5 Go d'échange) ;
  FBar n'y est pas lancé (il avait épuisé la mémoire le 11/10).
- Raffinement 2 sur le GPU : avec la règle à 27 points des hexaèdres, la
  mémoire du GPU ne suffit pas (11/10). La règle d'ordre 3 (8 points,
  option `-iro 3`) a été essayée : elle ne convient pas à ce cas. Sur
  8 cœurs, Newton échoue 8 fois et les pas concernés sont redécoupés
  (810 s contre 601 s) ; sur le GPU, Newton stagne à 1,1e-4 de
  résidu relatif au premier pas, puis le redécoupage échoue faute de
  mémoire (allocation de 451 Mo de hypre dans le triple produit). Ce
  raffinement ne tient pas sur ce GPU avec la règle complète ; il faudrait
  réduire la mémoire aux points (blocs de la tangente en simple précision
  ou calculés à la volée) ou utiliser plusieurs GPU.
- Au raffinement 0 (28 000 inconnues), le GPU reste plus lent que
  8 cœurs : le système est trop petit pour lui.
