# Mesures du lot L2 sur le polycristal (09/10/2026)

Polycristal de mm-opera-hpc : 5 grains, loi de monocristal de l'UO2
(grandes transformations, orthotrope, schéma implicite), 20 553
inconnues, 10 pas, PCG + BoomerAMG. mfem-mgis 4a68166ca9 (lot L2 relu),
loi `generic-parallel` compilée avec les options corrigées (patch 0093 de
TFEL). Machine calme (charge 0,9 au départ). Script :
`cases/measure-L2.sh polycristal` ; résumé par phase :
`cases/summary-profiling.py *.txt`.

Intégration : `NS::processNewUnknownsEstimate` (gradients élément par
élément puis appel de la loi). Temps en secondes, maximum sur les
processus.

| Cas | Total | Intégration | Résidu | Tangente | Résolution |
|---|---|---|---|---|---|
| 1 processus, point par point | 150,35 | 104,90 | 4,21 | 8,55 | 32,37 |
| 1 processus, lots, `generic` | 148,31 | 103,07 | 4,21 | 8,48 | 32,23 |
| 1 processus, lots, `openmp` 1 fil | 150,72 | 105,56 | 3,97 | 8,47 | 32,39 |
| 1 processus, lots, `openmp` 8 fils | 63,03 | 17,84 | 4,03 | 8,39 | 32,45 |
| 1 processus, lots, `stlpar` 8 cœurs | 59,13 | 13,94 | 4,08 | 8,33 | 32,47 |
| 2 processus, point par point | 85,64 | 53,35 | 2,27 | 5,65 | 24,09 |
| 2 processus, lots, `openmp` 4 fils | 49,44 | 17,23 | 2,19 | 5,53 | 24,22 |
| 8 processus, point par point | 30,52 | 14,51 | 0,60 | 1,72 | 13,45 |
| 8 processus, lots, `generic` | 30,38 | 14,33 | 0,62 | 1,69 | 13,50 |
| 8 processus, lots, `openmp` 1 fil | 30,73 | 14,72 | 0,58 | 1,70 | 13,49 |

Contraintes macroscopiques comparées au calcul point par point sur le
même nombre de processus : identiques avec la loi `generic` par lots ;
écart de 4e-9 au plus (rapporté au maximum de chaque colonne) avec la loi
`generic-parallel`, compilée autrement, au niveau de la tolérance du
solveur itératif.

## Constats

- La loi pèse 70 % du calcul sur un processus. Par lots, plusieurs fils
  la divisent par 5,9 (`openmp`, 8 fils) et par 7,5 (`stlpar`, TBB sur les
  16 fils matériels des 8 cœurs) ; le calcul complet passe de 150 s à
  59 s.
- Un fil par processus : la loi `generic-parallel` va aussi vite que la
  loi `generic`, une fois les options de compilation corrigées (elle
  était 17 % plus lente avant le patch 0093 de TFEL).
- À 8 processus d'un fil, la loi prend 14,5 s, comme 8 fils dans un
  processus : sur CPU, le gain des fils vient quand le nombre de processus
  est limité (mémoire, GPU par processus) ; le reste du calcul (résolution
  surtout) passe alors en tête.
- La phase d'intégration comprend aussi le calcul des gradients élément
  par élément, séquentiel : c'est la partie que le lot L3 remplace par des
  noyaux.

## Premier essai (abandonné, gardé pour mémoire)

Les premières mesures (dossier `scratch/`) donnaient `stlpar` séquentiel
et `openmp` à 8 fils seulement 4,3 fois plus rapide : `scripts/env.sh`
exporte `OMP_PROC_BIND=spread OMP_PLACES=threads`, et libgomp attachait le
fil principal au premier processeur (TBB ne voyait qu'un processeur, et
`spread` n'occupait que 4 cœurs). Le script de mesures retire maintenant
ces variables.
