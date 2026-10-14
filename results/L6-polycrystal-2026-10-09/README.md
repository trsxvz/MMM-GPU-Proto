# Polycristal sur GPU, tranche minimale (09/10/2026)

Polycristal de mm-opera-hpc : 5 grains, loi de monocristal de l'UO2
(orthotrope, grandes transformations, implicite), problème périodique,
traction uniaxiale, 10 pas. PCG + BoomerAMG de hypre. mfem-mgis
decb6fd2c1, MGIS a0b5140. Script : `cases/measure-L6.sh polycristal`. Une
RTX 5080, 16 cœurs. Machine calme (charge en baisse au départ, sortie
d'une relecture).

Calculs sur GPU : périphérique `cuda` de MFEM, hypre sur le GPU, tous les
noyaux (le polycristal n'utilise pas de prédiction), loi `generic`
(intégrée sur l'hôte) ou `generic-parallel` compilée pour `cuda` (intégrée
sur le GPU, rotations des grains comprises).

| Cas | Total | Intégration | Résidu | Tangente | Résolution |
|---|---|---|---|---|---|
| CPU, 8 processus, point par point | 31,65 | 14,86 (158) | 0,65 (158) | 1,78 (98) | 14,11 |
| CPU, 8 processus, tout par noyaux | 30,55 | 14,03 (158) | 0,10 (158) | 2,09 (98) | 14,08 |
| CPU, 1 processus, point par point | 151,34 | 105,62 (158) | 4,20 (158) | 8,70 (98) | 32,49 |
| GPU, 1 processus, loi `generic` | 106,70 | 99,16 (158) | 0,15 (158) | 0,35 (98) | 5,56 |
| GPU, 1 processus, loi `cuda` | **13,18** | 5,71 (158) | 0,07 (158) | 0,22 (98) | 5,60 |
| GPU, 2 processus, loi `cuda` | 91,57 | 12,18 (158) | 0,17 (158) | 0,95 (98) | 76,99 |

Temps en secondes (maximum sur les processus), nombre d'appels entre
parenthèses. Contraintes macroscopiques de chaque calcul égales à celles
du calcul point par point sur le même nombre de processus, à 2,4e-13 à
4,4e-13 près de l'échelle des contraintes (`compare-outputs.py
--groupes 0 1-3 4-6`).

## Constats

- **Un GPU contre 8 cœurs** : 13,18 s contre 31,65 s (×2,4) ; contre un
  cœur, ×11,5.
- **La loi domine sur CPU** (105 s sur 151 s à un processus) : sur le GPU,
  l'intégration passe à 5,71 s (36 ms par appel, contre 14,9 s sur 8
  cœurs), à égalité avec la résolution linéaire (5,60 s).
- **Loi sur l'hôte avec le GPU** : 106,7 s ; la loi de monocristal sur un
  seul cœur annule tout le gain des noyaux.
- **Deux processus sur un même GPU** : résolution 14 fois plus lente
  (77 s), comme pour la poutre ; voir `results/L6-beam-2026-10-09/`.
