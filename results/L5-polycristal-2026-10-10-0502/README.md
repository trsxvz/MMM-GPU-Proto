# Newton inexact sur le polycristal (10/10/2026)

Polycristal de mm-opera-hpc : 5 grains, loi de monocristal de l'UO2,
problème périodique, traction uniaxiale, 10 pas. mfem-mgis 4b72095ee0.
Script : `cas/mesures-L5-polycristal.sh` (options
`--non-linear-solver` et `--adaptive-linear-solver-tolerance` du
polycristal). Une RTX 5080, 16 cœurs, machine calme (charge 1,9).

Tous les noyaux. Newton : PCG de hypre et BoomerAMG reconstruit à chaque
itération. NewtonK0 : gradient conjugué préconditionné par BoomerAMG
construit une fois sur l'opérateur élastique. Tolérance du solveur
linéaire du polycristal (1e-12), ou Eisenstat-Walker 1 avec une tolérance
initiale de 0,1 (EW1).

| Cas | CPU, 8 processus | GPU |
|---|---|---|
| Newton | 29,60 s (98 it.) | 11,88 s (98 it.) |
| Newton, EW1 | 33,96 s (156 it.) | 11,57 s (162 it.) |
| NewtonK0 | 25,02 s (98 it.) | 11,22 s (98 it.) |
| NewtonK0, EW1 | **23,46 s** (156 it.) | **9,58 s** (164 it.) |

Temps totaux (`cas/resume-profilage.py`), itérations de Newton au total
(une résolution linéaire chacune).

Détail sur le GPU :

| Cas | Intégration | Résolution |
|---|---|---|
| Newton | 5,09 s (158) | 5,11 s (98) |
| NewtonK0 | 5,16 s (159) | 4,45 s (98) |
| NewtonK0, EW1 | 7,12 s (225) | 0,70 s (164) |

## Constats

- **NewtonK0 avec Eisenstat-Walker 1 est le meilleur sur les deux
  machines** : ×1,24 contre Newton sur le GPU, ×1,26 sur 8 cœurs. La
  résolution linéaire devient négligeable sur le GPU (0,70 s contre
  5,11 s).
- **La loi domine maintenant** : Eisenstat-Walker demande 66 itérations de
  Newton de plus, donc autant d'intégrations de la loi (7,1 s sur 9,6 s
  sur le GPU). La loi en simple précision ou un meilleur noyau de loi sont
  les leviers suivants sur ce cas.
- Avec Newton classique, Eisenstat-Walker ne paie pas sur CPU (chaque
  itération reconstruit BoomerAMG), à peine sur le GPU.
- Contraintes macroscopiques : écarts à Newton de 9,6e-7 à 3,5e-6 avec
  Eisenstat-Walker (2,4e-13 sans). Le Newton microscopique du polycristal
  converge à 1e-4 près en relatif : avec une tolérance linéaire de 1e-12,
  la convergence quadratique dépasse largement ce critère ; avec
  Eisenstat-Walker, la solution s'arrête plus près du critère. Les écarts
  sont dans la tolérance du calcul.
