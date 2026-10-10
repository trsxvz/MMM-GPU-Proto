# NewtonK0 sur le polycristal (10/10/2026)

Polycristal de mm-opera-hpc : 5 grains, loi de monocristal de l'UO2
(orthotrope, grandes transformations, implicite), problème périodique,
traction uniaxiale, 10 pas. mfem-mgis 9249013150. Script :
`cas/mesures-L5-polycristal.sh`. Une RTX 5080, 16 cœurs. Machine calme au
départ (charge 1,09).

Tous les noyaux (intégration par lots, gradients, forces internes et
tangente assemblée par noyaux). Sur CPU : 8 processus, loi `generic`. Sur
GPU : un processus, loi `generic-parallel` compilée pour `cuda`.

- **Newton** : PCG de hypre préconditionné par BoomerAMG, tolérance
  1e-12, nouvelle mise en place de BoomerAMG à chaque itération de Newton
  (sur la tangente cohérente).
- **NewtonK0** : option `--non-linear-solver NewtonK0` du polycristal,
  méthode `NewtonK0` du solveur `FixedPoint` de mfem-mgis. Gradient
  conjugué de MFEM (tolérance 1e-12) préconditionné par un cycle de
  BoomerAMG construit une seule fois, sur l'opérateur élastique K0 obtenu
  par une intégration de la loi (`ElasticOperatorComputation` à
  `Integration`, `ElasticOperatorUpdate` à `Once`).

| Cas | Total | Intégration | Tangente | Résolution |
|---|---|---|---|---|
| CPU, 8 processus, Newton | 29,60 | 13,84 (158) | 2,02 (98) | 13,42 (98) |
| CPU, 8 processus, NewtonK0 | **24,92** | 14,08 (159) | 2,04 (99) | 8,53 (98) |
| GPU, 1 processus, Newton | 11,88 | 5,09 (158) | 0,20 (98) | 5,11 (98) |
| GPU, 1 processus, NewtonK0 | **11,26** | 5,16 (159) | 0,20 (99) | 4,45 (98) |

Temps en secondes (maximum sur les processus), nombre d'appels entre
parenthèses. Le résidu (0,07 s sur CPU, 0,02 s sur GPU) et la fin des pas
sont négligeables. Contraintes macroscopiques de NewtonK0 égales à celles
de Newton dans la même configuration, à 2,4e-13 (CPU) et 1,8e-14 (GPU)
près de l'échelle des contraintes.

## Constats

- **Mêmes itérations de Newton** : 98 résolutions dans les deux cas.
  NewtonK0 fait une intégration et un assemblage de plus (K0, au premier
  pas).
- **CPU** : la résolution passe de 13,42 s à 8,53 s (-36 %), le calcul de
  29,60 s à 24,92 s (×1,19). La mise en place de BoomerAMG à chaque
  itération pèse donc environ un tiers de la résolution de Newton.
- **GPU** : la résolution passe de 5,11 s à 4,45 s (-13 %), le calcul de
  11,88 s à 11,26 s (×1,06). La mise en place de BoomerAMG est moins
  chère sur le GPU ; l'intégration de la loi (5,1 s, 32 ms par appel)
  domine désormais.
- Premier essai sur GPU (même jour, machine moins calme) : 11,90 s contre
  11,24 s, écart 1,3e-14 : reproductible.
- Comme sur la poutre (`results/L5-poutre-2026-10-10-0250/README.md`),
  NewtonK0 est le seul gain stable sur GPU parmi les variantes du point
  fixe ; il ne touche pas au coût de la loi.
