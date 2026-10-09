# Tangente sans matrice sur la poutre (10/10/2026)

Poutre r4 (111 843 inconnues, tétraèdres d'ordre 1), plasticité de von
Mises, 10 pas, prédiction élastique. mfem-mgis 1bbcb877d8. Script :
`cas/mesures-L4.sh`. Une RTX 5080, 16 cœurs, machine calme.

Deux façons d'appliquer la tangente, mêmes noyaux pour les gradients, le
résidu et la loi (loi `generic` sur CPU, `cuda` sur le GPU) :

- **assemblée** : matrice creuse assemblée par noyaux
  (`TangentOperatorComputation=Kernels`), PCG + BoomerAMG de hypre
  (`-ls pcg-amg`) ;
- **sans matrice** : action `Σ Bᵀ K B` par noyaux
  (`TangentOperatorComputation=MatrixFree`), gradient conjugué de MFEM,
  préconditionneur BoomerAMG construit sur la tangente assemblée, refait
  toutes les 1, 2 ou 4 itérations de Newton (`-ls cg-amg -pu n`,
  préconditionneur `AssembledOperator`).

Même tolérance relative du solveur linéaire (1e-4) partout.

| Cas | Total | Tangente | Préconditionneur | Résolution | Prédiction |
|---|---|---|---|---|---|
| CPU 8 processus, assemblée | 19,15 | 2,60 (29) | 0,01 | 12,35 | 3,75 |
| CPU 8 processus, sans matrice, pu 1 | 27,52 | 0,00 | 2,63 | 20,79 | 3,65 |
| CPU 8 processus, sans matrice, pu 2 | 20,84 | 0,00 | 1,26 | 16,19 | 3,00 |
| CPU 8 processus, sans matrice, pu 4 | 25,28 | 0,00 | 0,75 | 21,00 | 3,13 |
| GPU, assemblée | **4,08** | 0,24 (32) | 0,01 | 2,40 | 0,94 |
| GPU, sans matrice, pu 1 | 7,38 | 0,00 | 0,23 | 5,69 | 0,98 |
| GPU, sans matrice, pu 2 | 7,27 | 0,00 | 0,18 | 5,54 | 1,05 |
| GPU, sans matrice, pu 4 | 6,83 | 0,00 | 0,11 | 5,39 | 0,85 |

Temps en secondes (maximum sur les processus), nombre d'appels entre
parenthèses. « Préconditionneur » : assemblage de la tangente pour
BoomerAMG (la mise en place de BoomerAMG est comptée dans la résolution).
Même solution partout (norme 48,83704586), 29 à 32 itérations de Newton.

## Constats

- **La tangente sans matrice est juste, mais plus lente sur ce cas** :
  6,8 s contre 4,1 s sur le GPU, 20,8 s contre 19,2 s sur 8 cœurs au
  mieux. Le gain sur l'assemblage (0,24 s) ne compense pas le coût de
  l'action dans le gradient conjugué.
- **Garder BoomerAMG plusieurs itérations** ne sert guère : la résolution
  passe de 5,7 s à 5,4 s sur le GPU (pu 1 à 4) ; sur CPU, pu 4 est plus
  lent que pu 2 (le gradient conjugué fait plus d'itérations avec un
  préconditionneur ancien).

## Profil Nsight Systems (GPU, pu 1, 10 pas)

`nsys profile -t cuda,nvtx`, temps de calcul du GPU :

| Variante | Lancements | Temps GPU | Temps total |
|---|---|---|---|
| assemblée | 318 768 | 2,62 s | 4,0 s |
| sans matrice | 312 315 | 5,15 s | 6,4 s |

- La variante assemblée est limitée par la latence des lancements (2,6 s
  de calcul sur 4 s), la variante sans matrice par le calcul.
- Une action de la tangente sans matrice coûte environ 1,3 ms sur le GPU
  (1 788 actions), contre 0,1 à 0,2 ms pour un produit par la matrice
  assemblée. Détail par action : dérivées des inconnues aux points par
  `QuadratureInterpolator::Derivatives` de MFEM, 0,89 ms (version
  générique pour les simplexes, sans produit tensoriel : c'est le poste
  principal) ; forces intérieures, 0,22 ms ; gradients à partir des
  dérivées, 0,07 ms ; produit par les blocs de la tangente, 0,04 ms ;
  restriction aux éléments, 0,05 ms.
- Pistes : un noyau fusionné par élément (lecture des inconnues de
  l'élément, gradient aux points avec les gradients des fonctions de forme
  et l'inverse de la jacobienne, produit par K, `Bᵀ`, accumulation), qui
  supprime les tableaux intermédiaires et le passage par l'interpolateur
  de MFEM ; le même noyau sert au résidu. Pour des tétraèdres d'ordre 1,
  le volume de données d'une action fusionnée (blocs de K aux points)
  reste du même ordre que celui de la matrice assemblée : l'intérêt de la
  tangente sans matrice viendra surtout de l'ordre 2 et des gros
  maillages (mémoire de la matrice).
