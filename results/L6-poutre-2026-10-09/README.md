# Poutre sur GPU, tranche minimale (09/10/2026)

Poutre r4 (111 843 inconnues, tétraèdres d'ordre 1), plasticité de von
Mises, PCG + BoomerAMG (tolérance 1e-4), 10 pas. mfem-mgis 7159b6de4c,
MGIS a0b5140. Script : `cas/mesures-L6.sh poutre`. Une RTX 5080, 16 cœurs.

Calculs sur GPU : périphérique `cuda` de MFEM, hypre sur le GPU, tous les
noyaux (gradients, forces intérieures, tangente assemblée directement dans
la matrice creuse, prédiction), loi `generic` (intégrée sur l'hôte) ou
`generic-parallel` compilée pour `cuda` (intégrée sur le GPU, sur les
copies des tableaux du matériau).

Réserve : la machine sortait d'une compilation (charge moyenne de 10 au
départ, en baisse) ; le calcul de référence point par point sur 8
processus prend 17,5 s contre 15,7 s sur machine calme.

| Cas | Total | Intégration | Résidu | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|---|
| CPU, 8 processus, point par point | 17,54 | 0,58 (41) | 0,48 (41) | 2,38 (31) | 10,89 | 3,12 (10) |
| CPU, 8 processus, tout par noyaux | 14,68 | 0,20 (39) | 0,11 (39) | 1,91 (29) | 9,50 | 2,89 (10) |
| GPU, 1 processus, loi `generic` | 5,01 | 0,79 (42) | 0,05 (42) | 0,33 (32) | 2,55 | 0,99 (10) |
| GPU, 1 processus, loi `cuda` | **4,06** | 0,10 (42) | 0,03 (42) | 0,24 (32) | 2,47 | 0,87 (10) |
| GPU, 2 processus, loi `cuda` | 28,64 | 0,11 (41) | 0,04 (41) | 0,44 (31) | 22,90 | 4,91 (10) |

Temps en secondes (maximum sur les processus) ; entre parenthèses, nombre
d'appels. Même solution partout (norme 48,83704586) ; les nombres
d'itérations de Newton varient de 29 à 32 avec l'arrondi (ordre des
sommes, additions atomiques sur GPU).

## Constats

- **Un GPU contre 8 cœurs** : 4,06 s contre 17,54 s (×4,3), et contre
  14,68 s pour les mêmes noyaux sur CPU (×3,6).
- **Ce qui reste** : la résolution linéaire (PCG + BoomerAMG de hypre sur
  le GPU) fait 61 % du temps, la prédiction 21 % (dont sa propre
  résolution linéaire). La tangente coûte 7,5 ms par appel, l'intégration
  de la loi 2,4 ms par appel.
- **Loi sur le GPU** : l'intégration passe de 0,79 s à 0,10 s (plus de
  copies des tableaux entre l'hôte et le GPU).
- **Deux processus sur un même GPU** : la résolution devient 9 fois plus
  lente (22,9 s). Non analysé. Hypothèse la plus probable : sans le
  service MPS de CUDA, les contextes des deux processus se partagent le
  GPU par tranches de temps, et les nombreux petits noyaux de BoomerAMG
  s'exécutent l'un après l'autre avec des changements de contexte. À
  vérifier avec MPS (`nvidia-cuda-mps-control -d`). Un processus par GPU
  reste la configuration naturelle.
- **Mémoire** : chaque processus sur le GPU réserve une réserve de 4 Gio
  (Umpire) ; quatre processus sur le GPU de 16 Go échouent (vu en lançant
  les tests en parallèle).

## Profil Nsight Systems (3 pas, 1 processus, loi `cuda`, mfem-mgis f0239ab807)

`nsys profile -t cuda,nvtx` (Nsight Systems 2025.5 du HPC SDK 25.11 ; le
`nsys` de CUDA 13.0 ne relit pas ce rapport) sur la poutre r4, 3 pas (10
itérations de Newton, 13 résidus), environ 1,8 s :

| Copies | Nombre | Volume | Temps |
|---|---|---|---|
| hôte vers GPU | 835 | 146 Mo | 7,3 ms |
| GPU vers hôte | 5 814 | 9 Mo | 2,2 ms |
| dans le GPU | 14 704 | 3,0 Go | 13,7 ms |

- **Pas de copie cachée notable** : les grosses copies vers le GPU ne
  sont faites qu'une ou quelques fois (géométrie, motif, maillage) ; par
  itération, il reste un ou deux vecteurs (0,9 Mo) dans chaque sens (norme
  et ddl imposés de la prédiction, conditions de Dirichlet) et des
  scalaires (normes, statuts de hypre : 3 215 copies de 8 octets).
- **Temps de calcul du GPU : 0,71 s sur 1,8 s**, en 95 339 lancements de
  noyaux. Les produits matrice-vecteur creux de hypre (`csrmv`, PCG et
  BoomerAMG) en prennent 65 % (460 ms), le produit de matrices de la mise
  en place de l'AMG 78 ms, les noyaux de MFEM et de mfem-mgis
  (`CuKernel1D`, `CuKernel2D`) 69 ms, la loi (`apply`) 16 ms.
- Le reste du temps (environ 1,1 s) est passé sur l'hôte à lancer et
  attendre ces petits noyaux : la résolution est limitée par la latence
  des lancements de hypre, pas par le calcul. Pistes : moins d'itérations
  de PCG (préconditionneur réutilisé d'une itération à l'autre, options
  de BoomerAMG pour le GPU), problèmes plus gros, opérateur sans matrice.

## Tangente sans matrice (mfem-mgis 4d465460c3)

Option `TangentOperatorComputation=MatrixFree` : la tangente est appliquée
par noyaux sans assemblage (`Σ Bᵀ K B` aux points), sa diagonale aussi ;
gradient conjugué de MFEM avec un préconditionneur de Jacobi (`-ls
cg-jacobi`), même tolérance relative (1e-4) que PCG + BoomerAMG. Loi
`cuda`, 10 pas, 1 processus :

| Cas | Total | Tangente | Résolution | Prédiction |
|---|---|---|---|---|
| tangente assemblée, PCG + BoomerAMG | 3,61 | 0,22 (32) | 2,23 | 0,77 |
| tangente sans matrice, CG + Jacobi | 95,31 | 0,00 (30) | 93,30 | 1,58 |

Même solution (48,83704586), mêmes nombres d'itérations de Newton à
l'arrondi près. Le chemin sans matrice est juste, mais Jacobi ne suffit
pas : 3,1 s de gradient conjugué par itération de Newton (des milliers
d'itérations sur cette poutre élancée). Il faut un vrai préconditionneur
pour l'opérateur sans matrice : AMG construit sur la tangente assemblée
(ou sur l'opérateur élastique) et gardé plusieurs itérations, ou
multigrille polynomial avec des lisseurs de Chebyshev (lots L4 et L5 du
plan).
