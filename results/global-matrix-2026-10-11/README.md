# Matrice globale avec motif et échanges gardés (11/10/2026)

Banc `cases/assembly --parallel` sur 8 processus : assemblage local par
noyaux (avec les positions), puis matrice globale par le produit `Pᵀ A P`
de hypre (comme avant) et par `batch::GlobalStiffnessMatrix` (mfem-mgis
branche `globale`, 28f2367b39, installation séparée
`install/mfem-mgis-globale`) : motif de la matrice globale et échanges
entre processus calculés au premier assemblage, puis seules les valeurs
copiées, échangées (`MPI_Alltoallv`) et additionnées. Temps médians de 10
assemblages, maxima sur les processus. Machine calme. Sortie :
`cpu-np8.txt`.

| Cas | Assemblage local | Pᵀ A P de hypre | Motif gardé (premier assemblage) | Gain |
|---|---|---|---|---|
| poutre, tétraèdres P1 (r4) | 9,0 ms | 15,4 ms | 2,5 ms (32 ms) | ×6,3 |
| thermique, tétraèdres P1 (r5) | 14,1 ms | 36,8 ms | 1,7 ms (46 ms) | ×21,6 |
| grandes transformations, hexaèdres Q1 (r4) | 6,4 ms | 3,0 ms | 0,15 ms (7 ms) | ×20 |
| tétraèdres P2, mécanique (r3) | 15,2 ms | 22,0 ms | 5,8 ms (65 ms) | ×3,8 |

Mêmes matrices (écart relatif des produits par un vecteur au plus
1,4e-16). Le premier assemblage coûte deux à trois produits de hypre :
rentabilisé dès la deuxième tangente. La matrice globale ne coûte plus
qu'une fraction de l'assemblage local.
