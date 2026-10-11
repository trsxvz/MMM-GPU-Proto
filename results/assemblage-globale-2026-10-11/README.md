# Assemblage local et matrice globale (11/10/2026)

Banc `cas/assemblage` en parallèle (`--parallel`) : assemblage par noyaux
des matrices de rigidité dans la matrice locale de chaque processus, puis
construction de la matrice globale comme mfem-mgis
(`assembleTangentOperatorByKernels`) : matrice diagonale par blocs,
produit `Pᵀ A P` de hypre, élimination des degrés de liberté d'une
extrémité encastrée. mfem-mgis 4a007e85f2 (installation séparée
`install/mfem-mgis-registres`). Temps médians de 10 assemblages, maxima
sur les processus. Sorties : `cpu-np8.txt`, `gpu.txt`.

| Cas | Assemblage local | Pᵀ A P | Élimination |
|---|---|---|---|
| 8 processus, poutre r4 | 9,5 ms | 15,4 ms | 0,3 ms |
| 8 processus, thermique r5 | 14,8 ms | 36,5 ms | 0,2 ms |
| 8 processus, grandes transformations, hexaèdres Q1 r4 | 6,5 ms | 3,1 ms | 0,2 ms |
| GPU, 1 processus, poutre r4 | 2,3 ms | 1,1 ms | 0,2 ms |
| GPU, 1 processus, thermique r5 | 13,8 ms | 1,5 ms | 0,04 ms |
| GPU, 1 processus, grandes transformations r4 | 4,5 ms | 0,4 ms | 0,5 ms |

- Sur 8 processus, le produit `Pᵀ A P` coûte maintenant plus que
  l'assemblage local pour les éléments linéaires : hypre recalcule le
  motif et les communications à chaque tangente, alors qu'ils ne changent
  pas. Piste : motif de la matrice globale et schéma d'échange gardés, seules
  les valeurs étant envoyées (point du plan L4).
- Sur le GPU à un processus, ce produit est faible.
- Nsight Compute (un lancement du noyau des matrices) : thermique et
  mécanique limitées par le calcul (85 à 89 % du débit des
  multiprocesseurs, mémoire presque inactive), d'où le lot suivant
  (inverses des jacobiennes calculées une fois,
  `results/jacobiennes-inverses-2026-10-11`).
