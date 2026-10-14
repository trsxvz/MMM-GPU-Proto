# Bilan sur le GPU des lots du 11/10/2026 (poutre r4)

Script : `cases/measure-gpu-summary.sh`. Poutre r4 (111 843 inconnues, 10
pas), un processus et la RTX 5080, tout par noyaux, loi `generic-parallel`
compilée pour CUDA. Avant : mfem-mgis 8a6325335e (`install/mfem-mgis`) ;
après : e111df8454 (`install/mfem-mgis-globale` : positions des entrées,
matrices et action sans former B, inverses des jacobiennes, matrice
globale gardée ; cette dernière n'est pas utilisée avec hypre sur le GPU).
Machine calme.

Boucle sur les pas (« time » de `summary.txt`) et phases
(`cases/summary-profiling.py`), en secondes :

| Variante | Avant | Après | Tangente avant / après | Résolution avant / après |
|---|---|---|---|---|
| Newton, tangente assemblée | 3,17 | 2,92 | 0,22 / 0,14 | 2,25 / 2,12 |
| NewtonK0 + Eisenstat-Walker, tangente assemblée | 1,90 | 1,71 | 0,22 / 0,15 | 1,11 / 1,04 |
| NewtonK0 + Eisenstat-Walker, sans assemblage | 2,12 | 1,87 | 0 / 0 | 1,48 / 1,26 |

- Gains de 9 à 14 % : la résolution (gradient conjugué et BoomerAMG de
  hypre sur le GPU) domine. Sans assemblage, la résolution gagne 15 % par
  le noyau de l'action (B non formé, inverses des jacobiennes gardées).
- Mêmes solutions à la tolérance de Newton près (normes identiques à
  1e-9 près, itérations de Newton à une près pour les variantes à
  tolérance adaptée).
