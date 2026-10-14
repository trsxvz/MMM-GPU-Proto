# Inverses des jacobiennes : bout en bout sur le GPU (11/10/2026)

Script : `cases/measure-jacobians-gpu.sh`. Poutre r4 (10 pas) et
thermique r5 (5 pas), un processus et la RTX 5080, tout par noyaux, loi
`generic-parallel` compilée pour CUDA. Avant : mfem-mgis 4a007e85f2
(`install/mfem-mgis-registres`) ; après : bf7673f4f8
(`install/mfem-mgis-jacobiennes`, inverses des jacobiennes calculées une
fois ; la version fusionnée, 7af829b3a5, ne change que la documentation,
les noms de deux fonctions et un test). Machine calme.

Temps en secondes (`cases/summary-profiling.py`, total avec la mise en place ;
« time » de `summary.txt` : boucle sur les pas) :

| Cas | Boucle sur les pas | Intégration | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|
| poutre, avant | 3,18 | 0,05 | 0,18 | 2,33 | 0,56 |
| poutre, après | 3,04 | 0,05 | 0,14 | 2,23 | 0,55 |
| thermique, avant | 1,32 | 0,28 | 0,39 | 0,56 | |
| thermique, après | 1,08 | 0,26 | 0,18 | 0,56 | |

- Thermique : tangente ×2,2, boucle sur les pas ×1,22.
- Poutre : tangente ×1,3 ; la résolution (gradient conjugué et BoomerAMG
  de hypre) domine, gain total 4 %.
