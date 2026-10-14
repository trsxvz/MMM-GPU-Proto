# Temps de référence du lot L0 (09/10/2026)

Chemin actuel de mfem-mgis (branche `gpu-proto`, 373550f0d0), MFEM 4.10
de l'environnement `mmm-gpu`, sur CPU (pas de périphérique), un fil
OpenMP par processus. Machine : AMD Ryzen 7 9800X3D (8 cœurs, 16 fils),
30 Go, machine calme. Script : `cases/measure-reference.sh` ; résumé par
phase : `cases/summary-profiling.py *.txt`.

## Cas

| Cas | Inconnues | Processus | Pas | Newton (total) |
|---|---|---|---|---|
| poutre r3, PCG + BoomerAMG | 15 795 | 4 | 10 | 22 |
| poutre r3, PCG + BoomerAMG, tolérance 1e-4 | 15 795 | 4 | 10 | 25 |
| poutre r3, MUMPS | 15 795 | 4 | 10 | 22 |
| poutre r4, PCG + BoomerAMG, tolérance 1e-4 | 111 843 | 8 | 10 | 31 |
| poutre r4, MUMPS | 111 843 | 8 | 10 | 26 |
| ssna303 3D (ex4), FGMRES + BoomerAMG | 27 945 | 4 | 10 | 26 |
| ssna303 3D (ex4) r1 | environ 8 fois plus | 8 | 10 | 26 |
| polycristal, 5 grains (mm-opera-hpc) | 20 553 | 4 | 10 | 98 |
| polycristal r1 | 168 885 | 8 | 10 | 147 |
| thermomécanique (ex8), MUMPS | 9 516 | 4 | 2 | 17 |

La poutre donne la même solution avec tous les solveurs (norme finale
18,38573629 en r3, 48,83704586 en r4).

## Temps par phase (secondes, maximum sur les processus)

Entre parenthèses : nombre d'appels. Intégration : intégration des lois
(`NS::processNewUnknownsEstimate`) ; préconditionneur : construction du
solveur linéaire (factorisation de MUMPS ; BoomerAMG se construit au
premier produit, compté dans la résolution) ; prédiction : prédiction
élastique au début du pas (assemblage et résolution).

| Cas | Total | Intégration | Résidu | Tangente | Préconditionneur | Résolution | Prédiction | Fin de pas |
|---|---|---|---|---|---|---|---|---|
| poutre r3 PCG+AMG | 2,55 | 0,09 | 0,08 | 0,32 | 0,00 | 1,32 | 0,73 | 0,00 |
| poutre r3 PCG+AMG 1e-4 | 2,20 | 0,10 | 0,09 | 0,36 | 0,00 | 1,13 | 0,51 | 0,00 |
| poutre r3 MUMPS | 3,42 | 0,09 | 0,08 | 0,32 | 1,32 | 0,06 | 1,54 | 0,00 |
| poutre r4 PCG+AMG 1e-4 | 15,45 | 0,50 | 0,43 | 2,19 | 0,01 | 9,43 | 2,84 | 0,03 |
| poutre r4 MUMPS | 65,62 | 0,41 | 0,36 | 1,84 | 35,81 | 0,79 | 26,34 | 0,03 |
| ssna303 3D | 10,35 | 3,18 | 0,31 | 1,91 | 0,00 | 2,44 | 2,37 | 0,05 |
| ssna303 3D r1 | 66,86 | 12,61 | 1,26 | 8,37 | 0,01 | 28,95 | 14,87 | 0,48 |
| polycristal | 49,73 | 27,48 | 1,13 | 3,22 | 0,00 | 17,65 | - | 0,01 |
| polycristal r1 | 386,71 | 139,16 | 6,57 | 26,13 | 0,13 | 214,33 | - | 0,11 |
| thermomécanique | 5,07 | 1,62 | 0,37 | 0,85 | 1,70 | 0,46 | - | 0,01 |

## Constats

- **Ce que le GPU doit viser** : la résolution linéaire (jusqu'à 61 % du
  temps pour la poutre, 55 % pour le polycristal raffiné) et
  l'intégration des lois (55 % du polycristal, 31 % de ssna303). La
  tangente pèse 7 à 18 %. Le résidu est faible partout.
- **Prédiction** : avec MUMPS, elle refactorise l'opérateur élastique à
  chaque pas (26 s sur 66 pour la poutre r4). C'est le gain attendu de la
  prédiction avec `K0` factorisé une fois (plan, point fixe).
- **Comparaison avec l'ancienne pile** (MFEM 4.9, mfem-mgis master
  4c8a91d815, étude du 04/10), même cas, même Newton de mfem-mgis, poutre
  r4 à 8 processus, mêmes solutions et mêmes itérations :
  - PCG + BoomerAMG : 15,4 s contre 17,5 s (plus rapide) ;
  - MUMPS : 65,6 s contre 42,2 s (1,55 fois plus lent). OpenMPI compilé
    avec CUDA en explique environ 10 % : 59,3 s avec
    `mpirun --mca accelerator null`. Le reste n'est pas expliqué : MUMPS
    5.9.1 contre 5.8.1, ScaLAPACK 2.2.3 contre 2.2.2, OpenBLAS avec
    OpenMP contre sans, MFEM 4.10 contre 4.9. Les deux MUMPS sont
    construits sans METIS ni SCOTCH (renumérotation AMD ou PORD).
    À essayer lors de la reconstruction de MFEM : MUMPS `+metis
    +parmetis`.
- Le temps d'OpenMPI compilé avec CUDA sur les échanges du CPU (environ
  10 % ici) est à garder en tête pour les comparaisons CPU contre GPU :
  les mesures CPU peuvent être refaites avec `--mca accelerator null`.
