# cuDSS sur la poutre r4 (12/10/2026)

Poutre r4 du prototype (111 843 inconnues, von Mises à écrouissage
linéaire, 10 pas), tout par noyaux, loi cuda, un processus sur le GPU
(cœurs 0 à 3), hypre sur le GPU, prédiction élastique. mfem-mgis
45afaa0653 (`install/mfem-mgis-mumps`), construction `build/cas-cudss`
(seul l'exécutable `poutre`, lié à cette installation). Options ajoutées au
cas `poutre` : `-ls cudss` (solveur `CuDSSSolver`, symétrique défini
positif), `-mra` (réutilisation de l'analyse des solveurs directs) et
`-els cudss` (opérateur élastique du point fixe factorisé par cuDSS). Deux
passages par variante (`summary.txt`, commandes dans `measure.sh` et
`measure-k0.sh`). Temps de la boucle sur les pas, en secondes :

| Méthode | Solveur | Temps |
|---|---|---|
| Newton | PCG de hypre et BoomerAMG (tolérance 1e-4) | 3,12 et 3,04 |
| Newton | cuDSS, analyse réutilisée | 9,60 et 9,56 |
| NewtonK0 (K0 assemblé une fois, prédiction par K0) | un cycle de BoomerAMG sur K0, tolérance d'Eisenstat et Walker | 1,87 et 1,85 |
| NewtonK0 | K0 factorisé une fois par cuDSS, tolérance 1e-4 | 3,01 et 3,03 |
| NewtonK0 | K0 factorisé une fois par cuDSS, tolérance d'Eisenstat et Walker | 2,73 et 2,97 |

- À cette taille, le solveur direct ne paie pas. Avec Newton, une
  factorisation de cuDSS prend 0,235 s même avec l'analyse réutilisée
  (6,1 s pour 26 factorisations, et 3,0 s pour les 10 prédictions, qui
  factorisent l'opérateur élastique), contre 68 ms pour une résolution par
  PCG et BoomerAMG ; une résolution par cuDSS ne prend ensuite que 8 ms.
- Avec NewtonK0, K0 n'est factorisé qu'une fois (0,6 s), mais chaque
  résolution tangente par le gradient conjugué préconditionné par K0
  coûte 67 ms (descente et remontée triangulaires à chaque itération du
  gradient conjugué), contre 22 ms avec un cycle de BoomerAMG.
- La meilleure variante reste NewtonK0 avec un cycle de BoomerAMG et la
  tolérance d'Eisenstat et Walker (1,85 s ; 2,14 s le 10/10 avec mfem-mgis
  517c0c0b08, `results/L5-beam-2026-10-10-0429`).
- cuDSS est le bon choix pour les petits systèmes (ex7, 2 634 inconnues ;
  ex8, 9 516 et 28 548 inconnues : `results/cudss-2026-10-12`), où les
  résolutions itératives sur le GPU sont dominées par les lancements de
  noyaux, pas pour les grands systèmes tridimensionnels, où le remplissage
  de la factorisation coûte plus que l'AMG.
