# Tangente sans assemblage sur le GPU (10/10/2026)

Poutre r4 du prototype (111 843 inconnues, 10 pas), loi `cuda`, une RTX
5080, noyaux pour les gradients et les forces intérieures. mfem-mgis
4b72095ee0. Script : `cas/mesures-L5-sans-matrice.sh`. Machine calme
(charge 1,7 au départ, en baisse).

Tangente appliquée sans assemblage (`TangentOperatorComputation` à
`MatrixFree`, action par le noyau fusionné) avec trois préconditionneurs
possibles :

- le lisseur de Jacobi ;
- le lisseur de Chebyshev de MFEM (`OperatorChebyshevSmoother`, ordres 2 à 5) ;
- NewtonK0 : gradient conjugué préconditionné par BoomerAMG construit une
  fois sur K0, avec ou sans Eisenstat-Walker 1 et une tolérance initiale
  de 0,1.

Ces cas sont comparés à NewtonK0 avec la tangente assemblée par noyaux.
Tolérance du solveur linéaire : 1e-4.

| Cas | Temps | Itérations de Newton | Résolution |
|---|---|---|---|
| sans matrice, CG + Jacobi | 21,04 s | 30 | 19,50 s |
| sans matrice, CG + Chebyshev 2 | 23,13 s | 32 | 21,41 s |
| sans matrice, CG + Chebyshev 3 | 25,50 s | 30 | 23,72 s |
| sans matrice, CG + Chebyshev 4 | 28,22 s | 30 | 26,33 s |
| sans matrice, CG + Chebyshev 5 | 31,43 s | 32 | 29,45 s |
| sans matrice, NewtonK0 | 3,05 s | 33 | 2,24 s |
| sans matrice, NewtonK0 + Eisenstat-Walker 1 | 2,05 s | 37 | 1,24 s |
| assemblée, NewtonK0 + Eisenstat-Walker 1 | **1,99 s** | 37 | 0,97 s |

Temps de la boucle sur les pas ; résolution d'après
`cas/resume-profilage.py`.

## Constats

- **Chebyshev ne paie pas ici** : chaque ordre de plus ajoute une action de
  la tangente sans matrice par application du préconditionneur (tétraèdres
  d'ordre 1, action à 3 x 4 = 12 inconnues par élément), et le temps croît
  avec l'ordre. Les itérations du gradient conjugué ne sont pas
  enregistrées par ces mesures : la réduction qu'apporte chaque ordre n'est
  pas chiffrée ici (sur un Laplacien, 46 itérations avec Jacobi, 27 et 17
  avec Chebyshev d'ordre 2 et 4, test `OperatorSmoothersTest`). Le lisseur
  reste utile comme brique (lissage d'une multigrille géométrique ou p,
  ordres élevés) : il est dans la fabrique et testé.
- **Sans matrice, NewtonK0 avec Eisenstat-Walker 1 est à 3 % de la tangente
  assemblée** (2,05 s contre 1,99 s, une exécution par cas) : la tangente
  n'est plus assemblée pendant les itérations. La mémoire n'est pas réduite
  pour autant : K0 (copie d'une tangente assemblée), la matrice du problème
  et la matrice de la prédiction restent assemblées.
- Eisenstat-Walker apporte un peu plus sans matrice (3,05 s vers 2,05 s,
  -33 %) qu'avec la tangente assemblée (2,57 s vers 1,96 s, -24 %,
  `results/L5-poutre-2026-10-10-0429/`).
