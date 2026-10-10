# Newton inexact sur la poutre (10/10/2026)

Poutre r4 du prototype (111 843 inconnues, von Mises à écrouissage
linéaire, 10 pas), tous les noyaux, prédiction élastique. mfem-mgis
517c0c0b08. Script : `cas/mesures-L5-ew.sh`. Une RTX 5080, 16 cœurs.
Charge 3,7 au départ (fin de compilations), résultats de Newton et de
NewtonK0 égaux à ceux de `results/L5-poutre-2026-10-10-0250/` à 1 % près.

Tolérance du solveur linéaire fixe (1e-4, aussi pour les prédictions) ou
adaptée à la décroissance du résidu (option `-alst`, option
`AdaptiveLinearSolverTolerance` du Newton et de la méthode `NewtonK0` du
point fixe : choix 1 et 2 d'Eisenstat et Walker, réglages par défaut de
MFEM, tolérance initiale 0,5). Newton : PCG de hypre et BoomerAMG
reconstruit à chaque itération. NewtonK0 : gradient conjugué de MFEM
préconditionné par BoomerAMG construit une fois sur K0, prédiction par ce
gradient conjugué (`-fpm newton-k0 -els amg -k0 -pk0`).

| Cas | GPU | CPU, 8 processus |
|---|---|---|
| Newton, tolérance fixe | 3,12 s (32 it.) | 13,59 s (29 it.) |
| Newton, Eisenstat-Walker 1 | 3,01 s (53 it.) | 15,73 s (51 it.) |
| Newton, Eisenstat-Walker 2 | 3,32 s (63 it.) | 17,70 s (61 it.) |
| NewtonK0, tolérance fixe | 2,57 s (30 it.) | 11,56 s (29 it.) |
| NewtonK0, Eisenstat-Walker 1 | **2,14 s** (52 it.) | **10,93 s** (51 it.) |
| NewtonK0, Eisenstat-Walker 2 | 2,19 s (58 it.) | 12,30 s (61 it.) |

Temps de la boucle sur les pas, itérations de Newton au total. Pour
comparaison, le meilleur point fixe mesuré sur 8 cœurs (Anderson 8, un
cycle d'AMG) prend 9,5 s (`results/L5-poutre-2026-10-10-0250/`).

Résolution (temps et appels, `cas/resume-profilage.py`) :

| Cas | GPU | CPU, 8 processus |
|---|---|---|
| Newton, tolérance fixe | 2,01 s (32) | 8,73 s (29) |
| Newton, Eisenstat-Walker 1 | 1,78 s (53) | 9,37 s (51) |
| NewtonK0, tolérance fixe | 1,60 s (30) | 7,40 s (29) |
| NewtonK0, Eisenstat-Walker 1 | 1,05 s (52) | 5,21 s (51) |

## Constats

- **Sur le GPU, NewtonK0 avec Eisenstat-Walker 1 est le meilleur solveur
  mesuré** : 2,14 s, ×1,46 contre Newton et ×1,20 contre NewtonK0 à
  tolérance fixe. La résolution passe de 1,60 s à 1,05 s malgré 22
  itérations de Newton de plus : chaque gradient conjugué est beaucoup plus
  court, et une itération de plus ne coûte qu'une intégration et un
  assemblage par noyaux (BoomerAMG n'est pas reconstruit).
- **Avec Newton classique, Eisenstat-Walker ne paie pas** : chaque
  itération de plus reconstruit BoomerAMG et assemble la tangente ; sur
  8 cœurs, le calcul est 16 % (choix 1) à 30 % (choix 2) plus lent.
- **Choix 1 contre choix 2** : le choix 1 (résidu du système linéaire,
  un produit par la tangente de plus par itération) fait moins
  d'itérations de Newton (51 à 53 contre 58 à 63) et gagne partout.
- **Pas élastiques** : avec la tolérance fixe, la prédiction et une
  itération suffisent ; avec Eisenstat-Walker, la tolérance initiale (0,5)
  donne 4 ou 5 itérations aux pas élastiques. Une tolérance initiale plus
  basse (option `InitialTolerance`) est une piste pour gagner encore.
- **Sur CPU**, le point fixe accéléré reste le meilleur (9,5 s), NewtonK0
  avec Eisenstat-Walker 1 suit (10,9 s).

## Tolérance initiale (NewtonK0, Eisenstat-Walker 1)

Option `-alsti` de la poutre (`InitialTolerance`), mêmes réglages,
sorties dans `tolerance-initiale/` (exécutées à la main, charge 3,7 à
4,6 pendant les essais sur 8 cœurs).

| Tolérance initiale | GPU | CPU, 8 processus |
|---|---|---|
| 0,5 (défaut) | 2,19 s (52 it.) | 10,91 s (51 it.) |
| 0,1 | **1,96 à 2,02 s** (37 it.) | 9,66 s (38 it.) |
| 0,01 | 1,96 s (37 it.) | **9,46 s** (38 it.) |
| 0,001 | 2,07 s (36 it.) | |

- Une tolérance initiale de 0,1 ou 0,01 rend aux pas élastiques leurs
  deux itérations : sur le GPU, 1,96 s, soit ×1,59 contre Newton (3,12 s) ;
  sur 8 cœurs, 9,46 s, au niveau du meilleur point fixe (9,5 s). Sur le
  GPU, deux exécutions avec 0,1 donnent les mêmes itérations (1,96 s et
  2,02 s).
- La valeur par défaut (0,5, celle de MFEM) reste conservatrice ; le choix
  0,1 est un bon réglage pour ce cas, à confirmer sur le polycristal.
