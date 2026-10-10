# Point fixe sur la poutre, CPU et GPU (10/10/2026)

Poutre r4 (111 843 inconnues, tétraèdres d'ordre 1), plasticité de von
Mises à écrouissage linéaire, 10 pas, prédiction élastique, tolérance
relative du solveur non linéaire 1e-8. mfem-mgis d05475d265. Script :
`cas/mesures-L5.sh`. Une RTX 5080, 16 cœurs ; une relecture tournait en
même temps (charge moyenne de 4,5 au départ) : les temps sur CPU sont
majorés d'environ 20 % (Newton : 16,6 s sur machine calme, 20,2 s ici).

Tous les calculs utilisent les noyaux (gradients, résidu, tangente
assemblée) ; loi `generic` sur CPU (8 processus), loi `cuda` sur le GPU.

- **Newton** : meilleur Newton mesuré, PCG + BoomerAMG de hypre, tolérance
  linéaire 1e-4 ;
- **point fixe** sur l'opérateur élastique K0 (`FixedPoint`, méthode
  `ElasticOperator`, K0 assemblé une fois), accéléré par Anderson (profondeur
  8) ou par l'algorithme de Cast3M, K0 résolu par un cycle de BoomerAMG
  par itération (réglages par défaut de hypre, ou lisseur de Chebyshev,
  `RelaxationType=16`), ou par PCG + BoomerAMG ; prédiction résolue à 1e-8
  près ;
- **Newton modifié** accéléré par Anderson (tangente de la première
  itération du pas).

| Cas | CPU, 8 processus | GPU |
|---|---|---|
| Newton | 20,27 s (29 it.) | **3,43 s** (32 it.) |
| point fixe, Anderson 8, un cycle d'AMG | **13,33 s** (341 it.) | 4,80 s (1 231 it.) |
| idem, lisseur de Chebyshev | 13,25 s (330 it.) | **3,44 s** (438 it.) |
| point fixe, Cast3M, un cycle d'AMG | 14,36 s (477 it.) | 4,32 s (1 162 it.) |
| point fixe, Anderson 8, PCG + AMG | 31,52 s (94 it.) | 8,49 s (96 it.) |
| Newton modifié, Anderson 8 | 29,67 s (74 it.) | 7,82 s (76 it.) |

Temps total profilé (maximum sur les processus).

Répartition sur le GPU (secondes, nombre d'appels) :

| Cas | Intégration | Résidu | Résolution | Accélération | Prédiction |
|---|---|---|---|---|---|
| Newton | 0,05 (42) | 0,02 (42) | 2,05 (32) | | 0,78 (10) |
| point fixe, Chebyshev | 0,39 (449) | 0,23 (448) | 0,83 (438) | 0,12 (438) | 1,49 (10) |

## Constats

- **Sur CPU, le point fixe gagne** : 13,3 s contre 20,3 s pour Newton
  (×1,5, à comparer au ×1,8 du prototype du 04/10). Le point fixe évite
  l'assemblage de la tangente et la construction de BoomerAMG à chaque
  itération ; une itération (loi, résidu, un cycle d'AMG, Anderson) coûte
  environ 21 ms contre 560 ms pour une itération de Newton (temps hors
  prédiction divisés par le nombre d'itérations).
- **Sur le GPU, égalité** (3,44 s contre 3,43 s), et seulement avec le
  lisseur de Chebyshev : avec les réglages par défaut de hypre sur GPU
  (l1-Jacobi, grossissement PMIS), un cycle de BoomerAMG est bien plus
  faible (1 231 itérations au lieu de 438). Le prototype (×4 sur GPU)
  utilisait sa propre AMG par agrégation lissée avec les modes rigides.
  Chebyshev ne profite pas à Newton (3,95 s).
- **La prédiction pèse** : 1,49 s sur 3,44 s sur le GPU (6,2 s sur 13,3 s
  sur CPU), parce qu'elle est résolue à 1e-8 près par PCG + BoomerAMG, avec
  une nouvelle construction de BoomerAMG à chaque pas. Sans prédiction, le
  point fixe est plus lent (1 394 itérations sur le GPU, 23 s sur CPU) :
  Anderson et un cycle d'AMG résolvent mal le grand incrément élastique.
  Piste : faire la prédiction avec K0 et son AMG déjà construits (prévu
  au plan).
- Une itération du point fixe sur le GPU coûte environ 4 ms (loi 0,9 ms,
  résidu 0,5 ms, cycle d'AMG 1,9 ms, Anderson 0,3 ms) : la latence des
  lancements domine encore (voir `results/L6-poutre-2026-10-09/`).
- Les nombres d'itérations du point fixe sur le GPU varient d'une
  exécution à l'autre (de 1 128 à 2 133 pour une même configuration avec
  le lisseur par défaut) : hypre sur GPU et les additions atomiques ne
  sont pas reproductibles bit à bit, et Anderson amplifie ces écarts.

## Suite : prédiction par l'AMG de K0 (mfem-mgis 846551e0cb, 1ef05104d1)

Option `PredictionLinearSolver=ElasticOperator` du point fixe (`-pk0`) : la
prédiction élastique est résolue par le gradient conjugué préconditionné
par le cycle de BoomerAMG de K0, déjà construit, au lieu d'une nouvelle
construction de BoomerAMG à chaque pas. Point fixe, Anderson 8, un cycle
d'AMG de Chebyshev, K0 gardé, sur le GPU :

| Tolérance de la prédiction | Itérations | Total | Prédiction | Résolutions |
|---|---|---|---|---|
| 1e-4 | 534 | **2,80 s** | 0,77 s | 1,01 s |
| 1e-5 | 535 | 3,13 s | 1,10 s | 1,02 s |
| 1e-6 | 458 | 2,94 s | 1,18 s | 0,88 s |
| 1e-7 | 605 | 3,60 s | 1,31 s | 1,15 s |
| 1e-8 | 490 | 3,31 s | 1,41 s | 0,96 s |

Avec une prédiction à 1e-4, trois exécutions de chaque : Newton 3,15 s,
3,16 s, 3,18 s ; point fixe 3,32 s, 2,88 s, 2,97 s (566 à 674
itérations). **Sur le GPU, le point fixe est à égalité ou un peu plus
rapide que Newton** (médiane 2,97 s contre 3,16 s) ; sur 8 cœurs il
reste nettement plus rapide.

- La prédiction coûte environ 40 ms par pas en régime (assemblage par
  noyaux 3,3 ms, assemblage global 2,6 ms, gradient conjugué environ
  35 ms) ; le premier pas porte en plus les mises en place ponctuelles
  (géométrie des éléments 214 ms, motif de la matrice), que les moyennes
  par appel mélangent.
- Une itération du point fixe coûte environ 3,5 ms sur le GPU : intégration
  de la loi 0,96 ms (dont 0,65 ms de noyau de la loi), résidu 0,55 ms (le
  noyau des forces intérieures, un fil par nœud, inverse 4 fois chaque
  jacobienne : à écrire comme le noyau fusionné de l'action, un fil par
  élément), cycle d'AMG 1,9 ms, Anderson 0,3 ms.
- Réglages de BoomerAMG essayés en plus de Chebyshev (grossissement
  agressif, seuil 0,5, interpolation 18, niveaux limités) : aucun gain
  au-delà de la dispersion des nombres d'itérations.
