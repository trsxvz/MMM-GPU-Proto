# AmgX contre BoomerAMG pour NewtonK0 (10/10/2026)

Poutre r4 du prototype (111 843 inconnues, 10 pas), tous les noyaux,
mfem-mgis de2e76a913 (enveloppe d'AmgX propre à mfem-mgis, AmgX patché
pour CUDA 13). Script : `cas/mesures-L5-amgx.sh`. Une RTX 5080, 16 cœurs.
Le balayage et les répétitions du 11/10 (`balayage-k0`,
`repetitions-k0.txt`) ont été faits avec mfem-mgis 309849dff3, qui
téléverse la matrice d'un seul processus comme une matrice locale (et
non distribuée) : mêmes nombres d'itérations (110 et 111).

Paramètres d'AmgX passés par `-amgxo` avec leurs noms dans AmgX. « Réglé » :
`strength_threshold=0.8, interpolator=MULTIPASS` (meilleur réglage d'un
balayage sur K0, voir plus bas). « Simple » : calculs d'AmgX en simple
précision (`-amgxs`, mode `dFFI`).

## Opérateur K0 seul (`poutre-k0`, GPU)

Gradient conjugué préconditionné, tolérance relative 1e-8. Même
discrétisation et même loi élastique que la poutre. Temps d'une
résolution, mise en place à part.

| Solveur | Mise en place | Résolution | Itérations |
|---|---|---|---|
| CG de MFEM + BoomerAMG | dans la 1re résolution (+25 ms) | 70 ms (69 à 89 ms) | 70 |
| CG de MFEM + BoomerAMG agressif | dans la 1re résolution (+30 ms) | 70 ms | 80 |
| PCG de hypre + BoomerAMG | dans la 1re résolution (+22 ms) | 65 ms | 70 |
| CG de MFEM + AmgX par défaut | 190 ms | 1,95 s | 1078 |
| CG de MFEM + AmgX par défaut, simple | 176 ms | 1,10 s | 1110 |
| CG de MFEM + AmgX réglé | 71 ms | 92 ms (75 à 92 ms) | 111 |
| CG de MFEM + AmgX réglé, 2 balayages | 66 ms | 108 à 123 ms | 96 |
| **CG de MFEM + AmgX réglé, simple** | 67 ms | **53 à 61 ms (54 à 72 ms)** | 110 |
| PCG d'AmgX réglé | 63 ms | 78 ms | 116 |
| CG de MFEM + Jacobi | 1 ms | 227 ms | 1862 |

Entre parenthèses : trois nouvelles exécutions du 11/10 (`repetitions-k0.txt`).
D'une résolution à l'autre, le temps varie de 15 % environ sur ce GPU,
pour un même nombre d'itérations : seuls des écarts plus grands sont
significatifs. L'ordre ne change pas : AmgX réglé en simple précision
devant BoomerAMG et AmgX réglé en double.

Le réglage par défaut (celui de MFEM : seuil 0,25, interpolation D2) est
un mauvais préconditionneur pour cette élasticité vectorielle, traitée
comme une matrice scalaire. Le seuil 0,8 et l'interpolation MULTIPASS le
rendent comparable à BoomerAMG, dont la stratégie « System » tient compte
des trois composantes. Une itération d'AmgX réglé coûte 0,83 ms en double
et 0,48 à 0,55 ms en simple, contre 1,0 ms pour BoomerAMG : le FP64 de la
RTX 5080 est lent, la simple précision accélère le multigrille de 1,5 à
1,7 fois sans changer le nombre d'itérations.

Balayage sur K0 (`cas/balayage-amgx-k0.sh`, dossier `balayage-k0` :
configurations JSON, sorties brutes et `resume.txt`, refait le 11/10 avec
la version actuelle de `poutre-k0`). Itérations et temps d'une
résolution :

| Configuration | Itérations | Résolution |
|---|---|---|
| classique par défaut (seuil 0,25, D2) | 1078 | 1,8 à 1,9 s |
| classique par défaut, 2 balayages | 779 | 1,9 s |
| classique par défaut, agressif (1 niveau) | 1283 | 0,75 s |
| seuil 0,5 | 1071 | 1,0 à 1,2 s |
| seuil 0,5, 2 balayages | 769 | 1,15 à 1,18 s |
| seuil 0,6 | 128 | 165 ms |
| seuil 0,6, 2 balayages | 102 | 197 à 203 ms |
| seuil 0,7 | 116 | 130 ms |
| seuil 0,8 | 118 | 107 ms |
| seuil 0,85 | 124 | 99 ms |
| seuil 0,9 | 133 | 102 ms |
| **seuil 0,8, MULTIPASS** | **111** | **75 à 77 ms** |
| MULTIPASS, seuil 0,6 ou 0,7 | 129 à 130 | 90 à 93 ms |
| MULTIPASS, seuil 0,7, 2 balayages | 107 | 118 à 119 ms |
| MULTIPASS, seuil 0,8, 2 balayages | 96 | 101 à 104 ms |
| MULTIPASS, seuil 0,9 | 125 | 82 à 85 ms |
| seuil 0,8, 2 balayages | 98 | 140 ms |
| seuil 0,8, HMIS | 117 | 178 ms |
| seuil 0,8, agressif (1 niveau) | 252 | 170 ms |
| seuil 0,5, cycle F, W, LU dense ou max_row_sum 1 | 1055 à 1588 | 1,0 à 4,7 s |
| classique, interpolation D1 | 1602 | 4,5 s |
| agrégation, paires multiples | 925 | 1,15 s |
| agrégation SIZE_2, Gauss-Seidel multicolore | 458 à 483 | 5,7 à 6,0 s |
| agrégation, Jacobi par blocs (SIZE_2, SIZE_8) | non convergé | |
| classique, Gauss-Seidel multicolore | 806 à 6920 | 20 à 166 s |

Le seuil de connexion décide de tout : jusqu'à 0,5 le multigrille
classique est presque inutile, à partir de 0,6 le nombre d'itérations
est divisé par huit. Deux balayages du lisseur ôtent 13 à 28 %
d'itérations mais coûtent plus qu'ils ne rapportent. MULTIPASS, seuil
0,8, 2 balayages donne ici 101 à 104 ms, contre 108 à 123 ms dans le
tableau précédent (même configuration passée par `-amgxo`) : écart dans
la variabilité des temps. Avec le Jacobi par
blocs, l'agrégation rend un préconditionneur que le gradient conjugué ne
supporte pas (arrêt en 3 ou 4 itérations, résidu relatif 0,58 à 0,69).

## Calcul complet, NewtonK0 + Eisenstat-Walker 1

| Multigrille de K0 | GPU | Itérations | CPU, 8 processus | Itérations |
|---|---|---|---|---|
| BoomerAMG par défaut | 1,99 s | 33 | 9,23 s | 38 |
| BoomerAMG agressif | 1,77 à 1,87 s | 34 à 36 | **7,40 s** | 34 |
| AmgX par défaut | 51,5 s | 90 | 849 s | 92 |
| AmgX réglé | 1,95 à 2,04 s | 36 à 38 | 133 s | 40 |
| AmgX réglé, 2 balayages | 2,39 s | 38 | 167 s | 44 |
| **AmgX réglé, simple** | **1,58 à 1,60 s** | 37 | 128 s | 39 |

Temps de la boucle sur les pas ; sur le GPU, quatre exécutions (celle du
script et trois répétitions, `repetitions-gpu.txt`) pour les lignes en
intervalle et pour AmgX réglé en simple précision. Même solution dans tous les
cas (normes finales à moins de 1e-7 en relatif).

- Sur le GPU, AmgX réglé en simple précision est le meilleur résultat de
  la poutre : 1,58 à 1,60 s contre 1,77 à 1,87 s pour le meilleur
  BoomerAMG (-10 à -15 %).
- Sur 8 processus CPU, AmgX est inutilisable : chaque processus a ses
  ressources AmgX sur le même GPU. Cause probable, non mesurée : huit
  contextes CUDA se partagent le GPU par tranches de temps (pas de MPS),
  et les échanges de halo passent par l'hôte. AmgX n'a d'intérêt ici
  qu'avec un GPU par processus.

## Suites possibles

- Matrices par blocs 3×3 (ordre par composantes) : multigrille d'AmgX
  adapté aux systèmes, à comparer à la stratégie « System » de BoomerAMG.
- Simple précision pour BoomerAMG (hypre en précision mixte) et pour le
  multigrille de Ginkgo.
- Mise en place d'AmgX sans passer par l'hôte, réutilisation
  (`AMGX_solver_resetup`) pour un K0 mis à jour.
