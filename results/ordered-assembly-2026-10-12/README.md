# Assemblage ordonné de la tangente (12 et 13/10/2026)

Paramètre `TangentOperatorAssembly` de mfem-mgis (`Atomic` par défaut,
`Ordered`), variable `MFEM_MGIS_TANGENT_OPERATOR_ASSEMBLY`. Avec un
périphérique, les matrices de rigidité des éléments sont ajoutées à la
matrice creuse par additions atomiques, dans un ordre qui change d'un
passage à l'autre : les derniers chiffres des résultats changent aussi, et
le nombre d'itérations de Newton peut varier (32 ou 31 tangentes ci-dessous
avec l'assemblage atomique). L'assemblage ordonné stocke les matrices des
éléments, puis les ajoute dans un ordre fixe. Sans périphérique, le chemin
de l'hôte ajoute déjà les éléments l'un après l'autre : les deux valeurs du
paramètre y donnent le même assemblage.

Poutre r4 du prototype (111 843 inconnues, von Mises à écrouissage
linéaire, 10 pas), tout par noyaux, loi cuda, un processus sur le GPU
(cœurs 0 à 3), hypre sur le GPU. Installation séparée
`install/mfem-mgis-ordered`, construction `build/cas-ordered` (seul
l'exécutable `poutre`). Commandes dans `measure.sh`.

## Deux versions de l'assemblage ordonné sur le GPU

- Première version (mfem-mgis 456c298dd6) : un fil par ligne de la
  matrice, qui cherche la colonne de chacune de ses entrées dans la ligne.
  Profil Nsight Systems : 7,1 ms par assemblage pour ce noyau, contre
  0,48 ms pour le noyau des additions atomiques. Avec un seul fil par
  ligne (111 843 fils), chaque fil fait en série près de 200 recherches
  et lectures dispersées.
- Version retenue (mfem-mgis f57072220f) : un fil par valeur de la matrice
  (environ 9 millions de fils), qui ajoute dans un ordre fixe les deux à
  trois contributions de cette valeur. La liste des contributions est
  calculée une fois par agencement des lignes, à partir des positions des
  entrées (noyau existant), puis vérifiée colonne par colonne à chaque
  assemblage. Mêmes valeurs que la première version.

## Temps (secondes)

Temps total de la boucle sur les pas, deux passages par variante ; temps
cumulés de l'assemblage de la tangente (`NS::getJacobian`) et de
l'assemblage de la prédiction (`PRED::assemble`, qui contient le premier
assemblage, donc le calcul des contributions).

| Méthode | Assemblage | Total | Tangente | Prédiction |
|---|---|---|---|---|
| Newton, PCG et BoomerAMG | atomique | 3,09 et 3,09 | 0,137 (32) et 0,146 (33) | 0,104 et 0,106 |
| Newton, PCG et BoomerAMG | ordonné, par ligne | 3,20 et 3,20 | 0,354 (32) et 0,353 (32) | 0,178 et 0,175 |
| Newton, PCG et BoomerAMG | ordonné, par valeur | 3,09 et 3,67 | 0,189 (32) et 0,212 (32) | 0,234 et 0,247 |
| NewtonK0, Eisenstat et Walker | atomique | 1,89 et 1,95 | 0,123 (51) et 0,140 (53) | 0,106 et 0,109 |
| NewtonK0, Eisenstat et Walker | ordonné, par ligne | 2,40 et 2,61 | 0,493 (53) et 0,542 (53) | 0,173 et 0,183 |
| NewtonK0, Eisenstat et Walker | ordonné, par valeur | 2,22 et 2,35 | 0,218 (53) et 0,234 (53) | 0,235 et 0,244 |

Entre parenthèses, le nombre d'assemblages de la tangente. Les temps
totaux varient d'environ 10 % d'un passage à l'autre (autres sessions sur
la machine) ; les temps par assemblage sont plus stables.

- Par assemblage de la tangente : 4,3 ms en atomique, 11 ms par ligne,
  6 ms par valeur. Le premier assemblage coûte en plus 0,13 s par valeur
  (calcul des contributions, une fois par agencement des lignes ; 0,45 s
  quand les entrées étaient cherchées sur l'hôte).
- Mémoire en plus avec un périphérique : les matrices des éléments
  (8 n² octets par élément), les contributions (4 n² octets par élément,
  plus 4 octets par valeur de la matrice) et leurs positions, le temps de
  leur calcul (2 n² octets par élément).

## CPU

Même poutre sur 4 processus (cœurs 0 à 3), sans périphérique, loi openmp,
Newton, PCG et BoomerAMG. Passages alternés : atomique 17,96 et 16,36 s,
ordonné 17,84 et 16,71 s (tangente 0,76 et 0,70 s pour les deux). La
première version de l'assemblage ordonné, qui stockait les matrices aussi
sur l'hôte, prenait 22,2 s (tangente 4,1 à 4,3 s, sans les positions des
entrées) : fichiers `poutre-cpu4-newton-Ordered-lignes-*`.

## Reproductibilité

`FixedPointSolverTest` sur le GPU, noyaux forcés, assemblage ordonné :
quatre passages donnent des historiques d'itérations identiques bit à
bit, avec les deux versions ; avec l'assemblage atomique, trois passages
donnent trois historiques différents.
