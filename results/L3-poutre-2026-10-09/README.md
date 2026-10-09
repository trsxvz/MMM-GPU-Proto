# Mesures du lot L3 sur la poutre (09/10/2026)

Poutre r4 (111 843 inconnues), PCG + BoomerAMG (tolérance 1e-4), 10 pas,
sur CPU. Intégration par lots avec les gradients, les forces intérieures
(résidu) et les matrices élémentaires de la tangente calculés élément par
élément ou par noyaux (mfem-mgis ed91564643). Script :
`cas/mesures-L3.sh poutre`.

Réserve : une relecture indépendante exécutait des tests en même temps ;
le calcul point par point à 1 processus, lancé le premier, est plus lent
partout (résolution comprise) et n'est pas comparable. Les comparaisons se
font avec « lots ».

| Cas | Total | Intégration | Résidu | Tangente | Résolution |
|---|---|---|---|---|---|
| 8 processus, lots | 15,69 | 0,51 (41) | 0,45 (41) | 2,24 (31) | 9,56 |
| 8 processus, gradients par noyaux | 15,17 | 0,20 (41) | 0,45 (41) | 2,22 (31) | 9,37 |
| 8 processus, résidu par noyaux | 15,34 | 0,50 (41) | 0,21 (41) | 2,23 (31) | 9,48 |
| 8 processus, tangente par noyaux | 14,93 | 0,47 (38) | 0,40 (38) | 2,51 (28) | 8,61 |
| 8 processus, tout par noyaux | 14,32 | 0,19 (38) | 0,20 (38) | 2,51 (28) | 8,47 |
| 1 processus, lots | 65,71 | 3,18 (38) | 2,96 (38) | 11,79 (28) | 33,85 |
| 1 processus, gradients par noyaux | 62,48 | 0,94 (38) | 2,99 (38) | 11,64 (28) | 33,25 |
| 1 processus, résidu par noyaux | 64,04 | 3,21 (38) | 1,40 (38) | 11,71 (28) | 33,78 |
| 1 processus, tangente par noyaux | 68,04 | 3,22 (38) | 2,98 (38) | 14,62 (28) | 33,38 |
| 1 processus, tout par noyaux | 64,12 | 1,01 (38) | 1,35 (38) | 14,55 (28) | 33,39 |

Entre parenthèses : nombre d'appels. Même solution partout (norme
48,83704586). À 8 processus, la tangente par noyaux change l'arrondi et le
nombre d'itérations de Newton (28 au lieu de 31), d'où un temps total plus
bas qui ne vient pas des noyaux.

## Constats

- **Gradients par noyaux** : phase d'intégration divisée par 2,5 (8
  processus) à 3,4 (1 processus) ; il reste l'appel de la loi.
- **Résidu par noyaux** : divisé par 2,1, alors que la boucle de MFEM sur
  les éléments tourne encore (elle ne rend que des zéros). Elle est
  supprimée quand le formulaire n'a pas d'autre intégrateur (commit
  suivant).
- **Tangente par noyaux** : plus lente sur CPU (+24 % par appel) : calculer
  toutes les matrices élémentaires, les stocker puis les recopier dans
  l'assemblage de MFEM coûte plus que de les calculer dans sa boucle. Elle
  ne vaudra que sur GPU avec un assemblage sur le périphérique, à écrire
  (`ElementRestriction::FillJAndData` de MFEM ne traite que les espaces
  scalaires), ou sans matrice (option a du plan).
- **Premier essai, avant correction** (dossier `scratch/`) : les noyaux
  utilisaient `Mesh::GetGeometricFactors`, qui donne des nœuds au maillage
  (`EnsureNodes`) ; toutes les transformations d'éléments passaient
  ensuite par ces nœuds, et les phases restées élément par élément
  étaient 50 % plus lentes. Les facteurs géométriques sont maintenant
  calculés sur une copie des sommets.
