# Assemblage direct de la tangente par noyaux, poutre (09/10/2026)

Mesures intermédiaires du lot L3 : la tangente n'est plus formée de
matrices élémentaires stockées puis recopiées par la boucle de MFEM. Un
noyau calcule chaque bloc `N x N` d'une paire de nœuds d'un élément et
l'ajoute directement dans la matrice creuse locale (position cherchée dans
la ligne, addition atomique). Le motif de la matrice est celui du premier
assemblage de MFEM. La fin de `GetGradient` est celle de MFEM (bloc
diagonal et `PtAP` en parallèle, `RAP` en séquentiel).

Poutre r4 (111 843 inconnues, tétraèdres d'ordre 1), PCG + BoomerAMG
(tolérance 1e-4). Script : `cas/mesures-L3.sh poutre` (10 pas, sur CPU) ;
le dossier `essais-3-pas` contient des calculs de 3 pas à 1 processus sur
CPU et sur GPU (`--device cuda`), avec des chronomètres provisoires (lignes
`TMP`).

Réserve : les mesures sur 10 pas ont été faites avec un navigateur ouvert
(charge de 1 à 2) ; les temps totaux varient de 10 à 15 % entre calculs
identiques. Le résumé indique le commit de base (cf7f764628), les
modifications n'étant pas encore validées.

## Sur CPU (10 pas)

| Cas | Total | Intégration | Résidu | Tangente | Résolution |
|---|---|---|---|---|---|
| 8 processus, lots | 17,65 | 0,55 (41) | 0,48 (41) | 2,36 (31) | 10,99 |
| 8 processus, tangente par noyaux | 20,13 | 0,71 (39) | 0,64 (39) | 3,25 (29) | 11,70 |
| 8 processus, tout par noyaux | 18,46 | 0,22 (39) | 0,13 (39) | 3,04 (29) | 11,28 |
| 1 processus, lots | 66,61 | 3,22 (38) | 3,00 (38) | 11,72 (28) | 34,50 |
| 1 processus, tangente par noyaux | 70,29 | 3,26 (38) | 3,04 (38) | 14,27 (28) | 35,35 |
| 1 processus, tout par noyaux | 67,11 | 0,96 (38) | 0,65 (38) | 14,31 (28) | 36,23 |

L'assemblage direct ne gagne rien sur CPU : la tangente par noyaux reste
22 % plus lente qu'élément par élément, comme avec les matrices
élémentaires stockées. Les chronomètres montrent pourquoi : le noyau seul
prend 0,48 s par appel, contre 0,34 s pour toute la boucle de MFEM
(matrices élémentaires et `AddSubMatrix`). Un fil par paire de nœuds
recalcule l'inverse du jacobien et relit les blocs de la tangente pour
chaque paire : bonne granularité pour un GPU, mauvaise pour un cœur de
CPU.

Le résidu par noyaux est divisé par 4,6 à 1 processus (boucle de MFEM
sautée), les gradients par 3,4.

## Sur GPU (3 pas, 1 processus)

| Cas | Total | Intégration | Résidu | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|---|
| CPU, lots | 24,53 | 1,14 (13) | 1,05 (13) | 4,47 (10) | 13,36 | 4,32 |
| GPU, lots | 29,05 | 6,70 (13) | 6,91 (13) | 9,32 (10) | 0,89 | 4,96 |
| GPU, tout par noyaux, `EliminateRowsCols` | 9,05 | 0,47 (13) | 0,02 (13) | 2,38 (10) | 0,83 | 5,08 |
| GPU, tout par noyaux, `EliminateBC` | 7,84 | 0,49 (14) | 0,02 (14) | 1,26 (11) | 0,90 | 4,92 |

- hypre est compilé pour le GPU : la résolution (PCG + BoomerAMG) passe de
  13,4 s à 0,9 s.
- Sans noyaux, les boucles d'éléments de MFEM avec des vecteurs sur le GPU
  sont plus lentes que sur CPU (copies).
- Le noyau d'assemblage prend 7 ms par appel (0,48 s sur un cœur de CPU).
- Chronomètres de la fin de `GetGradient` (`cuda-np1-lots-noyaux-chronometres.txt`) :
  bloc diagonal 1 ms, `PtAP` 1,3 ms, élimination des ddl imposés 122 ms.
  `HypreParMatrix::EliminateRowsCols` recopie la matrice sur l'hôte,
  élimine en calculant la partie éliminée (inutilisée), puis la renvoie
  sur le GPU. `HypreParMatrix::EliminateBC` (`DIAG_ONE`) donne la même
  matrice, le motif étant symétrique, en place et sur le GPU : la tangente
  passe de 2,38 s à 1,26 s, dont environ 0,9 s pour le premier appel, qui
  passe encore par la boucle de MFEM pour construire le motif.
- Le nombre d'itérations varie d'un calcul à l'autre sur GPU (10 ou 11) :
  additions atomiques et réductions du GPU, ordre des sommes non
  déterministe. Même norme de la solution à 1e-10 près.

## Motif construit à partir des éléments (mfem-mgis 799941c8fd)

Le motif de la matrice locale est construit à partir des ddl des éléments
au premier appel : plus de boucle de MFEM. Dans chaque ligne, le terme
diagonal vient d'abord, puis les autres colonnes triées : c'est l'ordre
que `HypreParMatrix` impose en place aux lignes de la matrice locale qu'il
partage (`hypre_CSRMatrixReorder`). Le noyau cherche les colonnes par
dichotomie au lieu d'un parcours de la ligne.

| Cas (3 pas, 1 processus) | Total | Intégration | Résidu | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|---|
| CPU, lots | 24,53 | 1,14 (13) | 1,05 (13) | 4,47 (10) | 13,36 | 4,32 |
| CPU, tout par noyaux | 22,62 | 0,37 (13) | 0,22 (13) | 4,00 (10) | 13,32 | 4,53 |
| GPU, tout par noyaux | 6,84 | 0,50 (13) | 0,02 (13) | 0,16 (10) | 0,89 | 4,99 |

- Sur CPU, la tangente par noyaux devient plus rapide que la boucle de
  MFEM (0,40 s par appel contre 0,45 s) : le coût venait surtout du
  parcours linéaire des lignes, pas du recalcul par paire de nœuds.
- Sur GPU, 16 ms par tangente (contre 126 ms avec le premier motif de
  MFEM). Le calcul est dominé par la prédiction (73 %), qui assemble des
  formes bilinéaires élément par élément sur l'hôte.

## Suite

- prédiction par noyaux (lot L4) ;
- loi sur le GPU (bibliothèque `generic-parallel` CUDA) ;
- test MPI des noyaux dans mfem-mgis ;
- variante du noyau à un fil par élément, pour comparer sur CPU.

## Loi intégrée sur le GPU (mfem-mgis ee1f31fa1c, MGIS a0b5140)

Loi `generic-parallel` compilée pour `cuda`, intégrée par MGIS sur les
copies des tableaux du matériau sur le GPU (vue de MGIS sur les copies
de MFEM), propriétés et variables externes copiées sur le GPU à chaque
appel ; tous les noyaux. Poutre r4, 3 pas, 1 processus :

| Cas | Total | Intégration | Résidu | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|---|
| GPU, loi `generic` (hôte) | 6,51 | 0,47 (13) | 0,02 (13) | 0,15 (10) | 0,76 | 4,85 |
| GPU, loi `cuda` | 6,25 | 0,24 (13) | 0,01 (13) | 0,12 (10) | 0,77 | 4,80 |

Même solution. Il reste la prédiction (77 % du temps), qui assemble ses
opérateurs élément par élément sur l'hôte.
