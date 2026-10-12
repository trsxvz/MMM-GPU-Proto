# Réutilisation de l'analyse de MUMPS (12/10/2026)

Paramètre `ReuseAnalysis` du solveur `MUMPSSolver` de mfem-mgis
(branche mumps-reuse, c020804d68 et a114b6c870 après relecture) : l'analyse de MUMPS (ordonnancement et
factorisation symbolique) de la matrice précédente est réutilisée quand la
nouvelle matrice a le même motif sur tous les processus. Les matrices des
itérations successives de Newton ont le même motif : seule la première est
analysée. `mfem::MUMPSSolver::SetReorderingReuse` (MFEM 4.10) ne vérifie
pas le motif : mfem-mgis le compare (lignes, colonnes des blocs diagonal
et hors diagonal, indices globaux des colonnes hors diagonale), et refait
l'analyse sur tous les processus si le motif a changé sur l'un d'eux. La
fonction `getNumberOfMUMPSAnalyses` donne le nombre d'analyses faites. La
matrice étant distribuée, l'analyse de MUMPS ne dépend que du motif : la
réutiliser ne change pas la factorisation.
Option `-mra` d'ex8 (mfem-mgis-examples 2720017).

## Matrice d'élasticité d'ex8 seule

Programme `programme/mumps-reuse.cpp` (MFEM seul) : élasticité linéaire
sur le maillage d'ex8 (28 548 inconnues, éléments linéaires, rangement par
nœuds, tous les bords encastrés), MUMPS symétrique défini positif, six
factorisations de matrices de même motif et de valeurs différentes. Temps
moyen des factorisations après la première, en secondes
(`variantes-matrice-elasticite.txt`) :

| Variante | Sans réutilisation | Avec réutilisation |
|---|---|---|
| 1 processus, ordonnancement automatique | 0,258 | 0,207 |
| 1 processus, METIS, Scotch ou PORD | 0,255 à 0,256 | 0,205 à 0,206 |
| 1 processus, AMF | 0,380 | 0,359 |
| 1 processus, AMD | 0,496 | 0,476 |
| 1 processus, BLR à 1e-10, 1e-8 ou 1e-6 | 0,277 à 0,291 | 0,226 à 0,238 |
| 2 processus, ordonnancement automatique | 0,189 | 0,136 |

- La réutilisation fait gagner 20 % sur un processus et 28 % sur deux.
- L'ordonnancement automatique de MUMPS est le meilleur ; la factorisation
  approchée BLR est plus lente à cette taille.
- Les ordonnancements parallèles (ParMETIS, PT-Scotch) ne sont pas dans
  le MUMPS de l'environnement : MUMPS renvoie l'erreur -38 et MFEM
  s'arrête.

Seule la réutilisation de l'analyse est donc proposée dans mfem-mgis.

## ex8 complet

5 pas, maillage d'origine, tout par noyaux, MUMPS symétrique défini positif
pour la mécanique (`-msMc SymmetricPositiveDefinite`), GMRES et BoomerAMG
pour la thermique. Construction `build/examples-mumps`, liée à
`install/mfem-mgis-mumps`. Deux passages alternés par variante
(`resume-np2.txt`, `resume-gpu.txt`). Temps en secondes :

| Cas | Variante | Total | Couplage | Factorisations (29) |
|---|---|---|---|---|
| CPU, 2 processus (cœurs 0 à 3), lois generic | sans | 20,42 et 20,63 | 20,30 et 20,52 | 5,47 et 5,57 |
| | avec `-mra` | 19,15 et 19,21 | 19,03 et 19,08 | 4,08 et 4,07 |
| GPU, 1 processus (cœurs 4 à 7), lois cuda, hypre sur l'hôte | sans | 11,15 et 11,08 | 10,96 et 10,90 | 7,23 et 7,22 |
| | avec `-mra` | 9,73 et 9,73 | 9,55 et 9,55 | 5,87 et 5,86 |

- Les factorisations prennent 26 % de moins sur CPU (2 processus) et 19 %
  de moins sur le GPU (un processus, MUMPS sur l'hôte). ex8 sur le GPU
  passe de 11,1 à 9,7 s.
- Mêmes statistiques des champs (température, norme du déplacement) et
  mêmes nombres d'itérations de Newton (35 en thermique, 29 en mécanique)
  avec et sans réutilisation. Les résidus diffèrent aux arrondis, autant
  qu'entre deux passages identiques sans réutilisation.
- Sur le GPU, MUMPS reste le premier poste (60 % du temps).

## Tests

`MUMPSReuseAnalysisTest` (1 et 2 processus, et sur le GPU) : suite de
matrices dont le motif est gardé, change sur tous les processus, change
sur le dernier processus seulement, ne change que par les colonnes des
blocs hors diagonale, ou que par les indices globaux de ces colonnes,
avec et sans réutilisation, symétriques ou non. Après chaque matrice, le
test vérifie le résidu et le nombre d'analyses : une nouvelle analyse si
et seulement si le motif a changé. Sept mutations détectées (aucune
comparaison, pas de réduction sur les processus, indices globaux ou
colonnes hors diagonale non comparés, jamais de réutilisation, motif
oublié après une réutilisation, analyses non comptées). Relecture
indépendante Z8 traitée.

## Fichiers

- `ex8-np2-noyaux-{sans,avec}-{1,2}.txt`, `ex8-gpu-noyaux-{sans,avec}-{1,2}.txt` :
  sorties brutes avec les tableaux de profilage.
- `mesure-ex8-cpu.sh`, `mesure-ex8-gpu.sh` : commandes des mesures.
- `programme/` : programme de mesure sur la matrice d'élasticité
  (`make` avec `scripts/env.sh`, puis
  `./mumps-reuse assemblage_hexa.msh 6 [ordonnancement] [tolérance BLR]`).
