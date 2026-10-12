# cuDSS, solveur direct sur le GPU (12/10/2026)

Solveur `CuDSSSolver` de la fabrique de solveurs linéaires de mfem-mgis
(branche cudss, 840683d96f et 45afaa0653 après relecture), qui enveloppe
`mfem::CuDSSSolver` (MFEM 4.10, cuDSS 0.7.1). Paramètres `Symmetric`, `PositiveDefinite` (les modes
symétriques lisent le triangle inférieur, comme MUMPS) et `ReuseAnalysis`
(analyse réutilisée tant que le motif ne change sur aucun processus, comme
pour MUMPS, `results/mumps-analyse-2026-10-12`). Le solveur est refusé sans
périphérique sur un GPU CUDA. ex8 le choisit par `-lsMc CuDSSSolver`
(mfem-mgis-examples efeb8b8).

## Matrice d'élasticité d'ex8 seule

Programme `programme/cudss-bench.cpp` (MFEM seul) : élasticité linéaire
sur le maillage d'ex8 (28 548 inconnues, éléments linéaires, rangement par
nœuds, tous les bords encastrés), un processus, six matrices de même motif
et de valeurs différentes. Temps moyen des factorisations après la
première, en secondes (`matrice-elasticite-np1-hypre-gpu.txt`) :

| Solveur | Sans réutilisation | Avec réutilisation |
|---|---|---|
| MUMPS sur l'hôte, symétrique défini positif | 0,253 | 0,203 |
| cuDSS, matrice générale | 0,126 | 0,031 |
| cuDSS, symétrique défini positif (triangle inférieur) | 0,111 | 0,018 |

- Avec la réutilisation de l'analyse, cuDSS factorise 11 fois plus vite
  que MUMPS sur l'hôte.
- Résidus relatifs de l'ordre de 6e-15 dans tous les cas.
- Avec hypre sur l'hôte (`matrice-elasticite-np1-hypre-hote.txt`), mêmes
  ordres de grandeur (0,022 s pour cuDSS symétrique avec réutilisation) :
  MFEM copie la matrice vers cuDSS depuis la mémoire de l'hôte.
- Ce programme se bloque sur deux processus avant même la construction de
  l'espace d'éléments finis, avec la configuration par défaut de MFEM
  (mémoire de l'hôte par Umpire) : non analysé, sans lien avec cuDSS ; les
  tests de mfem-mgis sur deux processus et le GPU passent.

## ex8 complet sur le GPU

5 pas, maillage d'origine, un processus (cœurs 4 à 7), tout par noyaux,
lois cuda, hypre sur l'hôte (thermique par GMRES et BoomerAMG), mécanique
symétrique définie positive. Construction `build/examples-mumps`, liée à
`install/mfem-mgis-mumps` (branche cudss). Deux passages alternés
(`resume.txt`, commandes dans `mesure-ex8-gpu.sh`). Temps en secondes :

| Mécanique | Total | Couplage | Factorisations (29) |
|---|---|---|---|
| MUMPS sur l'hôte, analyse réutilisée | 9,80 et 9,73 | 9,60 et 9,55 | 5,87 et 5,83 |
| cuDSS, sans réutilisation | 7,12 et 7,14 | 6,93 et 6,95 | 3,47 et 3,47 |
| cuDSS, analyse réutilisée | 4,37 et 4,38 | 4,19 et 4,19 | 0,69 et 0,70 |

- ex8 sur le GPU passe de 11,1 s (MUMPS sans réutilisation, avant ce lot
  et le précédent) à 4,4 s, soit ×2,5.
- Après la relecture (45afaa0653 : solutions contrôlées par leur résidu,
  motif comparé sans les valeurs ; `resume-relecture.txt`, commandes dans
  `mesure-ex8-gpu-relecture.sh`, cœurs 0 à 3) : 4,52 et 4,41 s avec hypre
  sur l'hôte (le contrôle coûte au plus 0,1 s) ; 5,63 et 5,58 s avec hypre
  sur le GPU (thermique plus lente, factorisations 0,63 s, la comparaison
  des motifs ne lit plus que les indices) ; 8,38 et 8,45 s avec hypre sur
  le GPU sans réutilisation. Mêmes champs.
- Mêmes statistiques des champs (température, norme du déplacement) et
  mêmes nombres d'itérations (114) dans les trois cas.
- Les factorisations ne pèsent plus que 16 % du temps.

## ex8, thermique aussi par cuDSS

Avec les deux solveurs directs sur le GPU (`-lsTh CuDSSSolver -lsMc
CuDSSSolver -mra`, matrices thermiques données entières, mécaniques par
leur triangle inférieur ; `resume-thermique.txt`, commandes dans
`mesure-ex8-gpu-thermique.sh`, cœurs 0 à 3), en secondes :

| hypre | Total | Couplage | Factorisations thermiques (35) | Résolutions thermiques (35) | Factorisations mécaniques (29) |
|---|---|---|---|---|---|
| sur l'hôte | 2,78 et 2,77 | 2,60 et 2,59 | 0,22 et 0,21 | 0,036 et 0,034 | 0,68 et 0,68 |
| sur le GPU | 2,58 et 2,73 | 2,40 et 2,53 | 0,20 et 0,21 | 0,036 et 0,036 | 0,62 et 0,62 |

- Les résolutions thermiques par GMRES et BoomerAMG prenaient 1,83 s
  (hypre sur l'hôte) : cuDSS factorise la matrice thermique (9 516
  inconnues) en 6 ms et la résout en 1 ms.
- ex8 sur le GPU passe de 11,1 s ce matin (MUMPS sans réutilisation,
  GMRES et BoomerAMG) à 2,6 s, soit ×4,3, mêmes résultats et mêmes
  itérations (114).
- hypre ne fait plus que porter les matrices et les vecteurs parallèles :
  le garder sur le GPU évite les copies entre l'hôte et le GPU des
  matrices données à cuDSS.

## Limites de cuDSS distribué

Avec plusieurs processus sur un GPU, l'analyse distribuée de cuDSS 0.7.1
échoue (`CUDSS_STATUS_EXECUTION_FAILED`, code 5, ou erreur de
segmentation) pour certaines matrices, que MUMPS résout :

- quand des lignes n'ont que leur coefficient diagonal (sommets isolés du
  graphe de la matrice) : dès 2 telles lignes quand les processus sont
  couplés ; sans couplage, le seuil dépend du motif ;
- pour une matrice tridiagonale couplant les processus, dès 4 lignes par
  processus ;
- pour certaines répartitions inégales des lignes d'une matrice dont
  chaque ligne est couplée à une ligne d'un autre processus (4, 2 et 6
  lignes sur trois processus ; 4, 3 et 5 passent).

Un processus n'est pas touché. Reproducteur autonome :
`reproducers/cudss-distributed`. Les matrices de diffusion, d'appariement
et de couplage du test passent ; les matrices tangentes de mfem-mgis
gardent les coefficients des lignes des degrés de liberté imposés (mis à
zéro, pas retirés). Le test de mfem-mgis saute les matrices à lignes
diagonales pour cuDSS sur plusieurs processus, et déplace une seule ligne
entre processus.

## Tests

`DirectSolverReuseAnalysisTest`, 10 variantes : MUMPS sur 1, 2 et 3
processus et sur le GPU, cuDSS refusé sans GPU, cuDSS sur le GPU sur 1, 2
et 3 processus et avec hypre sur l'hôte sur 1 et 2 processus ; 10/10. Une
matrice définie négative déclarée définie positive donne avec cuDSS une
solution finie et fausse, rejetée par le contrôle du résidu. Mutations :
six détectées (décision collective court-circuitée, motif non comparé,
résidu non contrôlé, motif oublié après une réutilisation, tampons non
recopiés, indices globaux des colonnes hors diagonale non comparés) ; la
décision locale sans réutilisation ne l'est pas, et elle ne change rien
ici : la disposition comparée comprend le nombre global de lignes et la
première ligne, et un processus qui garde la sienne garde des
descripteurs valides. Relecture indépendante Z9 traitée.
