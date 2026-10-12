# Cas de la campagne de mesures

Cas propres au prototype, qui ne sont ni dans mfem-mgis-examples ni dans
mm-opera-hpc. Ils se construisent contre le mfem-mgis installé par
`build-mfem-mgis.sh` (`./build-mfem-mgis.sh cas`, dossier `build/cas`).

- `poutre/` : poutre encastrée en flexion, plasticité J2 en petites
  déformations, 3D (maillage `beam-tet.mesh` et loi
  `IsotropicLinearHardeningPlasticity` des tests de mfem-mgis). Résolue
  par le Newton de mfem-mgis, avec prédiction élastique. Solveurs
  linéaires : `pcg-amg` (PCG de hypre et BoomerAMG), `mumps`, `umfpack`
  (séquentiel). Options principales : `-r` (raffinements), `-n` (pas),
  `-ls`, `-lt` (tolérance du solveur itératif), `--profiling`, `--device`.
  Le bilan donne le nombre d'inconnues, la norme de la solution finale
  (identique pour tous les solveurs) et les itérations de Newton par pas.
  Les tests comparent cette norme à sa valeur de référence.
- `thermique/` : thermique stationnaire non linéaire dans la même poutre
  (conductivité qui dépend de la température, températures imposées aux
  deux extrémités, loi `StationaryNonLinearHeatTransfer` des tests de
  mfem-mgis). Résolue par le Newton de mfem-mgis, point par point ou par
  noyaux, sur l'hôte ou sur le GPU. La tangente thermique est définie
  négative et les températures imposées ont une diagonale à un : la
  matrice est indéfinie, d'où GMRES au lieu du gradient conjugué.
  Solveurs linéaires : `gmres-amg` (GMRES de hypre et BoomerAMG),
  `gmres-jacobi` (GMRES de MFEM et Jacobi, aussi sans matrice), `mumps`,
  `umfpack`. Mêmes options et même bilan que la poutre, plus `-dT`
  (échauffement imposé) et `-alst` (tolérance adaptative).

- `polycristal/` : lois du polycristal de mm-opera-hpc générées par
  l'interface `generic-parallel` (openmp, stlpar, cuda) ; l'exécutable est
  celui de mm-opera-hpc. Les tests comparent les contraintes
  macroscopiques au calcul point par point.
- `bulle/` : loi de la bulle de mm-opera-hpc (cavité sphérique sous
  pression dans une cellule périodique, élasticité, éléments
  quadratiques), générée par l'interface `generic-parallel` ; l'exécutable
  est `test-bubble` de mm-opera-hpc, sur le maillage de sa CI, chaque test
  dans son dossier. Les tests comparent ses sorties (contrainte
  hydrostatique moyenne, contraintes près des bulles) au calcul point par
  point.
- `cermet/` : lois du cermet de mm-opera-hpc (polycristal de l'UO2 et
  interface métallique de Norton), une bibliothèque pour les deux ; un
  pas sur deux processus, contraintes macroscopiques comparées au calcul
  point par point.
- `ssna303/` : lois de l'exemple ex4 de mfem-mgis-examples (barreau
  entaillé en 3D, plasticité en grandes transformations, déformation de
  Hencky) générées par l'interface `generic-parallel` ; l'exécutable est
  `ssna303_3d_hypre` de mfem-mgis-examples (options `--device`, `-l`,
  `-no-pv`, `-r`). Les tests comparent la force résultante du premier pas
  aux valeurs de référence de l'exemple.
- `ex1/`, `ex2/`, `ex3/`, `ex5/`, `ex6/`, `ex7/`, `ex8/` : autres exemples
  de mfem-mgis-examples (traction cyclique, ssna303 en déformations planes
  avec et sans FBar, cube à deux couches, plaque de Satoh, VER de deux
  couches, VER de MOX, plaque thermomécanique), avec leurs lois générées
  par l'interface `generic-parallel` ; l'exécutable est celui de
  l'exemple (options `-l` et `--device`, ou variable `MFEM_MGIS_DEVICE`
  pour ex1 et ex3). Chaque test garde la vérification de l'exemple
  (référence ou solution analytique) : point par point, par noyaux, par
  lots avec chaque backend et sur le GPU. Les dossiers de construction des
  exemples sont donnés par `MFEM_MGIS_EXAMPLES_EX<n>_DIR`.

Options de la poutre ajoutées par les lots suivants : point fixe et
NewtonK0 (`-nls fixed-point -fpm newton-k0 -k0 -pk0`), tolérance
adaptative (`-alst EisenstatWalker1 -alsti 0.1`), multigrille de K0
(`-els amg|amgx`, `-amgo` options de BoomerAMG, `-amgxo`, `-amgxc`,
`-amgxs` pour AmgX), verbosité du solveur non linéaire (`-vl`), solveur
direct cuDSS sur le GPU (`-ls cudss`, `-els cudss`, demande un mfem-mgis
qui le fournit) et réutilisation de l'analyse des solveurs directs
(`-mra`). L'outil
`poutre-k0` résout l'opérateur élastique K0 seul, mise en place et
résolutions chronométrées à part. Les options de mfem-mgis s'appliquent à
tous les cas : `--device`, `--hypre-device-pool-size` (`-hdps`), et les
variables `MFEM_MGIS_*` des calculs par lots et par noyaux.

Scripts de mesures (résultats dans `results/<lot>-<cas>-<date>/`) :
`mesures-commun.sh` (fonctions communes ; cas poutre, polycristal,
thermique, ssna303, ex8, bulle, cermet ; variable `refinement`), `mesures-L2.sh` à
`mesures-L7-thermique.sh` (un script par lot du plan),
`mesures-echelle-poutre.sh` (passage à l'échelle de la poutre, mémoire du
GPU), `mesures-ssna303.sh` (ssna303 sur 8 cœurs et sur le GPU),
`balayage-amgx-k0.sh` (configurations d'AmgX sur K0),
`profil-transferts.sh` (copies entre l'hôte et le GPU, nsys),
`mesures-memoire-hote.sh` (types de mémoire de MFEM sur le GPU). Outils :
`resume-profilage.py` (temps par phase des sorties `--profiling`),
`comparer-sorties.py` (comparaison de sorties numériques, `--csv` pour
les valeurs séparées par des virgules).

`mesures-reference.sh` mesure les temps de référence du chemin actuel de
mfem-mgis sur CPU (lot L0) : poutre, ssna303 en 3D (ex4), polycristal
(mm-opera-hpc), thermomécanique (ex8). Résultats dans
`results/reference-L0-<date>/`.
