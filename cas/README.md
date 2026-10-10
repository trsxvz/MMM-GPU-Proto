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

`mesures-reference.sh` mesure les temps de référence du chemin actuel de
mfem-mgis sur CPU (lot L0) : poutre, ssna303 en 3D (ex4), polycristal
(mm-opera-hpc), thermomécanique (ex8). Résultats dans
`results/reference-L0-<date>/`.
