# Post-traitements par noyaux (11/10/2026)

Bulle de mm-opera-hpc, maillage par défaut (`single_sphere.msh`, 111 580
tétraèdres quadratiques, 467 109 inconnues), un pas, PCG et BoomerAMG,
tout par noyaux pour la mécanique. Post-traitements de la bulle : export
des contraintes aux nœuds (`ParaviewExportIntegrationPointResultsAtNodes`),
contrainte principale et pression hydrostatique aux points, moyenne de la
pression hydrostatique (`computeMeasure`, `computeIntegral`), projection
L2 de la contrainte principale, export ParaView.

Comparaison de `MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=ElementByElement`
(boucles sur les éléments, défaut) et `Kernels` (noyaux de mfem-mgis,
branche post-traitements, 0434dfd1de, mémoire de l'hôte standard).
Exécutable : `build/mm-opera-hpc-post` (mm-opera-hpc 79c49a8, moyenne par
`computeMeasure` et `computeIntegral`). Sorties de la bulle identiques aux
chiffres écrits dans tous les cas.

## GPU (un processus, loi cuda, deux passages)

| Phase (s) | ElementByElement | Kernels |
|---|---|---|
| total | 14,77 et 14,86 | 5,32 et 5,33 |
| résolution | 1,48 et 1,48 | 1,47 et 1,48 |
| post-traitements | 11,32 et 11,40 | 1,86 et 1,87 |
| export aux nœuds | 4,46 et 4,46 | 0,85 et 0,85 |
| moyenne hydrostatique | 1,90 et 1,91 | 0,008 et 0,008 |
| projection L2 | 4,35 et 4,42 | 0,39 et 0,40 |
| dont opérateur | 2,42 et 2,46 | 0,037 et 0,037 |
| dont second membre | 1,42 et 1,47 | 0,002 et 0,002 |
| dont résolution | 0,49 et 0,46 | 0,32 et 0,33 |
| export de la contrainte principale | 0,49 et 0,49 | 0,49 et 0,49 |

## CPU (un cœur, un passage)

| Phase (s) | ElementByElement | Kernels |
|---|---|---|
| total | 33,25 | 25,19 |
| post-traitements | 12,16 | 5,35 |
| export aux nœuds | 4,47 | 1,47 |
| moyenne hydrostatique | 1,89 | 0,33 |
| projection L2 | 5,20 | 2,95 |
| dont opérateur | 2,51 | 1,67 |
| dont second membre | 1,41 | 0,01 |
| dont résolution | 1,25 | 1,07 |

- Sur le GPU, les post-traitements vont 6 fois plus vite et la bulle
  entière 2,8 fois plus vite. Avec le correctif de la mémoire de l'hôte
  (`results/host-memory-2026-10-11`), la bulle sur le GPU est passée de
  69 s à 5,3 s, contre 25 s sur un cœur du CPU.
- La matrice de masse des tétraèdres est calculée par un noyau et
  assemblée par MFEM (`AssemblyLevel::FULL`) : MFEM 4.10 n'assemble pas
  sur le périphérique la masse des simplexes.
- Ce qui reste sur le GPU : la résolution de la projection L2 (gradient
  conjugué préconditionné par ILU de hypre), les exports ParaView (sur
  l'hôte) et le rassemblement des valeurs aux points (sur l'hôte, où sont
  rangées les fonctions de l'utilisateur).
- Version finale après relecture (mfem-mgis affeddab21 : masse calculée
  par le noyau pour tous les éléments, maillage non modifié, décisions
  collectives), même machine calme : GPU 14,86 s élément par élément,
  5,32 s par noyaux (post-traitements 11,41 et 1,87 s ; fichiers
  `bulle-gpu-*-final.txt`), sorties identiques.
- Le chemin élément par élément de la moyenne passe maintenant par
  `computeMeasure` et `computeIntegral` (deux boucles) : 1,9 s contre 1,0 s
  pour l'ancienne boucle de la bulle.
