# cuDSS sur le polycristal (12/10/2026)

Polycristal de mm-opera-hpc (5 grains, loi de monocristal de l'UO2 en
grandes transformations, cellule périodique, traction uniaxiale, 20 553
inconnues, éléments linéaires, 10 pas), tout par noyaux, loi cuda, un
processus sur le GPU (cœurs 0 à 3). mfem-mgis 429f76ded3
(`install/mfem-mgis-mumps`), construction `build/mm-opera-hpc-cudss`
(seul `uniaxial-polycrystal`). Option `--reuse-analysis` et solveurs
directs `MUMPSSolver` et `CuDSSSolver` ajoutés à l'option
`--linear-solver` du polycristal (matrices non symétriques). Deux passages
(commandes dans `measure.sh`, contraintes macroscopiques dans `*.res`).
La charge de la machine valait 3,3 au premier passage et 2 au second (une
autre session déboguait sur les cœurs 12 à 15) ; les deux passages
concordent à 1 % près. Temps totaux, en secondes :

| Méthode | Solveur | Temps |
|---|---|---|
| Newton | HyprePCG et BoomerAMG (hypre sur le GPU) | 11,46 et 11,44 |
| NewtonK0, tolérance d'Eisenstat et Walker | gradient conjugué et BoomerAMG construit sur K0 | 8,78 et 8,76 |
| Newton | cuDSS, analyse réutilisée, hypre sur le GPU | 19,12 et 19,11 |
| Newton | cuDSS, analyse réutilisée, hypre sur l'hôte | 19,30 et 19,26 |

- cuDSS ne paie pas sur ce cas : une factorisation prend 132 ms (12,9 s
  pour 98), contre 57 ms pour une résolution par PCG et BoomerAMG (5,6 s
  pour 98). La réutilisation de l'analyse fonctionne (essai court sur deux
  pas : 134 ms par factorisation avec, 222 ms sans), mais la factorisation
  numérique de cette matrice coûte cher : la cellule périodique relie les
  faces opposées, et la matrice est non symétrique.
- À taille voisine, la matrice mécanique d'ex8 (28 548 inconnues, non
  périodique, symétrique définie positive) se factorise en 24 ms
  (`results/cudss-2026-10-12`).
- Mêmes contraintes macroscopiques avec cuDSS et avec PCG (écart relatif
  2,6e-9). NewtonK0 avec Eisenstat et Walker s'écarte de Newton d'environ
  1e-6 sur la contrainte axiale (quelques pascals sur les contraintes
  transverses, presque nulles).
- La meilleure variante reste NewtonK0 avec la tolérance d'Eisenstat et
  Walker (8,8 s ; 9,6 s le 10/10 avec mfem-mgis 4b72095ee0,
  `results/L5-polycrystal-2026-10-10-0502`).
