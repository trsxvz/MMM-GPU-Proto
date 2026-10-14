# ex8 : matrices symétriques données à MUMPS (11/10/2026)

Exemple ex8 de mfem-mgis-examples (plaque thermomécanique), maillage
d'origine, 5 pas, MUMPS pour la mécanique, tout par noyaux. Option
`-msMc` d'ex8 (mfem-mgis-examples 0782d73, puis 03c347e) : symétrie des matrices de la
mécanique données à MUMPS (paramètres `Symmetric` et `PositiveDefinite`
du solveur MUMPS de mfem-mgis). Construction `build/examples-post-ex8`,
liée à `install/mfem-mgis-post` (mémoire de l'hôte standard). Machine
calme.

La tangente de la mécanique d'ex8 est symétrique : mêmes itérations de
Newton (29 résolutions de la mécanique, donc 29 factorisations, et 35 de
la thermique) et mêmes statistiques des champs dans tous les cas.

## CPU, 2 processus, lois `generic` (un passage)

| Symétrie | Total (s) | Factorisation (s) |
|---|---|---|
| `Unsymmetric` (défaut) | 22,05 | 7,49 |
| `SymmetricIndefinite` | 20,54 | 5,99 |
| `SymmetricPositiveDefinite` | 19,85 | 5,31 |

## GPU, 1 processus, lois cuda (deux passages)

| Symétrie | Total (s) | Factorisation (s) | Résolutions, deux physiques (s) |
|---|---|---|---|
| `Unsymmetric` (défaut) | 16,50 et 16,55 | 11,37 et 11,41 | 3,54 et 3,57 |
| `SymmetricPositiveDefinite` | 12,34 et 12,48 | 7,30 et 7,33 | 3,47 et 3,56 |

- Sur le GPU, ex8 passe de 16,5 à 12,4 s. MUMPS (sur l'hôte) reste
  l'essentiel : 7,3 s de factorisations et 0,31 s de résolutions (29) sur
  12,4 s. Les 3,5 s de la colonne des résolutions sont surtout celles de
  la thermique : GMRES et BoomerAMG sur le GPU, 3,15 s pour 35
  résolutions (environ 1,5 s sur 2 processus du CPU), second poste après
  MUMPS.
- Avec les mesures du 11/10 (`results/ex8-behaviours-cuda-2026-10-11`,
  18,5 s) et le correctif de la mémoire de l'hôte
  (`results/host-memory-2026-10-11`, 16,9 s), ex8 tout sur le GPU tombe
  à 12,4 s.
