# hypre sur l'hôte avec le périphérique sur le GPU (12/10/2026)

Option `--hypre-on-host` de `mfem_mgis::initialize` (ou variable
`MFEM_MGIS_HYPRE_ON_HOST`, mfem-mgis 279b5583d4 et 3ad387b002, branche
hypre-hote) : les
lois, les noyaux de mfem-mgis et les calculs de MFEM restent sur le GPU,
hypre calcule sur l'hôte (MFEM lui donne ses données sur l'hôte). Pour les
petits systèmes linéaires, dont les résolutions sur le GPU sont dominées
par les lancements de noyaux.

Construction `build/examples-hypre` (mfem-mgis-examples, ex7 et ex8),
liée à `install/mfem-mgis-hypre`. GPU, un processus lancé par
`mpirun -n 1` sur les cœurs 4 à 7, tout par noyaux, lois cuda, deux
passages par variante. Mêmes itérations de Newton et mêmes résultats dans
les deux cas.

## ex7 complet (40 pas, 2 634 inconnues, GMRES et BoomerAMG)

| hypre | Total (s) | Résolution (s) | Solveur linéaire (s) | Jacobienne (s) |
|---|---|---|---|---|
| sur le GPU | 6,73 et 6,65 | 5,61 et 5,76 | 4,94 et 5,08 | 0,44 et 0,45 |
| sur l'hôte | 4,32 et 4,57 | 3,46 et 3,47 | 3,12 et 3,13 | 0,11 et 0,11 |

Sur le CPU (`results/examples-gpu-2026-10-12`) : 8,0 s sur un cœur et
6,3 s sur deux, par noyaux. Avec hypre sur l'hôte, le GPU passe devant.
La jacobienne est aussi plus rapide : la matrice globale est assemblée par
les noyaux de mfem-mgis au lieu du produit triple de hypre sur le GPU.

## ex8 (5 pas, thermique de 9 516 inconnues par GMRES et BoomerAMG, mécanique par MUMPS symétrique)

| hypre | Couplage (s) | Thermique (s) | Résolutions thermiques (s) | Mécanique (s) |
|---|---|---|---|---|
| sur le GPU | 12,70 et 12,46 | 4,01 et 3,94 | 3,18 et 3,13 | 8,65 et 8,49 |
| sur l'hôte | 11,76 et 11,40 | 2,89 et 2,71 | 1,97 et 1,85 | 8,84 et 8,65 |

La thermique gagne 1,2 s : autant que MUMPS pour la thermique
(`results/ex8-thermal-gpu-2026-10-11`, 11,5 s), sans changer de
solveur. La mécanique (MUMPS sur l'hôte) ne change pas.

## Tests

`ctest.txt` (mfem-mgis 3ad387b002, après relecture) : `DeviceOptionTest`
(où calcule hypre, option, variable, annulation, refus de
`--hypre-device-pool-size`, périphérique configuré par `MFEM_DEVICE` ou
par le programme), `HypreILUTestP-cuda-hypre-on-host`,
`BatchIntegrationTestP-cuda` et `BatchIntegrationTestP-cuda-hypre-on-host`
(63 s et 32 s), 18/18. Les deux réglages de `initialize` sont
nécessaires : sans le drapeau de MFEM (posé avant l'initialisation de
hypre), le test où le programme crée le périphérique échoue ; sans la
politique hôte (posée après le périphérique), quatre tests échouent,
hypre compilé pour le GPU calculant sur le GPU par défaut. Les mesures
ont été faites avec 279b5583d4, dont le comportement est le même pour un
périphérique créé par `initialize`.
