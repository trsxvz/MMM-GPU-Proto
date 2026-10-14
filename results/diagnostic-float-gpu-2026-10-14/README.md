# Divergence du float sur le GPU, poutre (14/10/2026)

Suite de results/precision-beam-gpu-2026-10-14-0604 : sur le GPU, en
float, au critère de Newton 1e-4, un calcul de la poutre sur trois divergeait
et les calculs convergés s'écartaient du double de 3e-4. Poutre r4 (111 843
inconnues), plasticité J2, PCG de hypre et BoomerAMG (tolérance 1e-4), un
processus sur le GPU, GPU NVIDIA RTX 5080.

## Conclusion

Cause : nvcc contracte par défaut les multiplications et les additions en
multiplications-additions fusionnées (FMA, `-fmad=true`). Dans les noyaux
des matrices de raideur, en float, ces contractions faisaient diverger le
Newton, ou converger vers des solutions moins précises. Le plus probable :
le Newton dérivait le long des modes mous de la poutre (flexion), que le
critère sur le résidu ne voit pas.

Le float travaille ici près de sa limite. Les tangentes par noyaux et
élément par élément ne diffèrent que de 2e-7 du plus grand terme, mais au
pas 8, la direction de Newton n'est connue qu'à 16 à 32 % près (écart entre
PCG et MUMPS). Les arrondis des FMA suffisaient à faire basculer les
calculs.

Correction dans mfem-mgis (efb06c1267) : en float, le code GPU des
noyaux qui calculent des matrices de raideur (standard, FBar, Faltus et
al.) est compilé sans contraction (`-fmad=false` pour nvcc,
`-ffp-contract=off` pour clang), par la fonction CMake
`mfem_mgis_disable_device_floating_point_contraction`. En double, les
contractions sont gardées : elles ne changeaient pas la convergence, et
l'assemblage de la poutre était 13 à 16 % plus lent sans elles. Avec nvcc,
le code compilé pour l'hôte est inchangé.

Après correction, 31 calculs sur 31 convergent (28 en assemblage atomique,
n = 8 à 14 pas, et 3 calculs identiques en assemblage ordonné), à 4,2e-5 au
plus du double. En float, l'assemblage coûte autant qu'avant (83 à 85 ms
pour 16 assemblages, avant comme après). Mesure de performance refaite :
results/precision-beam-gpu-2026-10-14-0843 (float 1,37 fois plus rapide
que le double, 3 calculs sur 3 convergés, à 3,4e-5 au plus du double). Test
reproducteur : `poutre-gpu-r4` dans cases/beam (assemblage ordonné, n = 10,
norme comparée au double à 5e-5 près).

## Étapes

1. Localisation (`scratch`, assemblage ordonné, n = 10). Avec PCG et
   BoomerAMG, seule la tangente élément par élément converge. Elle est
   alors calculée et assemblée sur l'hôte par la boucle de MFEM (425 ms par
   assemblage contre 5,3 ms par noyaux), les lois, les gradients et les
   forces intérieures restant sur le GPU. Gradients et forces intérieures
   par noyaux sont hors de cause. Toutes les tentatives par noyaux en
   assemblage ordonné à n = 10 divergent (PCG de hypre, gradient conjugué
   de MFEM, tolérance linéaire 1e-6). Un calcul avec MUMPS (sur le CPU)
   converge.
2. Tangentes (`tangents.txt`). Assemblées par noyaux (K) et élément par
   élément (E) au même état, en fin des pas 1, 5, 7 et 8 (le pas 8 étant
   non convergé, après 50 itérations) : écart maximal 2e-7 du plus grand
   terme de la matrice, écart de Frobenius 4e-8, soit l'arrondi du float.
   Même nombre de termes (4 702 185). Même asymétrie à partir du pas 5
   (5e-8 du terme maximal ; au pas 1, 0 pour K et 1,3e-8 pour E).
   L'assemblage ne modifie pas le résidu (identique au bit,
   `directions-pas8.txt`).
3. Directions de Newton (`directions-pas8.txt`), à l'état de la première
   itération du pas 8 de la trajectoire par noyaux (résidu 5,70e5), sur des
   copies des matrices migrées sur le GPU (les matrices d'origine donnent
   les mêmes rapports à 1 % près). Un seul état mesuré.

   | Matrice | Solveur | Itérations | Résidu linéaire vrai | Nouveau résidu |
   |---|---|---|---|---|
   | K, float | MUMPS | | 5,4e-3 | 1,15e5 (0,20) |
   | K, float | PCG 1e-4 | 65 | 7,2e-3 | 1,64e5 (0,29) |
   | K, float | PCG 1e-6 | 102 | 7,8e-3 | 1,64e5 (0,29) |
   | E, float | MUMPS | | 3,9e-3 | 7,99e4 (0,14) |
   | E, float | PCG 1e-4 | 62 | 4,1e-3 | 9,45e4 (0,17) |
   | K, double (sa trajectoire) | MUMPS | | 8,9e-14 | 1,69e4 (0,030) |
   | K, double (sa trajectoire) | PCG 1e-4 | 64 | 1,1e-4 | 1,69e4 (0,030) |

   En float, le résidu linéaire vrai plafonne entre 4e-3 et 9e-3, quels que
   soient le solveur et la tolérance. Les directions de PCG et de MUMPS
   diffèrent de 32 % avec K et de 16 % avec E, contre 6e-5 à 4e-4 en double
   (copies des matrices).
   Un pas de Newton réduit le résidu de 3,5 à 7 fois, contre 34 fois en
   double. Le critère de PCG sur la norme 2 du résidu ne change rien. Dans
   les calculs qui divergent, le résidu relatif atteint son minimum à la
   deuxième itération (1,0e-4 à 1,5e-4, juste au-dessus du critère), puis
   croît de 1,5 à 2 fois par itération jusqu'à 3e-3.
4. Sensibilité (`gaps.txt`), n = 8 à 14 pas, critère 1e-4, assemblage
   atomique (non déterministe, un calcul par n), écart de la norme des
   déplacements au double sur le GPU au même nombre de pas :

   | Configuration | Convergés | Écart maximal |
   |---|---|---|
   | float, GPU, tangente par noyaux | 4 sur 7 | 2,8e-4 |
   | float, tangente élément par élément (hôte), lois sur le GPU | 7 sur 7 | 4,9e-5 |
   | float, CPU, noyaux (deux processus) | 7 sur 7 | 7,6e-5 |
   | float, GPU, noyaux, critère 3e-4 | 7 sur 7 | 2,4e-4 |
   | float, GPU, noyaux, critère 1e-3 | 7 sur 7 | 5,3e-4 |

   Sur les calculs conservés du jour avant correction, au critère 1e-4, en
   assemblage atomique (n = 8 à 14, mesure 0604, profilage), 5 calculs sur
   11 divergent (4 sur 8 dans une première série).
5. Sans contraction. mfem-mgis entièrement compilé avec `-fmad=false`
   (`gaps-without-contraction.txt`) : 7 calculs sur 7 convergent par noyaux,
   à 3,7e-5 au plus du double. Élément par élément, les résultats sont
   identiques à 10 chiffres avec et sans contraction : sur ce cas, seuls
   les noyaux de raideur influent. Au premier état du pas 8 de cette trajectoire
   (`directions-pas8-without-contraction.txt`) : résidu linéaire vrai 5e-4, pas
   de MUMPS qui réduit le résidu de 20 fois, K et E équivalents, et
   translations rigides mieux respectées par K (4,2e-8 du terme maximal
   contre 5,1e-8 avec contraction). Correction ciblée sur les noyaux de
   raideur (`gaps-fix.txt`, première version puis version finale) :
   31 calculs sur 31, dont 3 identiques en assemblage ordonné.

## Méthode

Programme de diagnostic (brouillon, non suivi) : copie de cases/beam avec
deux crochets. Le premier agit après la résolution des pas listés dans
`COMPARE_TANGENT_STEPS`. Il assemble la tangente par noyaux
(`GetGradient`), puis élément par élément (`ParNonlinearForm::GetGradient`,
après `setStiffnessMatricesComputedByKernels(false)` et
`synchronizeOnHostForReading`), et compare les deux matrices terme à terme
sur l'hôte.

Le second arrête le pas `COMPARE_DIRECTIONS_STEP` après une itération de
Newton (`COMPARE_MAX_ITERATIONS=1`), puis le programme s'arrête. À cet état,
il calcule les directions de Newton avec MUMPS et avec PCG et BoomerAMG
configuré comme dans mfem-mgis (options de systèmes posées après
l'opérateur), sur des copies migrées sur le GPU puis sur les matrices
d'origine.

Pour chaque direction, il calcule le résidu linéaire vrai (accumulé en
double sur une copie de la matrice float) et le résidu de Newton obtenu
après intégration des lois. Les résolutions sont faites avant la copie sur
l'hôte, qui déplace la matrice de hypre vers l'hôte. L'assemblage est
ordonné dans tous ces calculs.

Autres essais : sans prédiction, le calcul échoue au pas 1 dans tous les
modes (cas qui exige la prédiction) ; cuDSS n'est pas disponible dans
l'environnement float.
