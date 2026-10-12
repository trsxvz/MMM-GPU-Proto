# Exemples ex1, ex2, ex3, ex5, ex6 et ex7 par lots et sur le GPU (12/10/2026)

Six exemples de mfem-mgis-examples qui n'étaient pas encore suivis par
le prototype (ex4 et ex8 l'étaient déjà), avec leurs lois générées par
l'interface `generic-parallel` :

| Exemple | Problème | Loi | Vérification de l'exemple |
|---|---|---|---|
| ex1 | traction uniaxiale cyclique d'un cube (sources installées par mfem-mgis) | plasticité à écrouissage linéaire | gradient, force et déformation plastique cumulée de chaque pas contre la référence |
| ex3 | cube à deux couches, six cas de chargement (sources installées par mfem-mgis) | élasticité linéaire | solution contre la solution analytique |
| ex7 | VER de MOX, une inclusion sphérique, conditions périodiques, tétraèdres quadratiques, GMRES et BoomerAMG | Norton à seuil (implicite) | contraintes moyennes par matériau contre la référence (tolérance relative 1e-4), 5 pas |
| ex6 | VER de deux couches, conditions périodiques, gradient de déformation macroscopique imposé, PCG et BoomerAMG | Saint Venant-Kirchhoff | gradient et contrainte de Piola-Kirchhoff en chaque point contre la solution analytique (1e-10) |
| ex2 | ssna303 en déformations planes, grandes déformations (Hencky), MUMPS, sans et avec FBar | plasticité à écrouissage linéaire | effort sur le bord supérieur contre la référence, 2 pas |
| ex5 | plaque de Satoh en déformations planes, température parabolique imposée aux points, problème séquentiel, UMFPack | thermoélasticité | température vue par la loi en chaque point (1e-9 K) et réactions des bords encastrés (1e-8) |

Options ajoutées aux exemples (mfem-mgis-examples 2d67da3 et b602495) :
`-l` (bibliothèque de la loi) pour ex2 et ex5, options de mfem-mgis
(`declareDefaultOptions`, dont `--device`) pour ex6 et ex5. ex1 et ex3,
dont les sources sont des tests de mfem-mgis, ont déjà `-l` ; le
périphérique leur est donné par la variable `MFEM_MGIS_DEVICE`.

Tests (`cas/ex1`, `cas/ex2`, `cas/ex3`, `cas/ex5`, `cas/ex6`, `cas/ex7`, construction
`build/cas-post`, mfem-mgis 9338554cd7) : point par point, tout par
noyaux avec la loi de l'exemple, par lots avec les backends openmp et
stlpar, et sur le GPU (cuda, 2 processus pour ex2, ex6 et ex7, un
processus pour ex1, ex3 et ex5).

| Exemple | CPU | GPU |
|---|---|---|
| ex1 | 4/4 | 1/1 |
| ex3 (six cas) | 24/24 | 6/6 |
| ex7 | 4/4 | 1/1 |
| ex6 | 4/4 | 1/1 |
| ex2 (sans et avec FBar) | 8/8 | 2/2 |
| ex5 | 4/4 | 1/1 |

Sur le GPU, le profil d'ex1 (nsys) montre les noyaux de mfem-mgis (forces
intérieures, matrices de rigidité) et l'intégration de la loi par lots
(`mfront::gpb::gpu`).

## Mesures : ex7 complet (40 pas)

VER de MOX, 656 tétraèdres quadratiques, 2 634 inconnues, 109 itérations
de Newton (149 évaluations du résidu) dans tous les cas. Processus placés
sur les cœurs 4 à 7 (un cœur par processus), un passage (deux sur le
GPU). Durées totales et commandes : `mesures/durees.txt` ; tables de temps
de mfem-mgis : `mesures/*.txt`.

| Variante | Total (s) | Résolution (s) | Jacobienne (s) | Solveur linéaire (s) |
|---|---|---|---|---|
| CPU, 1 processus, point par point | 15,1 | 14,3 | 4,69 | 3,26 |
| CPU, 1 processus, noyaux | 8,0 | 7,3 | 1,22 | 3,18 |
| CPU, 2 processus, point par point | 10,1 | 9,3 | 2,64 | 3,43 |
| CPU, 2 processus, noyaux | 6,3 | 5,5 | 0,60 | 3,39 |
| CPU, 2 processus, lots openmp | 6,1 | 5,5 | 0,60 | 3,38 |
| GPU, 1 processus, lots cuda | 6,5 à 6,8 | 5,6 à 5,7 | 0,44 | 4,96 à 5,07 |

- Les noyaux divisent par 2 le temps sur le CPU (jacobienne par 4).
- Sur le GPU, l'intégration et la jacobienne ne coûtent presque plus
  rien, mais GMRES et BoomerAMG sur le GPU prennent 5,0 s contre 3,2 s
  sur un cœur : le système (2 634 inconnues) est trop petit pour le GPU,
  comme la thermique d'ex8 (`results/ex8-thermique-gpu-2026-10-11`).
- Piège de mesure trouvé en route : avec `scripts/env.sh`
  (`OMP_PROC_BIND=spread`), `mpirun --bind-to none taskset -c ...` met
  tous les processus sur le même cœur (ex7 à 5 pas : 95 s au lieu de
  2,6 s). Les scripts de mesure du prototype
  (`cas/mesures-commun.sh`) ne l'utilisent pas ; ici, `taskset` est mis
  autour de `mpirun`, qui lie un processus par cœur.

## Variables statiques des quantités sur le GPU

La loi d'ex5 ne se compilait pas pour CUDA : sa variable statique
`@StaticVariable temperature Tref` (avec `@UseQt true`) est un membre
constexpr d'un type classe, et `T + dT - Tref` lui lie une référence dans
le code GPU, ce que nvcc refuse (« identifier Tref is undefined in device
code »). Correction dans TFEL (837338523, f1ff18259 et d22416824 après
relecture) : MFront écrit les variables statiques
`TFEL_STATIC_VARIABLE(v)`, macro de `TFELConfig.hxx` qui donne la
variable elle-même dans le code de l'hôte (inchangé) et une copie
initialisée par une expression constante, `(decltype(v){v})`, dans le
code GPU (`reproducers/mfront-static-quantity-cuda`). Seconde relecture
(bf143cf46) : rigidité recalculée par les algorithmes de Runge-Kutta et
constantes de rkCastem écrites avec la macro (les six algorithmes et
rkCastem avec `@UseQt true` échouaient en CUDA), et modificateurs de
Runge-Kutta renvoyés par valeur (une lambda statique gardait le nom de la
première loi traitée par un appel de MFront).

- Loi de test `StrainHardeningCreep3` de TFEL (0f58645c7) : une vitesse de
  déformation statique nulle utilisée à l'exécution dans la loi
  d'écoulement. Avec l'ancien MFront, la configuration CUDA échoue ; avec
  la correction, comparaisons openmp-host-gcc 10/10 et cuda-nvcc 11/11.
- Cinq lois de test de TFEL à variables statiques, `StrainHardeningCreep3`
  et la loi d'ex5 se compilent pour l'interface `generic` et pour CUDA.
  Les dix lois écrites par la relecture pour les cas limites (adresse,
  `std::ref`, références non constantes, fonction qui renvoie une
  référence) se compilent sur l'hôte comme avant.

## Fichiers

- `ctest-*.txt` : sorties des tests (`ctest-ex5-cpu.txt` avant la
  régénération des lois par le nouveau MFront, mêmes résultats après).
- `tfel-*-StrainHardeningCreep3.txt` : comparaisons de TFEL
  (`scratch/dev-test.sh`).
