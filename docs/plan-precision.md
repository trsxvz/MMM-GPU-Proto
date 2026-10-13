# Plan : précision float ou double unifiée sur toute la chaîne

Statut : plan du 13/10/2026 (tâche 9b), à valider avant de coder. Une
relecture indépendante (même jour) a vérifié chaque constat contre les
sources et en a corrigé ou complété plusieurs : elle est intégrée ici. Les
décisions à prendre sont regroupées à la fin (D1 à D9), avec une
recommandation pour chacune.

## Objectif

Un seul type réel pour tout un calcul, dicté par le MFEM installé
(`mfem::real_t`) : hypre, MGIS, les lois générées par MFront, mfem-mgis
et les exemples suivent. Un mélange (par exemple des lois en double et un
MFEM en float) doit être refusé avec un message clair, à la configuration
ou au chargement, jamais donner des résultats faux.

Intérêt : les cartes GeForce, dont celle de la machine (RTX 5080),
calculent en double beaucoup plus lentement qu'en simple (1/64 d'après les
caractéristiques publiées) ; les phases limitées par la mémoire lisent
moins d'octets (deux fois moins pour les valeurs, 1,5 fois moins pour une
matrice CSR à index 32 bits).

Hors du périmètre : précision mixte (lois en double et éléments finis en
float, ou l'inverse). Le schéma retenu ne l'empêche pas plus tard
(conversion des tableaux des matériaux).

## État actuel (constaté le 13/10/2026)

MFEM et ses dépendances :

- `mfem::real_t` vaut `float` si MFEM est compilé avec `MFEM_USE_SINGLE`
  (`config/config.hpp`), `double` sinon. `config.mk` l'écrit aussi, et
  `FindMFEM.cmake` de mfem-mgis en lit déjà toutes les variables.
- Recette spack (MFEM 4.7 et suivants) : variante `precision=single|double`,
  qui impose `hypre precision=single`, `petsc~double` et `mumps+float`.
  MFEM vérifie à la compilation hypre (`hypre.hpp`) et PETSc (`petsc.hpp`),
  et utilise `smumps` en simple. cuDSS accepte le simple.
- MFEM refuse en simple : SuiteSparse (UMFPACK), AmgX, GSLIB, SUNDIALS,
  SuperLU, STRUMPACK, Ginkgo, SLEPc, libCEED et quelques autres. Ce refus
  n'existe que dans son `makefile` (`PKGS_NO_SINGLE`), la recette spack n'a
  pas de `conflicts` : l'échec arrive tard, à la construction de MFEM. Le
  spec de l'environnement mmm-gpu a 8 de ces variantes (`+libceed +ginkgo
  +sundials +superlu-dist +strumpack +amgx +suite-sparse +gslib`).

MGIS :

- `mgis::real` vaut `mgis_real`, défini par la macro `MGIS_REAL_TYPE`
  (`MGIS/Config-c.h`), `double` par défaut. Aucune option de CMake, aucun
  en-tête généré : une application qui définirait la macro autrement que
  la bibliothèque compilée aurait des structures incompatibles, sans
  erreur.
- MGIS ne compile pas en float en l'état : la clé de `rdt` de
  `BatchDataView.hxx` (16 sources sur 47), `CUDAFieldOperations.hxx` et
  `HIPFieldOperations.hxx`. Avec `BatchDataView.hxx` corrigé, tout le cœur
  sur CPU compile.
- Les structures C de MGIS et de MFront sont identiques membre à membre,
  structures par lots du prototype comprises.
- Les paramètres des lois (`setParameter`, valeurs par défaut) sont en
  `double`, les bornes en `long double` : c'est l'interface binaire des
  bibliothèques de MFront, indépendante du type des tableaux.
- Liaisons Python, Fortran, Julia et FEniCS : près de 300 `double` (ou
  `c_double`, `float64`) en dur. Les liaisons C suivent `mgis_real`.

MFront et MTest :

- `mfront_gb_real` est `double` en dur (`MFront/GenericBehaviour/Types.h`),
  `mfront_gmp_real` aussi pour les propriétés matériau. Les lois sont des
  modèles paramétrés par leur type numérique, et le code généré par
  l'interface `generic` n'utilise que `mfront::gb::real`.
- L'interface `generic-parallel` code le facteur `rdt` sur 64 bits
  (`memcpy` de 8 octets dans `GenericParallelBehaviour/BatchData.h`,
  clé `toOrderedKey(const double)` dans `Integrate.hxx`), et MGIS partage
  cette clé : en float, le facteur rendu est faux, sans erreur (nvcc le
  signale par `-Wstringop-overread`). Changer le `typedef` ne suffit donc
  pas pour cette interface.
- Le corps des propriétés matériau (`generic`, `generic-parallel`) et des
  modèles est calculé en `double` en dur : juste, mais sans gain sur GPU ;
  l'écriture des sorties des modèles ne compile probablement pas en float
  (non essayé).
- Aucun symbole n'indique le type d'une bibliothèque : une loi compilée en
  double chargée par un MGIS en float lirait des tableaux de travers.
- MTest calcule en double (`MTest/Types.hxx`) et charge une loi sans
  vérifier son type : il lirait une loi en float de travers. Le test
  `u3si2_swelling` d'ex8 charge la bibliothèque des lois de mfem-mgis.
- Critères de convergence des lois implicites : 1e-8 par défaut, mais les
  lois des exemples les fixent à 1e-16 (ex7), 1e-14 (NortonCr), 1e-12
  (MonoCristal_UO2), 1e-11 (Mono_UO2_CosH_Jaco3). En float (précision
  relative 1,2e-7), ils sont hors d'atteinte. `numerical_jacobian_epsilon`
  vaut 0,1 × epsilon.
- `maximal_time_step_scaling_factor` vaut `DBL_MAX` : converti en float, il
  donne +inf en pratique (comportement indéfini selon la norme), sans effet
  sur les `std::min` où il sert.

mfem-mgis :

- `FindMFEM.cmake` exige SuiteSparse (erreur fatale sans
  `MFEM_USE_SUITESPARSE`) : mfem-mgis ne se configure pas avec un MFEM en
  float.
- `mfem_mgis::real` vaut `mgis::real`. Un mélange MGIS double et MFEM
  float ne compile déjà pas (`std::span<real>` sur `real_t*` dans
  `MaterialArrays.cxx`), mais sans message clair.
- `ParameterVariant` (`Parameter.hxx`) contient `real`, pas `double` : en
  float, `{"Tolerance", 1e-10}` ne compile plus, ni `is<double>` et
  `get<double>` (une vingtaine d'appels).
- AmgX : l'option `enable-amgx` de mfem-mgis est indépendante de MFEM ;
  en float, les enveloppes (mode `dDDI`) liraient des float comme des
  double, sans erreur.
- MPI : `MPI_DOUBLE` sur des données en `real` dans
  `ComputeResultantForceOnBoundary.cxx` (débordement silencieux) et
  `MeanThermodynamicForces.ixx` (ne compile pas). Aussi :
  `PartialQuadratureFunction.cxx` (`double Eval(...) override`, ne
  compile pas), `numeric_limits<double>::max()` rangé dans un `real`
  (`PeriodicNonLinearEvolutionProblem.cxx`). `mfem_mgis::mpi_type<real>`
  (`MPI.hxx`) existe déjà.
- UMFPACK et GSLIB sont conditionnels dans la bibliothèque et ses tests,
  mais ex2 choisit UMFPACK sans garde, et ex5 est entièrement sauté sans
  SuiteSparse.
- Tolérances en dur : une vingtaine dans les exemples et mm-opera-hpc, 46
  fichiers de tests de mfem-mgis. Une boucle de mm-opera-hpc
  (`MacroscropicElasticMaterialProperties.cxx`, `while (... > 1e-10)` en
  `real`) pourrait ne jamais finir en float (non vérifié).
- Lois construites hors de `mfem_mgis_add_behaviours_library` : la
  bibliothèque `BehaviourTest` des tests de mfem-mgis, les lois des tests
  de MGIS.

Essai sans rien modifier : code généré de lois compilé avec
`mfront_gb_real` en float (copie de `Types.h` placée avant celle de TFEL).

- Exemples et mm-opera-hpc, 14 lois : 9 compilent avec les mêmes
  avertissements qu'en double (Elasticity, IsotropicLinearElasticity,
  IsotropicLinearHardeningPlasticity, IsotropicLinearThermoElasticity,
  MonoCristal_UO2, SaintVenantKirchhoffElasticity,
  StationaryLinearHeatTransfer, U3SI2_IrradiationCreep,
  U3SI2_SolidSwelling). 5 échouent dans leur propre code : un littéral
  double avec un objet en float donne une expression en double, que TFEL
  refuse d'affecter ou de copier dans un objet en float (`2.*mu*...` dans
  NortonViscoplasticityWithThreshold d'ex7 et NortonCr, `273.15` avec
  `auto` dans les CoupledHeatTransfer d'ex8, une initialisation par copie
  dans Mono_UO2_CosH_Jaco3).
- Tests de mfem-mgis, 17 lois : Mazars échoue dans son code,
  MicromorphicDamageI_SpectralSplit dans TFEL lui-même
  (`StensorComputeIsotropicFunctionDerivative.hxx`, une fonction qui rend
  `0.5`).
- nvcc, interface `generic-parallel` : 3 sources sur 4 compilent, la
  quatrième (CoupledHeatTransfer) échoue comme sur le CPU.

Essai élargi (validation, point 2, fait le 13/10/2026) : les 856 fichiers
MFront de TFEL (tests, suite de tests, galerie), de MGIS, de mfem-mgis,
des exemples et de mm-opera-hpc, interface `generic`, sur le CPU, avec
`mfront_gb_real` en float. 99 ne passent pas la génération (DSL ou options
non prévus par l'essai), 46 ne compilent pas non plus en double (fichiers
anciens) : restent 711 fichiers, dont 580 compilent en float et 131
échouent seulement en float (18 %). Causes :

- `std::max` et `std::min` sur un double et un float : 39 (la déduction
  des paramètres échoue, aucune conversion à l'affectation n'y peut
  rien) ;
- construction d'un tenseur à partir d'une expression en double : 29 ;
- affectation ou conversion d'une expression en double vers un objet en
  float : une trentaine ;
- modèles (sorties en `map<double>`) : 9 ;
- dans les en-têtes de TFEL eux-mêmes : 18 (critères de Cazacu, de
  Gurson-Tvergaard-Needleman, de Rousselier, de Mohr-Coulomb, tenseur
  d'élasticité, nucléation, `StensorComputeIsotropicFunctionDerivative`).

Par origine : TFEL 97 échecs sur 641, MGIS 2 sur 28, mfem-mgis 2 sur 17,
exemples 5 sur 19, mm-opera-hpc 3 sur 6. L'option D3 (a) seule ne
corrigerait donc qu'une partie des lois (affectations et constructions) :
les appels mixtes à `std::max` et `std::min` et le code de TFEL demandent
de toute façon des corrections du code.

## Proposition

### Source de vérité

Le type est celui de MFEM (`MFEM_USE_SINGLE` de `config.mk`). Il n'y a pas
d'option de mfem-mgis pour le choisir : changer de précision veut dire
installer un autre MFEM (variante spack `precision`), qui fixe aussi hypre,
PETSc et MUMPS.

### TFEL (MFront, MTest)

- `Types.h` des lois et des propriétés matériau : le type vaut la macro
  `MFRONT_GB_REAL_TYPE` (`MFRONT_GMP_REAL_TYPE`) si elle est définie,
  `double` sinon. Rien ne change pour qui ne la définit pas.
- Une option de `mfront` (D2) écrit les deux macros en tête des sources
  générées (une bibliothèque peut contenir lois et propriétés, comme
  `libU3SI2.so` d'ex8).
- Chaque bibliothèque exporte la taille de son type réel, par fonction
  (symbole `<fonction>_real_size`, calculé à la compilation), écrit par les
  trois générateurs (`GenericBehaviourSymbolsGenerator`, propriétés
  `generic` et `generic-parallel`). Les versions d'API (`_api_version`,
  `_generic_parallel_api_version`) changent aussi en float : un MGIS
  ancien refuse alors la bibliothèque au lieu de la lire de travers.
- `generic-parallel` : le facteur `rdt` est codé et relu converti en
  double, des deux côtés (MFront et MGIS).
- Code de TFEL lui-même écrit pour tout type numérique
  (`StensorComputeIsotropicFunctionDerivative` et ce que l'essai élargi
  trouvera) ; littéraux double dans le code des lois : D3.
- Propriétés matériau : interface au type choisi, corps en double
  inchangé dans un premier temps. Modèles : refusés en float au
  chargement dans un premier temps (D9).
- MTest : vérifie `<fonction>_real_size` au chargement et refuse une loi
  qui n'est pas en double.

### MGIS

- Option de CMake `MGIS_REAL_TYPE` (`double` par défaut, `float`), écrite
  dans un en-tête installé (`MGIS/RealType.h`, inclus par `Config-c.h`) et
  exportée par la configuration CMake installée (comme `MGIS_HAVE_TFEL`).
  La macro définie par une application devient une erreur si elle
  contredit l'en-tête.
- `BatchDataView.hxx` (clé de `rdt` via double), `CUDAFieldOperations.hxx`
  et `HIPFieldOperations.hxx` corrigés.
- Vérification au chargement dans `LibrariesManager::checkBatchInterface`,
  par où passent les vingt accesseurs (lois, initialisations,
  post-traitements, rotations, fonctions par lots, propriétés) : si le
  symbole `<fonction>_real_size` existe, il doit valoir
  `sizeof(mgis::real)` ; sinon la bibliothèque est ancienne, donc en
  double, refusée si `mgis::real` ne l'est pas. Message : nom de la loi,
  des deux types, et comment reconstruire. Seul trou restant : une
  bibliothèque sans symbole compilée en float par un ancien `mfront` avec
  la macro en drapeau, prise pour double (rare, documenté).
- Paramètres et bornes : interface inchangée en double (D6).
- Liaisons (Python, Fortran, Julia, FEniCS) en float : refusées à la
  configuration dans un premier temps (D5).

### mfem-mgis

- Configuration : `MFEM_USE_SINGLE` de `config.mk`, `MGIS_REAL_TYPE` de la
  configuration de MGIS, erreur claire s'ils diffèrent ; `static_assert`
  dans un source de la bibliothèque (pas dans `Config.hxx`, qui n'inclut
  pas MFEM).
- `FindMFEM.cmake` n'exige plus SuiteSparse ; ce qui l'utilise est déjà
  conditionnel.
- `enable-amgx` refusé avec un MFEM en float (D9 peut ouvrir le mode
  `dFFI` plus tard).
- `ParameterVariant` garde des `double` (D8) ; conversion en `real` à la
  lecture.
- MPI par `mfem_mgis::mpi_type<real>` partout où les données sont en
  `real` ; `PartialQuadratureFunction::Eval` et les
  `numeric_limits<double>` corrigés.
- `mfem_mgis_add_behaviours_library` et la bibliothèque `BehaviourTest`
  passent le type à `mfront` (option de D2) ; la variable
  `<cible>_REAL_TYPE` le donne aux exemples.
- Lois implicites : `epsilon` et `numerical_jacobian_epsilon` bornés en
  float dès le départ (D7).
- Tolérances relatives des solveurs (D7) : une fonction
  `mfem_mgis::getDefaultRelativeTolerance(t)` qui rend `t` en double et
  `max(t, k ε)` en float ; les tolérances absolues ont une unité et restent
  à la charge de l'appelant.

### Exemples et mm-opera-hpc

- Tolérances relatives par la fonction ci-dessus, références comparées
  avec une tolérance déduite du type ; boucle de mm-opera-hpc bornée.
- ex2 : UMFPACK seulement si MFEM l'a.
- Lois corrigées selon D3 ; le `.mtest` d'ex8 lit la bibliothèque en
  double construite pour MTest (D4).

## Validation

1. En double : toutes les suites inchangées (CPU, GPU, exemples,
   mm-opera-hpc), mêmes résultats au bit près. C'est le chemin par défaut,
   il ne doit pas bouger.
2. Essai de compilation élargi avant de coder D3 : toutes les lois des
   tests de MGIS, de mfem-mgis, des exemples, de mm-opera-hpc et de TDLS,
   en float, sur le CPU et avec nvcc.
3. En float, sur le CPU : environnement spack séparé (MFEM
   `precision=single`, hypre, MUMPS ; sans les 8 variantes refusées, sans
   CUDA, compilation courte), MGIS en float ; suites de mfem-mgis et
   exemples, liste des tests sautés et pourquoi, refus des mélanges (loi
   en double chargée en float, MGIS et MFEM différents, MTest et loi en
   float, ancien MGIS et loi en float).
4. En float, sur le GPU : environnement spack séparé avec CUDA (compilation
   longue de MFEM, voir `deps/spack/README.md`), sans AmgX.
5. Mesures sur la poutre : temps et nombre d'itérations float contre
   double, sur le CPU et sur le GPU ; écart de la solution ; convergence
   des lois implicites avec leurs critères bornés.

## Décisions à prendre

- D1. Source de vérité : le MFEM installé, pas d'option propre à
  mfem-mgis. Recommandé (demandé).
- D2. Comment MFront reçoit le type. (a) Macros seules, passées par les
  drapeaux de compilation de chaque backend. (b) Macros, plus une option
  de `mfront` qui les écrit dans les sources générées. Recommandé : (b),
  un seul point d'entrée quel que soit le backend ; le symbole de taille
  protège dans les deux cas.
- D3. Littéraux double dans le code des lois. (a) TFEL accepte
  l'affectation et la construction d'un objet en float à partir d'une
  expression d'un type plus large, avec conversion ; les lois existantes
  compilent telles quelles, mais des calculs restent en double, même sur
  le GPU, sans rien signaler. (b) Réécrire les lois avec des littéraux
  typés (`real(2)`), TFEL reste strict. Recommandé : (b) pour les lois des
  exemples et des tests (mesures honnêtes, avertissements
  `-Wdouble-promotion` et de ptxas pour trouver le reste), et (a) en plus
  pour accepter les lois existantes des utilisateurs. (a) est un
  changement de conception de TFEL, à proposer à ses auteurs avant toute
  publication. D'après l'essai élargi, (a) ne couvre pas les appels mixtes
  à `std::max` et `std::min` (39 lois sur 131) : (b) est nécessaire pour
  une partie des lois de toute façon.
- D4. MTest en float. La fonction CMake construit aussi une bibliothèque
  en double pour les tests MTest (ex8), et MTest refuse une loi en float.
  Recommandé.
- D5. Liaisons de MGIS (Python, Fortran, Julia, FEniCS) en float :
  refusées à la configuration pour commencer. Recommandé.
- D6. Interface des paramètres des lois (double) et des bornes (long
  double) : inchangée. Recommandé : ce sont des données d'entrée,
  converties à l'usage, et la changer casserait toutes les bibliothèques
  existantes.
- D7. Critères : en float, `epsilon` des lois implicites borné à
  `k_l ε` (et `numerical_jacobian_epsilon` en proportion) par
  `setParameter` au chargement, avec un message ; tolérances relatives des
  solveurs bornées à `k ε`. `k` et `k_l` à mesurer (100 au départ).
- D8. Paramètres de mfem-mgis (`ParameterVariant`) : `double`, et non
  `real`, convertis à la lecture. Recommandé : l'API ne dépend pas de la
  précision, et `{"Tolerance", 1e-10}` compile toujours.
- D9. Périmètre de la première version en float : lois seulement ;
  propriétés matériau à l'interface en float mais calculées en double ;
  modèles refusés ; AmgX refusé. Recommandé.

Avant de construire : accord sur les deux environnements spack séparés en
float (CPU, puis GPU). L'environnement principal et ses MFEM et AmgX ne
sont pas touchés.

## Étapes

1. Essai de compilation élargi (validation, point 2) : fait sur le CPU
   (voir l'état actuel) ; reste nvcc, et la liste précise des corrections
   de TFEL et des lois.
2. TFEL : macros, option de `mfront`, symbole de taille, versions d'API,
   clé de `rdt`, code de TFEL pour tout type, D3, MTest ; tests de TFEL
   (une loi de chaque DSL compilée et intégrée en float, l'interface
   `generic-parallel` avec nvcc).
3. MGIS : option de CMake, en-tête et configuration installés,
   `BatchDataView.hxx` et opérations sur les champs, vérification au
   chargement, refus des liaisons ; tests de MGIS en float et en double,
   tests des refus.
4. mfem-mgis : configuration, SuiteSparse facultatif, AmgX refusé,
   `ParameterVariant`, MPI, critères et tolérances, fonction CMake et
   `BehaviourTest` ; double inchangé, float sur le CPU.
5. Exemples et mm-opera-hpc ; puis float sur le GPU.
6. Mesures, documentation (guide de mfem-mgis, pages de MGIS et de
   MFront).
