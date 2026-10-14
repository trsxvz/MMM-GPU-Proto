# Plan : précision float ou double unifiée sur toute la chaîne

Statut : plan du 13/10/2026 (tâche 9b). Une relecture indépendante (même
jour) a vérifié chaque constat contre les sources et en a corrigé ou
complété plusieurs : elle est intégrée ici. Les décisions D1 à D9 ont été
prises le 13/10 au soir, toutes selon la recommandation (D3 : (b) et (a)),
ainsi que les deux environnements spack séparés en float. La mise en œuvre
est décrite à la fin (section « Mise en œuvre »).

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
  ancien refuse alors la bibliothèque au lieu de la lire de travers. Les
  propriétés matériau `generic` n'exportent pas de version d'API : un
  appelant ancien (MGIS ou MTest d'avant ce travail) lirait une telle
  propriété en float de travers (trou documenté).
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

## Mise en œuvre

Fait le 13/10/2026 :

- TFEL : option de DSL `real_type` (`double` par défaut, ou `float`),
  macros `MFRONT_GB_REAL_TYPE` et `MFRONT_GMP_REAL_TYPE` écrites en tête des
  sources des interfaces `generic` et `generic-parallel`, symbole
  `<f>_real_size`, versions d'API décalées de 1000 dans un autre type que
  `double`, clé de `rdt` construite en `double`. Les autres interfaces et
  les modèles refusent `float`. MTest refuse une loi qui n'est pas en
  `double`. `IsAssignableTo` accepte un flottant vers un flottant plus
  étroit (D3 a). En double, le code généré ne change que par le symbole de
  taille.
- MGIS : option `MGIS_REAL_TYPE`, en-tête généré `MGIS/RealType.h`, type
  contrôlé au chargement de chaque point d'entrée, versions attendues
  décalées comme dans MFront. En `float`, le critère `epsilon` des lois
  implicites était relevé à 10 ε (1,2e-6) au chargement (première version
  de D7, remplacée par une borne de MFront à ε/10, voir plus bas). Tests :
  51/51 en double, 47/47 en float.
- Environnement `deps/spack/mmm-float-cpu` : MFEM 4.10
  `precision=single`, hypre et MUMPS en simple.

Ce que les essais ont appris :

- Le facteur du critère se mesure : avec 100 ε (1,2e-5), un incrément de
  déformation de 5e-5 de la loi Norton des tests est jugé convergé avant
  la première itération, avec une contrainte nulle et un succès. 10 ε donne
  le bon résultat sur ce test ; l'erreur sur la déformation viscoplastique
  reste alors bornée par le critère (2,2e-6 au début du fluage, puis moins
  de 1e-9). Le polycristal de mm-opera-hpc a montré ensuite que 10 ε reste
  trop lâche (voir plus bas) : la borne retenue est ε/10.
- Une loi écrite dans des unités où ses coefficients sortent de la plage
  du `float` ne peut pas tourner en `float` : la loi Norton des tests de
  MGIS (`A = 8e-67`, en Pa) donnait `A = 0` et un débordement de
  `seq^n`. Elle a été réécrite `A (seq/K)^n` avec `A = 1` et
  `K = A^(-1/n)`, ce qui ne change le résultat en double qu'à l'arrondi
  près. Les lois des utilisateurs devront souvent être écrites en MPa.
- Les tolérances des tests écrites pour le double (1e-14) n'ont pas de
  sens en `float` : elles sont bornées à 10 ε, ou données explicitement
  pour `float` quand l'erreur dépend du calcul.

Suite (13/10/2026, soir) :

- D7 retenue avec des valeurs mesurées : le critère donné par `@Epsilon`
  est borné à ε/10 (1,2e-8) par MFront à la génération, avec un
  avertissement, et `numerical_jacobian_epsilon` à ε (1,2e-7), y compris sa
  valeur par défaut (le dixième du critère, relevée sans avertissement).
  Le critère par défaut, sans `@Epsilon` (1e-8), proche de la borne, est
  gardé. Les avertissements écrits pour le double (« 1e-14 est un bon
  choix ») ne sont plus donnés en float ; tests de mfront-query sur les
  valeurs bornées (TFEL a84f8012e). Une première version bornait au
  chargement dans MGIS, par `setParameter` : elle écrasait les valeurs
  données à l'exécution (fichier de paramètres, `setParameter`), d'où le
  déplacement dans MFront. Seuls les paramètres `epsilon` et
  `numerical_jacobian_epsilon` sont bornés : un autre critère écrit pour
  le double (critère local d'une loi, borne relative d'une brique) passe
  tel quel.
- Pourquoi ε/10 et non 10 ε (première valeur) : les critères de ces DSL
  sont absolus et les résidus s'expriment en déformations. Avec des
  déformations élastiques de 1e-4, un critère de 1,2e-6 laisse environ
  1 % d'erreur sur les contraintes. Sur le polycristal de mm-opera-hpc
  (loi monocristalline en grandes transformations, `@Epsilon 1e-11`), la
  première itération de Newton de la structure ramenait le résidu à 0,45
  de sa valeur initiale avec 10 ε, contre 0,024 avec 1e-8 ou 1e-10, que
  le solveur local atteint en float ; Newton ne convergeait pas. ε/10 est
  l'ordre du critère par défaut du DSL implicite (1e-8), atteint par des
  résidus en déformations jusqu'à environ 0,1. La perturbation du
  jacobien numérique garde sa propre borne : une perturbation plus petite
  que ε se perd dans l'arrondi des inconnues. Tests de MGIS en float :
  49/49 avec ces bornes.
- Tolérances relatives des solveurs de mfem-mgis
  (`boundRelativeTolerance`) bornées à 100 ε (1,2e-5) : les résidus
  relatifs des tests stagnent vers 10 ε (1,2e-6), et une borne de 10 ε
  laissait Newton juste au-dessus. Une tolérance nulle (nombre
  d'itérations fixé) et les tolérances absolues ne sont pas touchées.
- Modèles : seuls ceux du DSL `Model` sont refusés en float ; les modèles
  écrits avec un DSL de lois (`RungeKuttaModelDSL`) se génèrent en float.
- mfem-mgis : configuration (type de MFEM, refus d'un MGIS d'un autre type
  et d'AmgX en float), paramètres stockés en `double`, types MPI, tests
  compilés en float ; double inchangé (216/216 tests sur le CPU ; les tests
  sur le GPU en double sont repassés à la fin, voir plus bas). Tests adaptés au float
  (`tests/Tolerance.hxx` : tolérance du double bornée à 10 ε, ou donnée
  explicitement pour le float, d'après l'erreur mesurée, avec une marge
  d'environ dix) : 204/204 en float sur le CPU, avec les lois régénérées
  par le mfront courant. Douze
  tests restent en double seulement : ParaSails (NaN dès la première
  itération : hypre en simple appelle `dpotrf_` et `dpotrs_` de LAPACK en
  double ; libHYPRE.so en simple importe sept routines LAPACK en double,
  défaut de hypre), UMFPACK, GSLIB, Umpire, cuDSS et CUDA absents
  du MFEM en float. Défaut trouvé dans les tests : quatre réductions MPI
  en `MPI_DOUBLE` sur des `real` (« stack smashing » en parallèle). Les
  sources de mfront sont régénérées quand mfront change (dépendance des
  commandes). `Tolerance.hxx` est installé avec ex1 et ex3, dont les
  sources sont des tests.
- Plancher du résidu en grandes transformations : le gradient de la
  transformation, proche de l'identité, est arrondi à environ 1e-7 en
  float, donc les contraintes ne sont connues qu'à environ ε E (2e4 Pa
  pour l'UO2). Le résidu du polycristal et du cermet de mm-opera-hpc
  stagne vers 6e3 quel que soit le chargement : une tolérance relative au
  résidu initial d'un pas ne peut pas l'atteindre. mm-opera-hpc (9b31d94)
  prend en float une tolérance absolue de dix fois ce plancher, gardée
  quand le critère des itérations de point fixe change ; types MPI de
  `real`, tolérance de la boucle de Kröner `max(1e-10, 10 ε)` et au plus
  1000 itérations, ligne la plus proche choisie par les comparaisons, à
  1e-6 près sur le temps. Double : 6/6, mêmes fichiers au bit près. Float : 6/6, contraintes
  axiales à 0,06 % à 0,15 % du double et à 0,005 % à 0,23 % de Cast3M.
  Le seuil du point fixe sur les contraintes latérales (1e4 Pa) vaut
  5e4 Pa en float depuis 2e4881c : les composantes latérales de F, juste
  sous 1, ne bougent que d'un ulp (6e-8), soit environ 1e4 Pa de
  contrainte. Avec la loi du cermet compilée par gcc 14 (CI simulée), le
  point fixe oscillait entre 1,5e4 et 1,7e4 Pa, le pas de 0,2 s était
  coupé et la contrainte à 0,2 s s'écartait de 2,7 % du double, effet du
  pas de temps (0,27 % du double aux mêmes pas). Avec 5e4 Pa : 0,24 % du
  double avec gcc 13, 0,25 % avec gcc 14 (results/cermet-float-2026-10-14).
- Exemples : 21/21 en float (satoh, qui demande SuiteSparse, et les
  variantes PETSc ne sont pas déclarées) et 24/24 en double contre les
  sources actuelles de mfem-mgis. Les tolérances de Newton et des comparaisons
  restent celles du double, et sont données explicitement pour le float
  d'après les résidus et les erreurs mesurés : le plancher du résidu suit
  l'échelle des efforts (environ 1e-5 de celle-ci), si bien qu'une
  tolérance relative au résidu initial d'un pas, petit en fin de
  chargement ou d'itérations de couplage, ne peut pas être atteinte : une
  tolérance absolue à environ dix fois ce plancher est nécessaire. ex8
  génère aussi ses lois en double pour MTest (D4, dans son propre
  `CMakeLists.txt`).
- MTest contrôle aussi le type des propriétés matériau et des évolutions
  `generic`.
- Deuxième relecture indépendante (scratch/plan-precision/relecture-9b-2.md,
  aucun défaut grave, 3 moyens, 24 remarques, 5 mutations) traitée :
  seuils de TwoLayerCube par solveur (un module faux de 2 % échappait à deux
  cas sur cinq), mesures et intégrales sommées en double dans mfem-mgis et
  comparées au volume calculé par MFEM (une réduction MPI de mauvais type
  dans `computeMeasure` passait inaperçue, avec des mesures de 1e37),
  marges des tolérances remesurées, contrôle du niveau de remplissage
  d'ILU par le résidu après un nombre fixe d'itérations, sources de mfront
  dépendant aussi de libTFELMFront, option de type pour la bibliothèque
  Kokkos des tests et les bancs d'essai.

- Mesures sur le CPU (results/precision-poutre-2026-10-14-0514, poutre
  r4, 111 843 inconnues, deux processus, environnements spack jumeaux
  `mmm-double-cpu` et `mmm-float-cpu`) : à critère de Newton égal (1e-4),
  le float est 1,11 fois plus rapide (20,7 s contre 23,1 s), avec les
  mêmes itérations ; le gain vient de la résolution linéaire (1,19 fois),
  les calculs élément par élément ne gagnant presque rien. La norme des
  déplacements s'écarte de 8e-5 de celle du double convergé. Le critère du
  cas (1e-8, borné à 1,2e-5) n'est pas atteint à r4 : le plancher du résidu
  relatif croît avec la taille du problème (1,3e-5 à r4, sous 1e-5 à r2) ;
  en float, la tolérance de Newton doit être choisie au-dessus.

- MGIS sur le GPU en float (85b1a2a) : les tests CUDA passent (51/51 avec
  CUDA, en double comme en float). Écarts mesurés entre l'hôte et le GPU
  (contractions des opérations flottantes) : 4,8e-7 sur les contraintes,
  2,5e-5 sur les tangentes (relatifs). Limite du float trouvée : la loi de
  Norton des tests (exposant 8,2) n'est plus intégrable vers 3e9 Pa, le
  résidu du schéma implicite ne pouvant pas atteindre son critère ; le
  cycle de FieldOperationsCUDATest fait trois pas en float au lieu de
  quatre.

- mfem-mgis en float sur le GPU (environnement spack `mmm-float-gpu`) :
  suite 270/270 après adaptation de quatre tests (tolérances écrites pour
  le double, écarts plus grands avec un MFEM compilé par nvcc, solveur
  séquentiel non reproductible au bit près sur le GPU sans solveur
  direct). Régularisation d'Anderson relevée à 10 ε en float.
- Poutre sur le GPU (results/precision-poutre-gpu-2026-10-14-0604) : le
  float est 1,31 fois plus rapide quand il converge, mais un calcul sur
  trois diverge, et les calculs convergés s'écartent du double de 3e-4.
- Diagnostic (results/diagnostic-float-gpu-2026-10-14) : nvcc contracte
  par défaut les multiplications et les additions en
  multiplications-additions fusionnées (FMA). Dans les noyaux de raideur, en float, ces
  contractions faisaient diverger le Newton, ou converger vers des
  solutions moins précises, probablement par une dérive le long des modes
  mous de la poutre (flexion). Les tangentes par noyaux et élément par
  élément ne diffèrent pourtant que de 2e-7 : le float travaille ici près
  de sa limite (direction de Newton connue à 16 à 32 % près au pas 8).
- Correction dans mfem-mgis : en float, le code GPU des noyaux de raideur
  (standard, FBar, Faltus et al.) est compilé sans contraction (fonction
  CMake `mfem_mgis_disable_device_floating_point_contraction`), au même
  coût d'assemblage. En double, les contractions sont gardées (assemblage
  13 à 16 % plus lent sans elles, convergence inchangée). Après correction, 31
  calculs sur 31 convergent, à 4,2e-5 au plus du double ; test
  reproducteur `poutre-gpu-r4` dans cas/poutre.
- Mesure refaite après correction
  (results/precision-poutre-gpu-2026-10-14-0843) : float 1,37 fois plus
  rapide que le double, mêmes itérations, à 3,4e-5 au plus du double.

- Polycristal en float (results/precision-polycristal-2026-10-14) : sur
  un pas, la contrainte axiale s'écarte du double de 6,1e-3 à 5e-5 de
  déformation et de 1,0e-3 à 2e-4 ; sur un pas plus grand (1 s), la loi
  échoue et le calcul s'arrête. Cause probable : F, stocké en float près
  de l'identité, ne bouge que par pas d'un ulp (1e-3 d'une déformation de
  5e-5, quantum de contrainte d'environ 1,5e4 Pa, du même ordre que le
  seuil du point fixe sur la contrainte). La cancellation de
  `cosh(x) - 1` dans la loi n'est pas en cause (forme en `expm1` essayée).

Reste : le plancher du résidu du float croît avec la taille du problème
(1,3e-5 à r4 sur le CPU) ; un critère de Newton adapté à la précision
serait utile sur les grands maillages. Pour les lois en grandes
transformations en float, une interface en H = F - I (décision à
prendre) éviterait la quantification de F.
