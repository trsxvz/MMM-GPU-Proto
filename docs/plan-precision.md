# Plan : précision float ou double unifiée sur toute la chaîne

Statut : plan du 13/10/2026 (tâche 9b), à valider avant de coder. Les
décisions à prendre sont regroupées à la fin (D1 à D7), avec une
recommandation pour chacune.

## Objectif

Un seul type réel pour tout un calcul, dicté par le MFEM installé
(`mfem::real_t`) : hypre, MGIS, les lois générées par MFront, mfem-mgis
et les exemples suivent. Un mélange (par exemple des lois en double et un
MFEM en float) doit être refusé avec un message clair, à la configuration
ou au chargement, jamais donner des résultats faux.

Intérêt : la carte de la machine (GeForce RTX 5080) calcule en double 64
fois moins vite qu'en simple ; les phases limitées par la mémoire
(assemblage, solveur) lisent deux fois moins d'octets en float.

Hors du périmètre : précision mixte (lois en double et éléments finis en
float, ou l'inverse). Le schéma retenu ne l'empêche pas plus tard
(conversion des tableaux des matériaux).

## État actuel (constaté le 13/10/2026)

MFEM et ses dépendances :

- `mfem::real_t` vaut `float` si MFEM est compilé avec `MFEM_USE_SINGLE`
  (`config/config.hpp`), `double` sinon. `config.mk` l'écrit aussi
  (`MFEM_USE_SINGLE = YES`), et mfem-mgis lit déjà ce fichier pour CUDA et
  HIP.
- Recette spack (MFEM 4.7 et suivants) : variante `precision=single|double`,
  qui impose `hypre precision=single`, `petsc~double` et `mumps+float`.
  MFEM vérifie à la compilation hypre (`hypre.hpp`) et PETSc (`petsc.hpp`),
  et utilise `smumps` en simple.
- MFEM refuse en simple : SuiteSparse (UMFPACK), AmgX, GSLIB, SUNDIALS,
  SuperLU, STRUMPACK, Ginkgo, SLEPc, libCEED et quelques autres
  (`makefile`, `PKGS_NO_SINGLE`). cuDSS accepte le simple.

MGIS :

- `mgis::real` vaut `mgis_real`, défini par la macro `MGIS_REAL_TYPE`
  (`MGIS/Config-c.h`), `double` par défaut. Aucune option de CMake, aucun
  en-tête généré : une application qui définirait la macro autrement que
  la bibliothèque compilée aurait des structures incompatibles, sans
  erreur.
- Les paramètres des lois (`setParameter`, valeurs par défaut) sont en
  `double`, les bornes en `long double` : c'est l'interface binaire des
  bibliothèques de MFront, indépendante du type des tableaux.
- Liaisons Python, Fortran et Julia : environ 300 `double` (ou
  `c_double`, `float64`) en dur.

MFront :

- `mfront_gb_real` est `double` en dur (`MFront/GenericBehaviour/Types.h`),
  `mfront_gmp_real` aussi pour les propriétés matériau. Le code généré
  (interfaces `generic` et `generic-parallel`) n'utilise que
  `mfront::gb::real`, et les lois sont des modèles paramétrés par leur type
  numérique : changer le `typedef` suffit pour le code généré.
- Aucun symbole n'indique le type d'une bibliothèque : une loi compilée en
  double chargée par un MGIS en float lirait des tableaux de travers.
- Le critère de convergence des lois implicites vaut 1e-8 par défaut ;
  MFront signale une valeur au-dessus de 1e-12, conseille 1e-14. En float
  (précision relative 1,2e-7), ces valeurs peuvent être inatteignables :
  à mesurer.

mfem-mgis :

- `mfem_mgis::real` vaut `mgis::real`, sans vérification de l'égalité avec
  `mfem::real_t` ; les tableaux des matériaux sont dans la mémoire de MFEM
  et passés à MGIS : les deux types doivent être les mêmes.
- 104 `double` dans `src` et `include`, la plupart légitimes (temps du
  profileur, paramètres donnés par l'utilisateur). À corriger : `MPI_DOUBLE`
  sur des données en `real` (`ComputeResultantForceOnBoundary.cxx`,
  `MeanThermodynamicForces.ixx`, distances de
  `ImposedDirichletBoundaryConditionAtClosestNode.cxx` et
  `PeriodicNonLinearEvolutionProblem.cxx` à vérifier) ;
  `Batch/GlobalStiffnessMatrix.cxx` utilise déjà `MPITypeMap<real>`.
- AmgX : les enveloppes de mfem-mgis sont en double ; MFEM refuse AmgX en
  simple, elles seront donc absentes.
- Tolérances en dur : une vingtaine de valeurs de 1e-8 ou moins dans les
  exemples et mm-opera-hpc, 41 fichiers de tests de mfem-mgis.

Essai fait sans rien modifier : le code généré de 14 lois des exemples et
de mm-opera-hpc, compilé avec `mfront_gb_real` en float (copie de
`Types.h` placée avant celle de TFEL) :

- 9 lois compilent, avec les mêmes avertissements qu'en double :
  Elasticity, IsotropicLinearElasticity, IsotropicLinearHardeningPlasticity,
  IsotropicLinearThermoElasticity, MonoCristal_UO2,
  SaintVenantKirchhoffElasticity, StationaryLinearHeatTransfer,
  U3SI2_IrradiationCreep, U3SI2_SolidSwelling ;
- 5 lois échouent, toutes dans leur propre code et pour la même raison :
  un littéral double multiplié par un tenseur en float donne une
  expression en double, que TFEL refuse d'affecter à un objet en float
  (`dfeel_ddeel += 2.*mu*...` dans
  NortonViscoplasticityWithThreshold d'ex7, `= 2 * ...` dans les
  CoupledHeatTransfer d'ex8, NortonCr et Mono_UO2_CosH_Jaco3 de
  mm-opera-hpc) ;
- le code généré pour le GPU (`generic-parallel`, nvcc) n'a pas été essayé.

## Proposition

### Source de vérité

Le type est celui de MFEM, lu dans `config.mk` par la configuration de
mfem-mgis (`MFEM_USE_SINGLE`). Il n'y a pas d'option de mfem-mgis pour le
choisir : changer de précision veut dire installer un autre MFEM (variante
spack `precision`), qui fixe aussi hypre, PETSc et MUMPS.

### MFront (TFEL)

- `Types.h` des lois et des propriétés matériau : le type vaut la macro
  `MFRONT_GB_REAL_TYPE` (`MFRONT_GMP_REAL_TYPE`) si elle est définie,
  `double` sinon. Rien ne change pour qui ne la définit pas.
- Une option de `mfront` (D2) écrit cette macro en tête des sources
  générées : le type est fixé dans le code produit, sans passer de
  drapeaux de compilation à chaque backend (CUDA, HIP, SYCL, Kokkos).
- Chaque fonction générée exporte la taille de son type réel (symbole
  `<fonction>_real_size`, calculé à la compilation, donc juste même si le
  drapeau a été oublié).
- Littéraux double dans le code des lois (D3).

### MGIS

- Option de CMake `MGIS_REAL_TYPE` (`double` par défaut, `float`), écrite
  dans un en-tête installé (`MGIS/RealType.h`, inclus par `Config-c.h`) :
  une application voit toujours le type de la bibliothèque qu'elle utilise.
  La macro actuelle, définie par l'application, devient une erreur si elle
  contredit l'en-tête.
- Au chargement d'une loi ou d'une propriété matériau : si le symbole
  `<fonction>_real_size` existe, il doit valoir `sizeof(mgis::real)` ;
  sinon la bibliothèque est ancienne, donc en double, refusée si
  `mgis::real` ne l'est pas. Message : nom de la loi, des deux types, et
  comment reconstruire.
- Paramètres et bornes : interface inchangée en double (D6).
- Liaisons en float : refusées à la configuration dans un premier temps
  (D5).

### mfem-mgis

- Configuration : lecture de `MFEM_USE_SINGLE` dans `config.mk`,
  vérification que MGIS a le même type (lecture de `MGIS/RealType.h`),
  erreur claire sinon ; `static_assert(std::is_same_v<mfem::real_t,
  mgis::real>)` dans `Config.hxx` en filet de sécurité.
- `mfem_mgis_add_behaviours_library` passe le type à `mfront` (option de
  D2) ; la variable `<cible>_REAL_TYPE` le donne aux exemples.
- `MPI_DOUBLE` remplacé par `mfem::MPITypeMap<real>::mpi_type` partout où
  les données sont en `real`.
- AmgX, UMFPACK et GSLIB déjà conditionnels à leur présence dans MFEM :
  vérifier que tout ce qui en dépend se compile sans eux (le chemin
  « float » les retire tous).
- Tolérances (D7) : une fonction `mfem_mgis::getDefaultTolerance(t)` qui
  rend `t` en double et `max(t, k ε)` en float, utilisée par les valeurs
  par défaut des solveurs et par les exemples ; les tests comparent aux
  références avec une tolérance déduite de ε.

### Exemples et mm-opera-hpc

- Tolérances par la fonction ci-dessus, références comparées avec une
  tolérance déduite du type.
- Lois corrigées selon D3.
- MTest ne charge que des lois en double (D4).

## Validation

1. En double : toutes les suites inchangées (CPU, GPU, exemples,
   mm-opera-hpc), mêmes résultats au bit près. C'est le chemin par défaut,
   il ne doit pas bouger.
2. En float, sur le CPU : environnement spack séparé (MFEM
   `precision=single`, hypre, MUMPS ; sans CUDA, compilation courte), MGIS
   en float ; suites de mfem-mgis et exemples, refus des mélanges (loi en
   double chargée en float, MGIS et MFEM différents).
3. En float, sur le GPU : environnement spack séparé avec CUDA (compilation
   longue de MFEM, voir `deps/spack/README.md`), sans AmgX.
4. Mesures sur la poutre : temps et nombre d'itérations float contre
   double, sur le CPU et sur le GPU ; écart de la solution ; convergence
   des lois implicites avec leur critère par défaut.

## Décisions à prendre

- D1. Source de vérité : le MFEM installé, pas d'option propre à
  mfem-mgis. Recommandé (demandé).
- D2. Comment MFront reçoit le type. (a) Macro seule, passée par les
  drapeaux de compilation de chaque backend. (b) Macro, plus une option de
  `mfront` qui l'écrit dans les sources générées. Recommandé : (b), un
  seul point d'entrée quel que soit le backend ; le symbole de taille
  protège dans les deux cas.
- D3. Littéraux double dans le code des lois. (a) TFEL accepte l'affectation
  d'une expression d'un type plus large dans un objet en float (`=`, `+=`,
  `-=`), avec conversion à l'affectation ; les lois existantes compilent
  telles quelles. (b) Réécrire les lois avec des littéraux typés
  (`real(2)`), TFEL reste strict. Recommandé : (a) dans TFEL, parce que les
  lois des utilisateurs sont écrites ainsi ; les 5 lois des exemples
  passent alors sans changement. C'est un changement de conception de
  TFEL, à proposer à ses auteurs avant toute publication.
- D4. MTest en float. (a) La fonction CMake construit aussi une
  bibliothèque en double pour les tests MTest (ex8). (b) Ces tests sont
  sautés en float. Recommandé : (a), seulement là où un test MTest existe.
- D5. Liaisons de MGIS (Python, Fortran, Julia) en float : refusées à la
  configuration pour commencer. Recommandé.
- D6. Interface des paramètres des lois (double) et des bornes (long
  double) : inchangée. Recommandé : ce sont des données d'entrée,
  converties à l'usage, et la changer casserait toutes les bibliothèques
  existantes.
- D7. Tolérances : valeur demandée en double, bornée par `k ε` en float,
  avec `k` à mesurer (100 au départ). Le critère des lois implicites
  (paramètre `epsilon`) n'est pas modifié tant que la mesure (validation,
  point 4) ne montre pas de non-convergence.

Avant de construire : accord sur les deux environnements spack séparés en
float (CPU, puis GPU). L'environnement principal et ses MFEM et AmgX ne
sont pas touchés.

## Étapes

1. TFEL : macro dans `Types.h`, option de `mfront`, symbole de taille,
   D3 ; tests de TFEL (une loi de chaque DSL compilée et intégrée en float
   par un petit programme, l'interface `generic-parallel` avec nvcc).
2. MGIS : option de CMake, en-tête installé, vérification au chargement,
   refus des liaisons ; tests de MGIS en float et en double, tests des
   refus.
3. mfem-mgis : configuration, `static_assert`, `MPITypeMap`, fonction CMake,
   tolérances ; double inchangé, float sur le CPU.
4. Exemples et mm-opera-hpc ; puis float sur le GPU.
5. Mesures, documentation (guide de mfem-mgis, pages de MGIS et de MFront).
