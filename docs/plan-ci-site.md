# Plan : intégration continue et site du prototype

Statut : plan du 13/10/2026 (tâches 9f et 9g), décisions prises le
13/10 (voir plus bas) ; ébauches de la CI sur le CPU et du site le
14/10, déclenchées à la main, pas encore lancées. Les deux tâches
dépendent de la visibilité du dépôt (décision E1).

## Contraintes

- Le dépôt est privé. Avec un compte GitHub gratuit (conditions à
  vérifier au moment du choix) : environ 2 000 minutes d'Actions par mois
  pour les dépôts privés, environ 500 Mo de stockage pour les paquets
  privés (le cache de spack de la CI de mfem-mgis pèse près de 2 Go), et
  GitHub Pages n'est pas disponible pour un dépôt privé. Un dépôt public
  n'a pas ces limites, mais publie tout son contenu.
- Les machines de GitHub n'ont pas de GPU : la CI peut compiler pour CUDA
  et pour HIP, pas exécuter sur un GPU.
- Coûts mesurés sur la CI de mfem-mgis (machines à 4 cœurs) : dépendances
  sans cache 70 min, avec cache 2 min. Ici s'ajoutent TFEL (le plus long,
  compilé depuis les sources du dépôt), MGIS, les exemples et mm-opera-hpc,
  et pour CUDA un MFEM compilé avec CUDA (plusieurs heures la première
  fois, voir `deps/spack/README.md`).
- L'environnement spack mmm-gpu est trop gros pour la CI (NVIDIA HPC SDK,
  LLVM, Kokkos, PETSc, Ginkgo, SUNDIALS...). La CI a ses propres
  environnements, réduits.

## Proposition pour la CI (9f)

Trois niveaux, chacun dans son job :

1. CPU, à chaque poussée : environnement spack `ci-cpu` (OpenMPI, hypre,
   MFEM avec MUMPS et SuiteSparse), et `ci-cpu-float` pour la simple
   précision, TFEL des sources du dépôt (ccache),
   MGIS appliqué depuis `patches/mgis`, mfem-mgis, exemples, mm-opera-hpc ;
   tests sur CPU de chacun. Vérifie aussi que les séries de patchs
   s'appliquent (`scripts/verify-patches.sh`).
2. CUDA, compilation seule, chaque nuit : environnement `ci-cuda` (CUDA de
   spack, MFEM avec CUDA, sans GPU), noyaux de mfem-mgis et backends
   `cuda` de TFEL compilés, tests sur CPU.
3. HIP, compilation seule, chaque nuit : ROCm par paquets Python et
   `scripts/build-hip-amd.sh` (MFEM HIP, mfem-mgis, exemples), tests sans
   GPU.

Cache : dépendances de spack sur GHCR comme pour mfem-mgis (actions
composites, index unique, vidage mensuel), plus ccache pour TFEL. Les
leçons de la CI de mfem-mgis sont reprises : versions épinglées de
spack-packages (sinon une mise à jour d'autoconf invalide tout le cache),
pas de `continue-on-error` sur l'envoi au cache, jetons en lecture seule
pour les demandes de fusion venues de forks.

Les tests sur GPU restent sur la machine de développement
(`scripts/mfem-mgis-gpu-suite.sh`), sauf décision E3.

## Proposition pour le site (9g)

Un seul site GitHub Pages pour le dépôt (un site par dépôt), construit et
poussé par la CI :

- accueil : ce qu'est le prototype, liens vers les trois parties ;
- TFEL, MGIS, mfem-mgis : leurs documentations construites depuis les
  sources du prototype (pages web de TFEL et de MGIS, Sphinx et Doxygen de
  mfem-mgis) ;
- un onglet « What's new with this prototype » : différences avec les
  versions officielles, ajouts, choix et justifications, corrections
  trouvées, pièges rencontrés, et le journal en dernier sous-onglet. Ces
  pages sont tirées de `README.md`, `docs/` et `JOURNAL.md`, en anglais
  comme les documentations des trois bibliothèques (traduction à faire :
  les documents du prototype sont en français).

## Décisions à prendre

- E1. Visibilité. (a) Dépôt public : minutes et paquets gratuits, Pages
  possible, tout le contenu publié. (b) Dépôt privé avec un compte payant
  (Pages exige alors un site public quand même). (c) Dépôt privé, CI sur
  une machine à soi (runner auto-hébergé), pas de site. Recommandé : (a)
  si la publication est voulue ; sinon (c) pour la CI, sans site.
- E2. Niveaux de la CI : les trois, ou seulement le niveau CPU au début.
  Recommandé : CPU d'abord, puis CUDA et HIP la nuit.
- E3. Tests sur GPU par un runner auto-hébergé sur la machine de
  développement (partage du GPU avec les autres travaux, sécurité : jamais
  pour des demandes de fusion venues de forks). Recommandé : non pour
  l'instant.
- E4. Langue du site : anglais (comme les documentations des
  bibliothèques), le dépôt gardant ses documents en français. Recommandé.

## Décisions prises (13/10/2026)

E1 (a) : dépôt public, après un audit de tout l'historique (fait le
13/10, rien à réécrire) ; le changement de visibilité reste à faire par le
propriétaire du dépôt. E2 : CPU d'abord, puis CUDA et HIP la nuit. E3 :
pas de runner sur la machine de développement. E4 : site en anglais.

## Ébauche du niveau CPU (14/10/2026)

`.github/workflows/cpu.yml`, `deps/spack/ci-cpu/spack.yaml` et
`deps/spack/ci-cpu-float/spack.yaml`, pas encore lancés. Déclenchement
manuel seulement tant que le dépôt est privé (minutes des machines) ;
poussées, demandes de fusion et lancement de nuit à ajouter une fois le
dépôt public.

- Job `patches` : chaque série de `patches/` est appliquée sur le dépôt
  amont à sa base, et l'arbre obtenu est comparé aux sources de la racine
  (`diff -r`). `scripts/verify-patches.sh` fait la même vérification en
  local, sur les dépôts imbriqués, qui n'existent pas sur la CI. Simulée en
  local pour les cinq séries : toutes s'appliquent et redonnent l'arbre de
  la racine.
- Job `cpu` : dépendances par spack (paquets épinglés au même commit que
  les environnements du prototype), TFEL des sources de la racine (backends
  hôte `stlpar` et `openmp` seulement, sans TDLS), MGIS appliqué depuis
  `patches/mgis`, puis mfem-mgis, les exemples et mm-opera-hpc, avec leurs
  tests sur le CPU ; ccache par les actions composites de mfem-mgis.
  OpenMPI 5.0.10 : `OMPI_MCA_btl=self,tcp` contre la collision des
  identifiants de travaux des `mpirun` concurrents, et surallocation
  permise (un test lance trois processus).
- Job `cpu` simulé en local le 14/10, de bout en bout, dans un
  environnement vide : même spack (1.2.2) et même commit de
  spack-packages, copie de `ci-cpu` avec gcc 14, puis TFEL, MGIS (49/49),
  mfem-mgis (208/208), exemples (22/22) et mm-opera-hpc (6/6). La
  simulation a trouvé un défaut, corrigé : la bibliothèque de lois des
  tests de mfem-mgis (cible `BehaviourTest`) n'était pas construite, d'où
  126 échecs sur 208.
- Sans cache des dépendances, chaque job de la matrice reconstruit ses
  propres dépendances, OpenMPI, hypre, MUMPS et MFEM (et SuiteSparse en
  double seulement : MFEM le refuse en float), puis TFEL : plusieurs
  heures sur une machine à deux cœurs. Le cache sur GHCR (comme
  mfem-mgis) vient ensuite.
- Float (9b) : le job `cpu` est une matrice sur la précision (`double`,
  `float`). L'environnement `deps/spack/ci-cpu-float` est le jumeau de
  `mmm-float-cpu`, et MGIS reçoit `MGIS_REAL_TYPE` de la matrice. Lancé à
  la main, le job ne prend que le double par défaut (entrée
  `precisions`) ; les deux précisions tourneront à chaque poussée une fois
  le dépôt public. Simulé en local en float : MGIS 49/49, mfem-mgis
  204/204, exemples 21/21, mm-opera-hpc 6/6. La simulation a trouvé un
  défaut de mm-opera-hpc en float : le point fixe sur les contraintes
  latérales du premier pas du cermet ne descendait pas sous son seuil de
  1e4 Pa (stagnation vers 1,1e4 Pa ou oscillation entre 1,5e4 et
  1,7e4 Pa, selon l'exécutable et le nombre de processus). Les
  composantes latérales de F, juste sous 1, ne sont connues qu'à un ulp
  près (6e-8), et la correction du point fixe s'arrondit à zéro sous
  environ 1,5e4 Pa de résidu. Le pas était coupé, et la contrainte axiale
  du cermet à 0,2 s s'écartait de 2,7 % du double, effet du pas de temps
  et non de l'arrondi (0,27 % du double aux mêmes pas). Le seuil vaut
  maintenant 3e4 Pa en float : écart au double de 0,10 à 0,16 % à 0,2 s
  (surtout l'arrondi de la déformation imposée, 1 + 1e-4 en float), et le
  test vérifie les mêmes instants dans les deux précisions
  (`results/cermet-float-2026-10-14`).
- Les lois générées par mfront dans les tests sont compilées avec
  `tfel-config --oflags`, qui contient `-march=native` même en
  construction portable : leurs résultats dépendent du processeur de la
  machine de la CI. Sur le cermet en float, les lois compilées par gcc 13
  et par gcc 14 donnent les mêmes contraintes à 1e-5 près ; le point
  d'arrêt du point fixe dépend de l'arrondi, ce qui déplace la contrainte
  axiale d'au plus 0,08 % avec le seuil de 3e4 Pa, loin de la tolérance
  de 1 % du test.

## Ébauche du site (14/10/2026)

- Générateur : Sphinx avec MyST (Markdown), thème `sphinx_rtd_theme`,
  comme la documentation de mfem-mgis. Les pages du prototype, en
  anglais, sont dans `site/` : accueil et onglet « What's new with this
  prototype » (vue d'ensemble, différences avec les versions officielles,
  ajouts, choix, défauts corrigés, pièges, jalons). Elles résument les
  documents en français du dépôt, qui restent la référence.
- `scripts/build-site.sh` construit le site complet dans
  `build/site/html` : ces pages à la racine, les pages web de TFEL et de
  MGIS (cible CMake `website`, pandoc 3.9 et pandoc-crossref 0.3.24, qui
  doivent être construits avec la même version de pandoc), et la
  documentation Sphinx et Doxygen de mfem-mgis. Rien n'est compilé : TFEL
  et MGIS sont seulement configurés. Les sources de MGIS patchées doivent
  être dans `mgis/` (la CI applique `patches/mgis`). Paquets Python :
  `site/requirements.txt`. Construit en local sans avertissement des
  pages du prototype (145 Mo).
- `.github/workflows/pages.yml` installe les outils, applique les patchs
  de MGIS, lance le script et publie sur GitHub Pages. Déclenchement
  manuel seulement : Pages demande un dépôt public.

## Étapes, une fois E1 à E4 décidées

1. Environnements `ci-cpu` et `ci-cpu-float`, job CPU sans cache, puis
   avec le cache GHCR.
2. Jobs CUDA et HIP de nuit.
3. Construction du site (Sphinx, voir « Ébauche du site »), publication
   par la CI.
4. Documentation de la CI et du site dans `docs/`.
