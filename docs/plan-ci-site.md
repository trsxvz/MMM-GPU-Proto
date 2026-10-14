# Plan : intégration continue et site du prototype

Statut : plan du 13/10/2026 (tâches 9f et 9g), à valider avant de coder.
Les deux tâches dépendent d'un même choix, la visibilité du dépôt
(décision E1).

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
   MFEM avec MUMPS et SuiteSparse), TFEL des sources du dépôt (ccache),
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

`.github/workflows/cpu.yml` et `deps/spack/ci-cpu/spack.yaml`, pas encore
lancés. Déclenchement manuel seulement tant que le dépôt est privé
(minutes des machines) ; poussées, demandes de fusion et lancement de nuit
à ajouter une fois le dépôt public.

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
- Sans cache des dépendances, chaque lancement reconstruit OpenMPI,
  hypre, MUMPS, SuiteSparse et MFEM, puis TFEL : plusieurs heures sur une
  machine à deux cœurs. Le cache sur GHCR (comme mfem-mgis) vient ensuite.
- Le float (9b) n'est pas encore couvert : un second job avec le jumeau de
  `deps/spack/mmm-float-cpu` viendra après le premier lancement réussi.

## Étapes, une fois E1 à E4 décidées

1. Environnement `ci-cpu`, job CPU sans cache, puis avec le cache GHCR.
2. Jobs CUDA et HIP de nuit.
3. Construction du site (générateur à choisir : Sphinx, déjà utilisé par
   mfem-mgis, peut porter l'accueil et l'onglet du prototype), publication
   par la CI.
4. Documentation de la CI et du site dans `docs/`.
