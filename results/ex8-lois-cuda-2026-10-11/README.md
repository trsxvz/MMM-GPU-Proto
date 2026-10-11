# ex8 entièrement sur le GPU, lois compilées pour le GPU (11/10/2026)

Script : `cas/mesures-ex8.sh`. Exemple ex8 de mfem-mgis-examples (plaque
thermomécanique : thermique avec condition de Robin, mécanique avec
gonflement, couplage itératif), maillage d'origine (8424 hexaèdres), 5 pas,
MUMPS pour la mécanique, GMRES et BoomerAMG pour la thermique.
Construction séparée `build/examples-pointwise-ex8` (ex8 seul, mfem-mgis-examples
6505b78) liée à mfem-mgis 8011fe9eaa (`install/mfem-mgis-pointwise` :
modèle aux points intégré sur le périphérique). Lois de l'exemple compilées
par l'interface `generic` (intégrées sur l'hôte), ou par l'interface
`generic-parallel` pour le GPU (backend cuda, cas/ex8, intégrées par lots
sur le GPU, modèle de gonflement compris). CPU : 2 processus d'un fil ;
GPU : un processus. Machine calme (charge 0,95), GPU libre.

Couplage sur les 5 pas, en secondes :

| Cas | Couplage |
|---|---|
| CPU, 2 processus, point par point | 30,51 |
| CPU, 2 processus, tout par noyaux | 22,13 |
| GPU, point par point | 165,99 |
| GPU, tout par noyaux, lois `generic` | 37,78 |
| GPU, tout par noyaux, lois `generic-parallel` (cuda) | 18,45 (répété : 18,53) |

Phases (`cas/resume-profilage.py`, deux problèmes cumulés) :

| Cas | Intégration | Résidu | Tangente | Factorisation (MUMPS) et préconditionneur | Résolution |
|---|---|---|---|---|---|
| GPU, lois `generic` | 19,57 | 1,59 | 1,40 | 11,63 | 3,56 |
| GPU, lois cuda | 0,66 | 1,51 | 1,11 | 11,61 | 3,52 |

- Mêmes statistiques des champs que point par point sur CPU (température,
  norme du déplacement, gonflement), mêmes nombres d'itérations de Newton
  (64 itérations, 114 évaluations du résidu).
- Avec les lois compilées pour le GPU, l'intégration passe de 19,6 à
  0,66 s (×30) ; ex8 tourne en 18,5 s sur le GPU, contre 166 s point par
  point sur le GPU et 22,1 s tout par noyaux sur 2 processus CPU.
- Il reste MUMPS (factorisation et résolutions sur l'hôte) : 15,1 s sur
  18,6 s.
- Sans la synchronisation des gradients mécaniques ajoutée à
  FieldUpdaterModel (mfem-mgis-examples 6505b78), la thermique lit sur
  l'hôte des gradients périmés et les statistiques s'écartent de la
  référence d'environ 3e-5 en relatif (essai sur le maillage grossier).
