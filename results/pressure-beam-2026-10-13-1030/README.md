# Pression imposée sur la poutre, lot L7.6 (13/10/2026)

Poutre r4 du prototype (111 843 inconnues), pression de 1e8 Pa sur
l'extrémité, 10 pas, NewtonK0 (tolérance d'Eisenstat et Walker, un cycle
de BoomerAMG sur K0), tous les calculs par noyaux, sur le GPU (un
processus, hypre sur le GPU) et sur 8 processus CPU. Script :
`cases/measure-pressure-beam.sh`.

- Avant : exécutable du 11/10 (`scratch/nuit-2026-10-11/l76/poutre-avant`,
  version de mfem-mgis antérieure au lot L7.6) : la pression est un
  intégrateur de bord de la forme non linéaire, évalué sur l'hôte à chaque
  résidu ; sa présence oblige aussi à assembler la tangente par MFEM sur
  l'hôte, à partir des matrices des éléments calculées par noyaux.
- Après : mfem-mgis 54ff8e6c8a (installation séparée
  `install/mfem-mgis-54ff8e6c8a`, construction `build/cas-54ff8e6c8a`) :
  la pression, indépendante des inconnues, est assemblée une fois pour une
  pression unité et ajoutée au résidu sur le GPU, sans intégrateur dans la
  forme ; la tangente reste assemblée par noyaux, ou appliquée sans
  assemblage.

Temps de la boucle sur les pas (secondes) :

| | Sans pression, avant | Avec pression, avant | Sans pression, après | Avec pression, après |
|---|---|---|---|---|
| GPU | 1,81 | 17,63 | 1,50 | 1,70 |
| 8 processus | 7,13 | 9,88 | 5,70 | 5,68 |

Avec la tangente appliquée sans assemblage et la pression : 1,86 s sur le
GPU.

- Sur le GPU, la pression coûtait 15,8 s (tangente 14,3 s, assemblée par
  MFEM sur l'hôte ; résidu 1,8 s) ; elle ne coûte plus que 0,2 s : le
  calcul avec pression est 10 fois plus rapide.
- Sur 8 processus, le surcoût de la pression (2,75 s) disparaît.
- Sans pression, la nouvelle version est aussi plus rapide (1,81 → 1,50 s
  sur le GPU, 7,13 → 5,70 s sur 8 processus) : lots du 11/10 postérieurs à
  l'exécutable d'avant (positions des entrées, matrices sans former B,
  inverses des jacobiennes).
- Un premier passage, fait pendant des compilations d'une autre session,
  donnait des temps CPU faux (23 s au lieu de 7 s) ; il a été refait sur
  une machine calme.
