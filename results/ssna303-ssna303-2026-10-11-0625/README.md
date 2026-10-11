# Barreau entaillé ssna303 (11/10/2026)

Campagne de mesures, cas 2. Exemple ex4 de mfem-mgis-examples, variante
résolue par hypre (Newton, FGMRES de hypre préconditionné par BoomerAMG),
plasticité en grandes transformations (déformation logarithmique de
Hencky), chargement complet en 10 pas, maillage raffiné zéro, une et deux
fois. Script : `cas/mesures-ssna303.sh`. mfem-mgis f0324f1d9c installé
(FBar par noyaux, avant les corrections des relectures F et G). Une RTX
5080 (16 Go), 8 cœurs (16 fils), 30 Go de mémoire. Machine calme pendant
ces calculs.

Série INCOMPLÈTE : arrêtée à 7 h 18 pendant le calcul
`ssna303-r2-fbar-np8-noyaux`, la mémoire de la machine étant presque
épuisée (6 Gio d'échange utilisés après coup). Manquent ce calcul et
`ssna303-r2-fbar-gpu-noyaux`.

| Raffinement | Inconnues | 8 cœurs, point par point | 8 cœurs, noyaux | GPU, noyaux | Gain du GPU |
|---|---|---|---|---|---|
| 0 | 27 945 | 5,97 s | 6,08 s | 6,10 s | ×1,0 |
| 1 | 196 047 | 66,8 s | 65,5 s | 27,1 s | ×2,4 |
| 2 | 1 461 915 | 683 s | 670 s | échec | |

FBar (`--use-fbar`, par noyaux depuis le lot L7.9) :

| Raffinement | 8 cœurs, point par point | 8 cœurs, noyaux | GPU, noyaux | Gain du GPU |
|---|---|---|---|---|
| 0 | 7,40 s | 6,97 s | 6,34 s | ×1,1 |
| 1 | 85,5 s | 78,7 s | 29,5 s | ×2,7 |
| 2 | 1420 s | arrêté | non lancé | |

Temps de la boucle sur les pas ; le gain compare le GPU aux noyaux sur 8
cœurs.

- Sur le GPU, la résolution (FGMRES et BoomerAMG de hypre) domine : 4,5 s
  sur 6,1 s à zéro raffinement, 18 s sur 27 s à un raffinement. Le gain
  vient surtout de l'intégration (1,2 s contre 11,2 s), de la tangente
  (3,0 s contre 10,6 s) et de la prédiction.
- Sur 8 cœurs, les noyaux gagnent sur l'intégration et le résidu, mais la
  tangente par noyaux est plus lente que la tangente élément par élément
  (10,6 s contre 8,4 s à un raffinement), comme déjà noté pour la
  thermique (recherche des colonnes dans la matrice creuse).
- FBar coûte 20 % de plus que la formulation standard sur 8 cœurs, 9 % sur
  le GPU, à un raffinement.
- Deux raffinements, FBar, point par point : 1420 s, avec un pas
  redécoupé (48 résolutions contre 25 sans FBar) et peut-être de
  l'échange mémoire : temps à ne pas comparer.
- Deux raffinements sur le GPU, sans FBar : échec dès le premier pas (10
  redécoupages, puis arrêt), sans message dans la sortie. ex4 affiche
  depuis la cause de l'échec de chaque résolution (mfem-mgis-examples
  db0ce2a) : à relancer pour la connaître.

Répartition à un raffinement (`cas/resume-profilage.py`, secondes, total
avec la mise en place) :

| Calcul | Total | Intégration | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|
| 8 cœurs, point par point | 66,8 | 12,6 | 8,4 | 29,4 | 14,5 |
| 8 cœurs, noyaux | 65,5 | 11,2 | 10,6 | 28,6 | 14,3 |
| GPU, noyaux | 27,2 | 1,2 | 3,0 | 18,1 | 3,6 |
| FBar, GPU, noyaux | 29,6 | 1,2 | 3,2 | 20,0 | 3,8 |
