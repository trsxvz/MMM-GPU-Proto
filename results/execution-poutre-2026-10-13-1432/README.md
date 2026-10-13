# Placements des calculs par l'option --execution, poutre (13/10/2026)

Poutre r4 du prototype (111 843 inconnues), 10 pas, Newton, PCG de hypre
préconditionné par BoomerAMG (tolérance 1e-4), loi de von Mises à
écrouissage linéaire, construite pour le CPU (interface `generic`) et pour
le GPU (interface `generic-parallel`, backend `cuda`). mfem-mgis 8728b7979f
(installation séparée `install/mfem-mgis-9a`, construction
`build/cas-9a`). Script : `cas/mesures-execution-poutre.sh`. Machine calme
(la session voisine en pause, GPU libre).

Tout sur le CPU : `--execution cpu`, 8 processus. Les autres calculs : un
processus et le GPU, `--execution gpu`, ou `--execution mixed` avec le
placement donné des éléments finis (EF), des lois et du solveur linéaire
(la poutre n'a pas de post-traitements).

| Placement | EF | Lois | Solveur | Temps (s) | Copies hôte vers GPU | Copies GPU vers hôte |
|---|---|---|---|---|---|---|
| `cpu`, 8 processus | CPU | CPU | CPU | 12,17 | | |
| `gpu` | GPU | GPU | GPU | 2,96 | 120, 598 Mo, 0,035 s | 120, 598 Mo, 0,036 s |
| `mixed` | GPU | GPU | GPU | 2,87 | 120, 598 Mo, 0,034 s | 120, 598 Mo, 0,036 s |
| `mixed` | GPU | CPU | GPU | 3,73 | 184, 2 803 Mo, 0,139 s | 82, 387 Mo, 0,027 s |
| `mixed` | GPU | GPU | CPU | 42,74 | 120, 598 Mo, 0,033 s | 120, 598 Mo, 0,036 s |
| `mixed` | GPU | CPU | CPU | 44,66 | 176, 2 671 Mo, 0,146 s | 78, 368 Mo, 0,024 s |
| `mixed` | CPU | GPU | GPU | 30,62 | 624, 3 108 Mo, 0,173 s | 372, 4 307 Mo, 0,242 s |
| `mixed` | CPU | CPU | GPU | 31,17 | | |
| `mixed` | CPU | GPU | CPU | 66,00 | 576, 2 869 Mo, 0,156 s | 348, 3 998 Mo, 0,226 s |
| `mixed` | CPU | CPU | CPU | 66,90 | | |

- Les dix placements donnent la même solution (norme des déplacements
  48,83704586). Le nombre d'itérations de Newton varie de 28 à 33 : le
  BoomerAMG de hypre n'est pas le même sur le CPU et sur le GPU, et la
  tolérance de 1e-4 du PCG rend Newton inexact.
- Les copies des tableaux des matériaux coûtent au plus 0,42 s (éléments
  finis sur le CPU et lois sur le GPU) : le surcoût des placements mixtes
  vient du processeur qui calcule, pas des transferts. Avec un seul
  processus, une phase placée sur le CPU n'a qu'un cœur : le solveur
  linéaire sur le CPU ajoute environ 40 s, les éléments finis élément par
  élément environ 28 s, l'intégration des lois sur le CPU 0,86 s (dont
  0,17 s de copies).
- Tout sur le GPU, 120 copies dans chaque sens (0,07 s, 2,5 % du calcul)
  restent, soit 12 par pas : hypothèse à vérifier, la mise à jour des états
  à la fin de chaque pas, faite sur l'hôte, qui recopie les tableaux sur
  l'hôte, le pas suivant les renvoyant sur le GPU. Piste : mettre à jour
  les états sur le GPU quand les lois y sont intégrées.
- Les compteurs ne voient que les tableaux des matériaux : les propriétés
  matériau et les variables externes copiées vers le GPU à chaque
  intégration, et les recopies faites par MFEM (vecteurs, matrices données
  à hypre) ne sont pas comptées. Le surcoût total d'un placement se lit
  donc dans les temps, les compteurs en donnant la part des tableaux.
- Tout sur le CPU, les éléments finis sont maintenant calculés par noyaux
  sur l'hôte par défaut : 12,17 s sur 8 processus, contre 16,47 s élément
  par élément (première série, `results/execution-poutre-2026-10-13-1418`,
  faite avant ce changement ; 12,43 s par noyaux avec les variables
  d'environnement). Sur 2 processus, élément par élément : 37,38 s.
