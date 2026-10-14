# Autres intégrateurs de la forme dans les calculs par noyaux (11/10/2026)

Script : `cases/measure-other-integrators.sh`. Thermique r5 (279 873
inconnues, 5 pas), 8 processus d'un fil, tout par noyaux, tangente
assemblée, GMRES et BoomerAMG. Seconde extrémité : échange convectif
(option `-hc 1`, intégrateur de faces de bord de la forme non linéaire,
comme la condition de Robin de l'exemple ex8) ou température imposée (sans
autre intégrateur, pour comparaison). Avant : mfem-mgis e111df8454
(`install/mfem-mgis-globale`, boucles de MFEM sur tous les éléments dès
qu'un autre intégrateur est présent) ; après : 517e33cde4
(`install/mfem-mgis-autres`). Machine calme, trois passages (`summary.txt`,
`repetition-2`, `repetition-3`).

Boucle sur les pas (ligne `time`), en secondes, médiane et extrêmes des
trois passages :

| Cas | Avant | Après |
|---|---|---|
| échange convectif | 4,17 (4,10 à 4,17) | 2,70 (2,66 à 2,96) |
| température imposée | 2,98 (2,86 à 3,04) | 2,86 (2,85 à 2,91) |

Phases du premier passage (`cases/summary-profiling.py`, nombre d'appels
entre parenthèses) :

| Cas | Résidu | Tangente | Résolution |
|---|---|---|---|
| échange convectif, avant | 0,49 (20) | 1,62 (15) | 1,53 |
| échange convectif, après | 0,25 (20) | 0,49 (15) | 1,68 |
| température imposée, après | 0,23 (23) | 0,47 (18) | 1,57 |

- Avec l'échange convectif, la tangente passe de 108 à 33 ms par appel
  (×3,3) et le résidu de 24,5 à 12,5 ms (×2) ; boucle sur les pas ×1,54
  (médianes). Par appel, il reste environ 25 % de plus que sans autre
  intégrateur (tangente 33 contre 26 ms, résidu 12,5 contre 10 ms) :
  évaluation des faces de bord sur l'hôte et prolongation des inconnues.
- Mêmes normes de la solution et mêmes itérations de Newton avant et
  après (15 pour l'échange convectif, 18 pour la température imposée).
- Sans autre intégrateur, le chemin ne change pas : les écarts de la
  température imposée (moins de 4 %) sont du bruit de mesure.
- La résolution linéaire domine maintenant le calcul (1,5 à 1,7 s sur
  2,7 à 3 s).
