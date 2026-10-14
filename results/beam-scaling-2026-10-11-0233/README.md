# Passage à l'échelle de la poutre, premier passage (11/10/2026)

Premier passage de `cases/measure-beam-scaling.sh`, avec la réserve de
mémoire de hypre par défaut (4 Gio) et sans la variante sans matrice.
Les calculs sur 8 cœurs (et leur reprise, fichiers `*-reprise.txt`) sont
plus lents de 13 à 43 % que ceux du passage suivant, pour une cause non
établie (des compilations d'une autre session ont été vues pendant ces
calculs, mais la reprise, à une charge de 1,1, a donné les mêmes temps).
Résultats de référence et analyse :
`../beam-scaling-2026-10-11-0247`.
