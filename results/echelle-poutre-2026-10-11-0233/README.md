# Passage à l'échelle de la poutre, premier passage (11/10/2026)

Premier passage de `cas/mesures-echelle-poutre.sh`, avec la réserve de
mémoire de hypre par défaut (4 Gio) et sans la variante sans matrice.
Une compilation tournait en même temps : les calculs sur 8 cœurs (et leur
reprise, fichiers `*-reprise.txt`) sont plus lents de 15 à 40 % que sur
une machine calme. Résultats de référence et analyse :
`../echelle-poutre-2026-10-11-0247`.
