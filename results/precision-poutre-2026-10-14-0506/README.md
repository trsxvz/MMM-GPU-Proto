# Précision float contre double, premier passage (14/10/2026)

Même montage que `results/precision-poutre-2026-10-14-0514`, avec le
critère de Newton du cas (1e-8) et 1,2e-5 pour les deux précisions. Les six
calculs en float échouent : la tolérance relative, bornée à 100 ε = 1,2e-5
par mfem-mgis en float, est sous le plancher du résidu relatif à r4
(environ 1,3e-5 à partir du sixième pas, voir
`float-critere-du-cas-1.txt`), et Newton s'arrête au bout de 50 itérations.
Les calculs en double servent de référence au second passage.
