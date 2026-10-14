# Polycristal, contrôle après les jacobiens par noyau (10/10/2026)

Mêmes cas que `results/L5-polycrystal-2026-10-10-0502/`, avec le calcul
des jacobiens par un noyau (mfem-mgis 8b0b47acfb ; `summary.txt` indique le
commit de base). Charge 4,8 au départ (fin des mesures de la poutre).

| Cas | CPU, 8 processus | GPU |
|---|---|---|
| Newton | 29,65 s | 11,86 s |
| Newton, EW1 | 34,08 s | 11,60 s |
| NewtonK0 | 25,12 s | 11,25 s |
| NewtonK0, EW1 | 23,50 s | 9,58 s |

Temps totaux (`cases/summary-profiling.py`). Tout est égal à la mesure
précédente à 0,5 % près : le maillage est petit, la géométrie n'y pesait
pas. Contraintes macroscopiques égales à Newton à 2,4e-13 près sans
Eisenstat-Walker, de 9,6e-7 à 3,5e-6 avec (dans la tolérance du Newton
microscopique, voir la mesure précédente).
