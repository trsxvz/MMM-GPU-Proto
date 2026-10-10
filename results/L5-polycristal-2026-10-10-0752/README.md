# Polycristal, contrôle après la refonte de la géométrie (10/10/2026)

Mêmes cas que `results/L5-polycristal-2026-10-10-0502/`, avec mfem-mgis
3e3cc5b0f2 : jacobiens calculés par un noyau à partir des nœuds ou des
sommets des éléments, matrices de rigidité calculées par élément sur
l'hôte (62908a890d). Charge 2,6 au départ.

| Cas | CPU, 8 processus | GPU |
|---|---|---|
| Newton | 29,35 s | 11,85 s |
| Newton, EW1 | 33,46 s | 11,57 s |
| NewtonK0 | 24,69 s | 11,23 s |
| NewtonK0, EW1 | **22,83 s** | **9,56 s** |

Temps totaux (`cas/resume-profilage.py`).

- Contraintes macroscopiques égales à Newton à 1,1e-14 près sans
  Eisenstat-Walker, de 9,6e-7 à 3,5e-6 avec : comme avant.
- Sur CPU, l'assemblage de la tangente passe de 2,01 s à 1,60 s (Newton,
  98 tangentes) grâce au calcul des matrices par élément sur l'hôte ;
  NewtonK0 avec Eisenstat-Walker 1 de 23,50 s à 22,83 s.
- Sur le GPU, rien ne change (le noyau des matrices y garde un fil par
  paire de nœuds).
