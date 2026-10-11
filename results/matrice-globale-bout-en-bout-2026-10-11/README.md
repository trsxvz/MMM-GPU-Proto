# Matrice globale avec motif gardé : bout en bout sur 8 processus (11/10/2026)

Script : `cas/mesures-matrice-globale.sh`. Thermique r5 (5 pas) et poutre
r4 (10 pas), 8 processus d'un fil, tangente assemblée par noyaux. Avant :
mfem-mgis 7af829b3a5 (`install/mfem-mgis-jacobiennes`) ; après :
e111df8454 (`install/mfem-mgis-globale`, `batch::GlobalStiffnessMatrix`).
La ligne « mfem-mgis » de `resume.txt` donne la version des sources, pas
celle des installations. Machine calme.

Temps en secondes (`cas/resume-profilage.py`) :

| Cas | Boucle sur les pas | Tangente | Résolution | Prédiction |
|---|---|---|---|---|
| thermique, avant | 3,38 | 1,04 (18) | 1,55 | |
| thermique, après | 2,74 | 0,43 (18) | 1,53 | |
| poutre, avant | 12,30 | 0,80 (29) | 8,77 | 2,44 |
| poutre, après | 11,84 | 0,44 (29) | 8,74 | 2,38 |

- Tangente ×2,4 (thermique) et ×1,8 (poutre) ; boucle sur les pas ×1,23
  et ×1,04 (la résolution domine la poutre).
- Mêmes itérations de Newton et mêmes normes de la solution.
- Depuis le début de la journée (tangente point par point, 3,12 s pour la
  thermique, 2,25 s pour la poutre, `results/positions-2026-10-11-1008`),
  la phase tangente sur 8 processus est 7 et 5 fois plus rapide.
