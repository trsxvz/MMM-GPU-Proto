# Passage à l'échelle de la poutre (11/10/2026)

Campagne de mesures, cas 5. Poutre du prototype raffinée trois, quatre et
cinq fois, 10 pas, meilleure méthode du lot L5 : NewtonK0 (gradient
conjugué préconditionné par le multigrille de K0, construit une fois),
prédiction par K0, tolérance d'Eisenstat-Walker 1, tous les noyaux.
Script : `cases/measure-beam-scaling.sh`. mfem-mgis e03fd927f7 (option de
la réserve de hypre dans la poutre, reprise depuis par fd93974cdc). Une
RTX 5080 (16 Go, 15,9 Gio utilisables), 8 cœurs (16 fils). Machine calme
pendant ces calculs.

Sur le GPU, la réserve de mémoire de hypre est réduite à 0,5 Gio
(`-hdps 0.5`, voir plus bas). Mémoire : maximum relevé par `nvidia-smi`
toutes les 200 ms, 686 Mio occupés par le bureau compris ; en Gio
(1 Gio = 1024 Mio).

| Raffinement | Inconnues | GPU, BoomerAMG | GPU, AmgX simple | GPU, sans matrice | 8 cœurs | Gain du GPU |
|---|---|---|---|---|---|---|
| 3 | 15 795 | 0,67 s | **0,48 s** | 0,70 s | 0,85 s | ×1,8 |
| 4 | 111 843 | 1,81 s | **1,54 s** | 2,01 s | 7,02 s | ×4,5 |
| 5 | 839 619 | **9,86 s** | 11,98 s | 12,12 s | 107,9 s | ×10,9 |

Temps de la boucle sur les pas ; le gain compare le meilleur temps du GPU
au calcul sur 8 cœurs. Itérations du point fixe : 33 à 38, presque
indépendantes de la taille. Mémoire du GPU au plus : 2,29 Gio à trois et
quatre raffinements (2,73 et 2,87 Gio avec AmgX), 6,43 Gio à cinq
(7,76 Gio avec AmgX, 6,72 Gio sans matrice).

- Le gain du GPU croît avec la taille : ×1,8 à 16 000 inconnues, ×10,9 à
  840 000. Le GPU passe de 4 à 5 raffinements en 5,5 fois plus de temps
  pour 7,5 fois plus d'inconnues (mieux occupé), le CPU en 15 fois plus
  de temps (la résolution passe de 3,4 s à 58 s : plus de données que les
  caches, hypothèse non mesurée).
- AmgX réglé en simple précision est le meilleur jusqu'à 110 000
  inconnues, BoomerAMG agressif au-delà (résolution 4,4 s contre 6,3 s à
  cinq raffinements).
- La tangente appliquée sans assemblage ne paie pas : l'assemblage évité
  (2,0 s à cinq raffinements) coûte plus en résolution (8,0 s contre
  4,4 s), et elle n'économise pas de mémoire (K0 reste assemblé).

Répartition (`cases/summary-profiling.py`, secondes, total avec la mise en
place) :

| Calcul | Total | Intégration | Tangente | Résolution | Prédiction |
|---|---|---|---|---|---|
| 5 raffinements, GPU, BoomerAMG | 12,0 | 1,1 | 2,0 | 4,4 | 2,6 |
| 5 raffinements, 8 cœurs | 108,1 | 3,1 | 22,0 | 58,4 | 23,6 |

## Six raffinements

6 502 275 inconnues : ne tient pas dans la mémoire du GPU. Avec la
réserve de hypre par défaut, l'allocation de la matrice de K0 (2,14 Gio)
est refusée à 13,46 Gio occupés (`poutre-r6-gpu-amg-agressif-assemblee`,
même refus sans matrice pour la tangente) ; avec une réserve de 0,5 Gio,
une allocation de 1,07 Gio est refusée à 15,34 Gio, la tangente étant
assemblée ou non (`poutre-r6-gpu-amg-agressif-reserve-hypre-0.5-*`). Les
sorties commencent par la ligne de commande et donnent le message
d'erreur, sauf `poutre-r6-gpu-amg-agressif-sans-matrice` (premier essai,
message seul). Mémoire de l'hôte : 18,3 Gio.

## Réserve de mémoire de hypre

hypre 3.2 alloue sur le GPU, par Umpire, une réserve de 4 Gio dès sa
première allocation, quelle que soit la taille du problème. Avec une
réserve de 0,5 Gio (qui grandit si besoin), la mémoire occupée à trois et
quatre raffinements passe de 5,79 Gio à 2,29 Gio (3,5 Gio de moins). Le
temps change de quelques % : avec BoomerAMG, 0,671 s contre 0,644 s à
trois raffinements, 1,807 s contre 1,776 s à quatre, 9,860 s contre
9,469 s à cinq (une itération de plus), soit 2 à 4 % plus lent ; avec
AmgX, moins de 1 %. Premier passage avec la réserve par défaut :
`../beam-scaling-2026-10-11-0233`. Ses calculs sur 8 cœurs, et leur
reprise, étaient plus lents de 13 à 43 %, pour une cause non établie (des
compilations d'une autre session ont été vues pendant ces calculs, mais
la reprise, à une charge de 1,1, a donné les mêmes temps) : seule la
série de ce dossier sert de référence pour le CPU.
L'option est devenue une option de mfem-mgis (`--hypre-device-pool-size`
ou variable `MFEM_MGIS_HYPRE_DEVICE_POOL_SIZE`).
