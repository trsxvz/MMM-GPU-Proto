# Mémoire de l'hôte de MFEM sur le GPU (11/10/2026)

Script : `cas/mesures-memoire-hote.sh <cas>`, un processus et le GPU,
machine calme, un passage.

## Le défaut

MFEM compilé avec Umpire (c'est le cas ici, comme hypre) prend Umpire
pour la mémoire de l'hôte dès qu'un périphérique est configuré
(`HOST_UMPIRE`). Umpire enregistre chaque allocation dans sa table. Les
boucles par élément faites sur l'hôte allouent et libèrent des
temporaires à chaque point (fonctions de forme, matrices élémentaires) :
chacun passe par cette table, sous verrou. Diagnostic par des piles
d'appels prises avec gdb pendant les post-traitements de la bulle :
toutes dans `umpire::...::AllocationMap`, appelé depuis `CalcShape` et
`CalcDShape`. Les appels CUDA ne comptaient que pour 1 s sur 61.

Correctif (mfem-mgis 0c61fbfcb4, branche memoire-hote) : avant de créer
le périphérique, `mfem_mgis::initialize` choisit la mémoire standard pour
l'hôte et garde Umpire pour le GPU, sauf si la variable `MFEM_MEMORY` de
MFEM est donnée.

## Variantes

- `umpire` : `MFEM_MEMORY=umpire`, Umpire pour l'hôte et le GPU (choix
  de MFEM, celui de toutes les mesures GPU précédentes du prototype) ;
- `host` : `MFEM_MEMORY=host`, mémoire standard pour l'hôte et le GPU ;
- `defaut` : choix de mfem-mgis après le correctif, mémoire standard
  pour l'hôte et Umpire pour le GPU.

Exécutables : mm-opera-hpc (`build/mm-opera-hpc-memoire`) et poutre
(`build/cas-memoire`) liés à `install/mfem-mgis-memoire` (correctif) ;
ex8 de `build/examples-pointwise-ex8`, lié à une installation sans le
correctif (sa variante `defaut` vaudrait `umpire`, elle n'est pas
mesurée). Lois cuda de `build/cas-memoire`, sauf mention point par point
(lois `generic`, intégrées sur l'hôte).

Résultats identiques entre variantes : bulle aux chiffres écrits,
cermet à 5,4e-15, polycristal à 2,5e-13 (70 valeurs).

## Temps (secondes, total du profil)

| Cas | umpire | host | defaut |
|---|---|---|---|
| bulle (467 109 inconnues, tout par noyaux, post-traitements compris) | 69,17 | 13,83 | 13,97 |
| cermet, 1 pas, tout par noyaux | 3,52 | 1,53 | 1,53 |
| cermet, 1 pas, point par point | 28,10 | 14,81 | 14,80 |
| ex8, 5 pas, tout par noyaux, lois cuda | 18,64 | 16,89 | |
| ex8, 5 pas, point par point | 164,57 | 57,30 | |
| polycristal, 10 pas, tout par noyaux | 13,26 | 11,76 | 12,32 |
| poutre r4, 10 pas, tout par noyaux | 3,55 | 3,22 | 3,31 |

- Le gain suit la part du travail fait sur l'hôte par des boucles par
  élément de MFEM. Bulle : post-traitements (export des contraintes aux
  nœuds, projection L2) de 60,8 s à 10,4 s (`host`) et 10,5 s
  (`defaut`), contre 11,3 s sur un cœur du CPU. ex8
  point par point : intégration, résidu et tangente (assemblés sur
  l'hôte) de 149 s à 41,7 s. ex8 tout par noyaux : la condition de Robin,
  évaluée sur l'hôte, fait passer le résidu de 1,47 à 0,31 s et la
  tangente de 1,10 à 0,42 s.
- Quand tout est calculé par noyaux sur le GPU (polycristal, poutre),
  l'écart reste de l'ordre du bruit d'un passage (de 3 à 11 %).
- `host` et `defaut` ne se distinguent pas ici : la mémoire du GPU
  d'Umpire n'apporte ni ne coûte rien de visible. Le correctif la garde,
  comme MFEM.
- Bulle sur le GPU après correctif : 14 s, contre 31 s sur un cœur du
  CPU, tout par noyaux. Avant, le GPU était deux fois plus lent que le
  CPU.
- Les mesures GPU antérieures du prototype ont toutes tourné avec
  `umpire`. Celles où le GPU calcule tout par noyaux changent peu.
  Celles qui gardent des boucles sur l'hôte étaient pénalisées : intégration
  point par point sur le GPU, post-traitements, autres intégrateurs
  évalués sur l'hôte.
