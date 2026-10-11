# ex8 point par point et tout par noyaux, CPU et GPU (11/10/2026)

Script : `cas/mesures-ex8.sh` (cas `ex8` de `cas/mesures-commun.sh`).
Exemple ex8 de mfem-mgis-examples (plaque thermomécanique : thermique avec
condition de Robin en faces de bord, mécanique avec gonflement, couplage
itératif), maillage d'origine (8424 hexaèdres), 5 pas, MUMPS pour la
mécanique, GMRES et BoomerAMG pour la thermique. Construction séparée
`build/examples-autres` liée à mfem-mgis 517e33cde4
(`install/mfem-mgis-autres`), mfem-mgis-examples 3a907f1. CPU : 2
processus d'un fil, comme le test de l'exemple ; GPU : un processus,
`--device cuda`, lois de l'exemple compilées par l'interface `generic`
(intégrées sur l'hôte). Deux passages (`resume.txt`, `repetition-2`) ;
au premier, la charge moyenne valait encore 6 (fin d'un essai précédent),
au second 1 ; les deux concordent à 2 % près.

Couplage sur les 5 pas, en secondes :

| Cas | Point par point | Tout par noyaux | Gain |
|---|---|---|---|
| CPU, 2 processus | 30,97 et 30,62 | 22,35 et 22,25 | ×1,38 à ×1,39 |
| GPU, 1 processus | 167,66 et 166,04 | 38,32 et 37,59 | ×4,4 |

Phases du premier passage (`cas/resume-profilage.py`, deux problèmes
cumulés) :

| Cas | Intégration | Résidu | Tangente | Factorisation (MUMPS) et préconditionneur | Résolution |
|---|---|---|---|---|---|
| CPU, point par point | 13,01 | 2,66 | 5,77 | 7,74 | 1,76 |
| CPU, noyaux | 10,19 | 0,33 | 2,26 | 7,75 | 1,80 |
| GPU, point par point | 66,87 | 48,47 | 36,33 | 11,94 | 4,00 |
| GPU, noyaux | 19,53 | 1,61 | 1,44 | 11,79 | 3,92 |

- Mêmes statistiques des champs (température, norme du déplacement,
  gonflement) dans les quatre cas, et mêmes nombres d'itérations de
  Newton (64 : 35 en thermique, 29 en mécanique, pour 114 évaluations du
  résidu).
- Sur CPU, l'intégration de la loi mécanique (U3Si2) et la factorisation
  de MUMPS dominent ; le résidu (×8) et la tangente (×2,6) par noyaux ne
  pèsent plus.
- Sur le GPU, le point par point est bien plus lent que sur CPU (résidu
  48,5 s contre 2,7 s, cause non analysée) ; tout par noyaux, il reste
  l'intégration des lois sur l'hôte (19,5 s) et MUMPS (11,8 s). Le levier suivant est de
  compiler les lois d'ex8 pour le GPU (interface `generic-parallel`).
- Le maillage raffiné une fois est hors de portée de cette mesure : plus
  de 20 min par pas avec MUMPS sur 8 processus, et 1 100 s de PCG avec
  BoomerAMG (options de systèmes) pour 5 pas sur 2 processus.
