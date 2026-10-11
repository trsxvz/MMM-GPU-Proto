# ex8 : BoomerAMG avec les options d'élasticité (11/10/2026)

Exemple ex8 de mfem-mgis-examples (plaque thermomécanique), maillage
d'origine, 5 pas, 2 processus du CPU, tout par noyaux pour la mécanique.
Construction `build/examples-post-ex8` (ex8 avec les options `-oMc` et
`-amgMc`, liée à `install/mfem-mgis-post`).

| Solveur de la mécanique | Couplage sur les 5 pas (s) |
|---|---|
| MUMPS (référence) | 23,52 |
| PCG + BoomerAMG, options d'élasticité, rangement byVDIM (`-lsMc HyprePCG -pcMc HypreBoomerAMG -oMc byVDIM -amgMc Elasticity`) | 201,16 |

- Mêmes statistiques des champs (températures et normes des
  déplacements) dans les deux calculs.
- Les options d'élasticité (interpolation par modes rigides) font pire
  que les options de systèmes essayées auparavant (environ 240
  itérations, 93 s, mesure du 11/10 après-midi notée dans `JOURNAL.md`,
  sans sortie brute conservée) : 441 itérations de PCG en moyenne par résolution
  (29 résolutions). La plaque est mince et ses matériaux contrastés :
  BoomerAMG ne convient pas tel quel (le premier niveau grossit peu :
  28 548 puis 22 482 lignes, complexité d'opérateur 3,6 ; l'explication
  par la géométrie et les matériaux n'est pas vérifiée). MUMPS reste le
  choix pour la mécanique d'ex8.
