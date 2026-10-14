# Coût des post-traitements du polycristal, lot L8 (13/10/2026)

Polycristal de mm-opera-hpc (5 grains, monocristal de l'UO2 en grandes
transformations, 20 553 inconnues), 10 pas, tous les calculs par noyaux,
Newton, sur le GPU (un processus) et sur 8 processus CPU, sans puis avec
les post-traitements du cas (export ParaView des résultats et moyennes des
contraintes à chaque pas), calculés sur l'hôte puis par noyaux
(`MFEM_MGIS_POST_PROCESSINGS_COMPUTATION=Kernels`). mfem-mgis 54ff8e6c8a
(installation séparée `install/mfem-mgis-54ff8e6c8a`, construction
`build/mm-opera-hpc-54ff8e6c8a`). Script : `cases/measure-post-processings.sh`.

Temps total (secondes), et temps des post-traitements
(`NLEP::PP::execute`) :

| | Sans | Sur l'hôte | Par noyaux |
|---|---|---|---|
| GPU, total | 11,00 | 12,48 | 12,53 |
| GPU, post-traitements | | 1,48 | 1,53 |
| 8 processus, total | 28,68 | 30,28 | 30,02 |
| 8 processus, post-traitements (moyenne) | | 0,19 | 0,18 |

- Les post-traitements coûtent 1,5 s sur le GPU (12 % du calcul), et ne
  gagnent rien à être calculés par noyaux : ce cas fait surtout l'export
  ParaView (écriture des fichiers sur l'hôte), alors que la bulle, où les
  noyaux divisent le temps total par 2,8 (`results/post-processings-kernels-2026-10-11`),
  fait surtout des projections L2 et des valeurs aux nœuds.
- La recopie des matériaux vers l'hôte avant les post-traitements
  (`NLEP::PP::synchronizeOnHost`) coûte 0,02 s.
- Mêmes contraintes macroscopiques dans toutes les variantes (écart relatif
  au plus 2,4e-13 sur le GPU, nul sur 8 processus).
