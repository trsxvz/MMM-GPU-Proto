# AmgX contre BoomerAMG sur le polycristal (10/10/2026)

Polycristal de mm-opera-hpc : 5 grains, loi de monocristal de l'UO2,
problème périodique, traction uniaxiale, 20 553 inconnues, sur le GPU.
mfem-mgis de2e76a913, tous les noyaux. Script :
`cas/mesures-L5-polycristal-amgx.sh` (option
`--elastic-operator-linear-solver AmgX` du polycristal, paramètres
d'AmgX par `--amgx-parameters`, simple précision par
`--amgx-single-precision`). Une exécution par cas, charge 3.

NewtonK0 avec Eisenstat-Walker 1 (tolérance initiale 0,1), le multigrille
de l'opérateur élastique étant BoomerAMG ou AmgX. « Réglé » :
`strength_threshold=0.8, interpolator=MULTIPASS`, le meilleur réglage
de la poutre.

| Multigrille de K0 | Résolution | Intégrations | Résolutions linéaires |
|---|---|---|---|
| Newton (référence, PCG de hypre et BoomerAMG) | 11,66 s | | 98, 5,71 s |
| **BoomerAMG** | **9,04 s** | 225 | 164, 0,78 s |
| AmgX par défaut | 13,60 s | 274 | 213, 3,34 s |
| AmgX réglé | 10,77 s | 245 | 184, 1,41 s |
| AmgX réglé, simple | 10,70 s | 245 | 184, 1,19 s |

Contraintes macroscopiques identiques à celles de Newton (écart relatif
au plus 3,6e-6).

L'intégration de la loi domine (75 % du temps) : un multigrille moins
efficace ajoute des itérations de Newton, donc des intégrations, et AmgX
perd 1,7 s, alors que les résolutions linéaires elles-mêmes ne pèsent que
0,8 à 1,4 s. Sur ce cas, BoomerAMG reste le meilleur choix ; sur la
poutre, où la résolution linéaire domine, AmgX réglé en simple précision
gagne 12 à 15 % (results/L5-poutre-2026-10-10-2259).
