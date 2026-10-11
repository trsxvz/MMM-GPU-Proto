# ex8 : solveur de la thermique sur le GPU (11/10/2026)

Exemple ex8 de mfem-mgis-examples (plaque thermomécanique), maillage
d'origine (9 516 inconnues pour la thermique, 28 548 pour la
mécanique), 5 pas, tout par noyaux. Mécanique : MUMPS avec
`-msMc SymmetricPositiveDefinite` (`results/ex8-mumps-symetrique-2026-10-11`).
Construction `build/examples-post-ex8`, liée à `install/mfem-mgis-post`.

Après le passage de la mécanique à MUMPS symétrique, la thermique est le
second poste d'ex8 sur le GPU. Variantes du solveur de la thermique
(`-lsTh`, `-pcTh`) : GMRES de hypre (défaut) avec BoomerAMG (défaut),
DiagScale, ILU ou ParaSails, et MUMPS sur l'hôte.

## GPU, 1 processus, lois cuda (un passage)

| Thermique | Total (s) | Thermique (s) | Résolutions (35) | Préparation |
|---|---|---|---|---|
| GMRES et BoomerAMG (défaut) | 12,62 | 4,26 | 3,49 | 0,003 |
| GMRES et DiagScale | 12,80 | 4,56 | 3,83 | 0,000 |
| MUMPS (`-lsTh MUMPSSolver`) | **11,50** | **3,03** | 0,06 | 2,17 (factorisations) |
| GMRES et ILU | arrêté à 600 s | | | |
| GMRES et ParaSails | plantage | | | |

La mécanique ne change pas (8,2 à 8,4 s, dont 7,3 s de factorisations).
Avec MUMPS, BoomerAMG ou DiagScale, mêmes itérations de Newton et mêmes
statistiques des champs.

- MUMPS sur l'hôte est le plus rapide pour la thermique d'ex8 sur le
  GPU : 4,26 à 3,03 s pour la thermique, 12,6 à 11,5 s pour ex8. Les 35
  factorisations (2,17 s) coûtent moins que les 35 résolutions de GMRES
  sur le GPU (3,49 s).
- Sur le GPU, une résolution de GMRES et BoomerAMG prend 100 ms pour
  9 516 inconnues, contre 42 ms sur 2 processus du CPU. DiagScale, sans
  préparation, met presque autant de temps que BoomerAMG. Le système est
  probablement trop petit pour le GPU : le temps tiendrait surtout au
  coût fixe de chaque itération (lancements de noyaux, réductions), le
  nombre d'itérations n'ayant pas été relevé.
- ILU : sur le GPU, mfem-mgis impose ILU(0) (les niveaux de remplissage
  supérieurs plantent dans hypre, `reproducers/hypre-ilu-gpu`). GMRES ne
  converge pas dès la deuxième itération de Newton du premier pas. Le pas
  de temps est alors redécoupé, et le calcul n'est qu'à 12 700 s sur
  20 000 s du premier pas après 600 s. Sur le CPU (2 processus, ILU(1)
  par blocs, défaut de MFEM), ILU converge (tableau suivant). Ces essais
  ne disent pas si l'échec vient du niveau 0 ou du GPU.
- ParaSails : deux défauts, corrigés dans mfem-mgis (branche `parasails`,
  607f895a0d, détails plus bas).

## CPU, 2 processus, lois generic (un passage, cœurs 0 à 3)

| Thermique | Thermique (s) | Résolutions (35) | Itérations de Newton |
|---|---|---|---|
| GMRES et BoomerAMG (défaut) | 6,04 | 1,49 | 114 |
| GMRES et ILU(1) par blocs | 8,09 | 3,53 | 114 |
| GMRES et ParaSails, avant correction | plantage à la 2e résolution | | |
| GMRES et ParaSails, après correction | échec : pas de temps redécoupés jusqu'à la limite | | |

Les totaux (28 à 30 s) ne sont pas comparables aux 19,9 s de
`results/ex8-mumps-symetrique-2026-10-11` : les deux processus étaient
liés aux cœurs 0 à 3, sur une machine partagée avec une autre session.
Les résolutions de la thermique (1,49 s) et les factorisations de la
mécanique (5,32 s) sont les mêmes que dans cette mesure. Seules les
colonnes de la thermique servent à comparer les préconditionneurs.

## ParaSails

1. Sur le CPU comme sur le GPU, `mfem::HypreParaSails::SetOperator`
   (MFEM 4.10, aussi sur la branche principale de MFEM) crée le nouveau
   préconditionneur de hypre avec le communicateur de l'opérateur
   précédent. Les solveurs non linéaires de mfem-mgis détruisent cet
   opérateur avant de donner le nouveau. Le communicateur lu est alors
   invalide, et la deuxième résolution plante (`MPI_ERR_COMM` ou erreur de
   segmentation). Reproducteur : `reproducers/hypre-parasails-operator`.
   `HypreEuclid` lit le communicateur du nouvel opérateur et n'a pas ce
   défaut. Correction dans mfem-mgis : le préconditionneur de la fabrique
   remplace l'opérateur par le nouveau avant d'appeler cette méthode.
   Correction amont possible : une ligne de MFEM (`*A` remplacé par
   `*new_A`).
2. ParaSails ne calcule que sur l'hôte. Quand hypre calcule sur le GPU, sa
   préparation lit la matrice rangée sur le GPU et plante dans
   `MatrixSetRow`, dès la première préparation. mfem-mgis le refuse
   maintenant avec un message (`gpu-gmres-parasails-refuse.txt`).

Après correction, ParaSails (options par défaut de MFEM) reste un mauvais
préconditionneur pour la thermique d'ex8. GMRES met 1 944 à 5 626
itérations quand il converge, et sinon s'arrête à 10 000 itérations sur
un résidu relatif de 0,5 à 0,7 (`np2-gmres-parasails-v1.txt`). Le pas de
temps est redécoupé jusqu'à la limite du nombre d'échecs.

## Fichiers

- `gpu-*.txt` : GPU, un processus (`gpu-gmres-parasails.txt` avant la
  correction, `gpu-gmres-parasails-refuse.txt` après).
- `np2-*.txt` : CPU, 2 processus (`np2-gmres-parasails.txt` avant la
  correction, `np2-gmres-parasails-corrige.txt` après,
  `np2-gmres-parasails-v1.txt` après avec `-v 1`, arrêté à 40 s).

Commande (GPU) :

```
./rjh_plate -m assemblage_hexa.msh -r 0 -ns 5 -no-pp --device cuda \
  -lU .../libGenericParallelU3SI2-cuda.so \
  -lA .../libGenericParallelALFENI-cuda.so \
  -msMc SymmetricPositiveDefinite [-pcTh HypreDiagScale|HypreILU|HypreParaSails] \
  [-lsTh MUMPSSolver]
```
