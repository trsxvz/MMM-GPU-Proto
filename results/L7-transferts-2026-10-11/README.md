# Transferts entre l'hôte et le GPU (lot L7.3, 11/10/2026)

Profil nsys (`cas/profil-transferts.sh`) de la thermique stationnaire r5
(279 873 inconnues, 5 pas, 18 itérations de Newton) et de la poutre r4
(111 843 inconnues, 10 pas, NewtonK0 et Eisenstat-Walker 1, BoomerAMG
agressif), tous les noyaux, loi `cuda`, mfem-mgis e03fd927f7. Copies
regroupées par sens et par taille ; leur origine est donnée par gdb
(points d'arrêt sur `mfem::CuMemcpyDtoH` et `mfem::CuMemcpyHtoD` filtrés
par taille, `gdb-copies.txt`, thermique sur un pas).

| Cas | Durée | Hôte vers GPU | GPU vers hôte |
|---|---|---|---|
| Thermique r5 | 1,46 s | 60 ms | 58 ms |
| Poutre r4 | 1,97 s | 33 ms | 30 ms |

Thermique : les copies représentent 8 % du temps.

| Copie | Nombre | Temps | Origine |
|---|---|---|---|
| états, 3 composantes par point, GPU vers hôte | 40 | 40 ms | fin de pas (`update`), tableaux d'état ramenés sur l'hôte |
| états, 3 composantes par point, hôte vers GPU | 30 | 30 ms | intégration suivante, tableaux renvoyés au GPU |
| température aux points, hôte vers GPU | 92 | 21 ms | intégration : variables externes en début et en fin de pas |
| température aux points, GPU vers hôte | 46 | 14 ms | valeurs des inconnues aux points recopiées dans la variable externe |

Poutre : contraintes et déformations (6 composantes) et variables
internes (7 composantes), 80 et 40 copies dans chaque sens, 56 ms en tout
(3 %), même origine (fin de pas puis intégration suivante).

## Décisions

- Fin de pas : la mise à jour est faite sur le GPU, puis les tableaux
  sont ramenés sur l'hôte et leur copie sur le GPU invalidée. C'est le
  contrat documenté de `MaterialArrays` : après une mise à jour, l'état
  peut être modifié directement sur l'hôte (cas « write after update » de
  `BatchIntegrationTest`). Une synchronisation en lecture seule, qui
  garderait la copie du GPU valide, économiserait le renvoi (2 % sur la
  thermique, 1,4 % sur la poutre) mais casse ce contrat : essai annulé
  (le test échoue sur les périphériques `debug` et `cuda`). Une option
  qui garderait les états sur le GPU au prix d'écritures sur l'hôte
  explicitement synchronisées reste possible.
- Température : la garder sur le GPU (source sur le GPU pour la variable
  externe de fin de pas, copie de début de pas gardée d'une itération à
  l'autre, copie sur l'hôte seulement à la fin du pas) économiserait au
  plus 2,3 % de la thermique, pour une gestion de validité de plus dans
  `BehaviourIntegratorBase`. Non fait pour l'instant.
