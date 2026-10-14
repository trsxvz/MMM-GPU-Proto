# Inverses des jacobiennes calculées une fois (11/10/2026)

Les noyaux n'utilisaient les jacobiennes que par leurs inverses, recalculées
par chaque fil à chaque point (pour chaque paire de nœuds dans le noyau des
matrices de rigidité) : 9 divisions et une trentaine de multiplications en
double précision, lente sur ce GPU. mfem-mgis (branche `jacobiennes`,
bf7673f4f8) : le noyau de la géométrie range les inverses à la place des
jacobiennes (même mémoire), les noyaux les lisent.

Banc `cases/assembly`, installations séparées `install/mfem-mgis-registres`
(4a007e85f2, avant) et `install/mfem-mgis-jacobiennes` (après), mêmes
commandes, un cœur (`taskset -c 4`) ou la RTX 5080, machine calme. Temps
médians de 10 assemblages. Sorties : `gpu.txt`, `host.txt`.

GPU (sans positions, chemin par défaut avec un GPU) :

| Cas | Avant | Après | Gain |
|---|---|---|---|
| grandes transformations, hexaèdres Q1 27 points (r5) | 31,8 ms | 19,0 ms | ×1,67 |
| petites déformations, hexaèdres Q1 27 points (r5) | 29,2 ms | 16,3 ms | ×1,79 |
| poutre, tétraèdres P1 (r5) | 18,1 ms | 8,1 ms | ×2,2 |
| grandes transformations, tétraèdres P2 (r4) | 52,8 ms | 26,4 ms | ×2,0 |
| thermique, tétraèdres P1 (r5) | 13,8 ms | 2,3 ms | ×5,9 |

Écarts entre les deux versions de 1e-17 à 6e-17 pour des valeurs de 4e-3 à
0,27 : ordre des additions atomiques.

Un cœur (avec les positions) : 52,3 contre 51,5 ms (grandes
transformations), 75,6 contre 72,3 ms (poutre), 116,7 contre 109,7 ms
(thermique) ; valeurs identiques au bit près. Le chemin de l'hôte
n'inversait les jacobiennes qu'une fois par point.

Registres (cuobjdump) : matrices de rigidité de la thermique 3D 56 au lieu
de 84 (64 avant le lot des matrices sans former B) ; gradients 40 à 46 ;
le reste inchangé.
