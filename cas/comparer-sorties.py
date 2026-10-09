#!/usr/bin/env python3
"""Comparaison de deux fichiers de résultats numériques.

usage : cas/comparer-sorties.py référence fichier [--rtol r] [--atol a]

Compare un à un les nombres des deux fichiers, rangés en lignes de
colonnes (par exemple : temps, composantes du gradient et des contraintes
macroscopiques). La différence de chaque valeur doit être inférieure à
atol + rtol * m, m étant la plus grande valeur absolue de sa colonne dans
la référence : une composante proche de zéro (une contrainte imposée
nulle) n'est pas comparée à elle-même. Affiche la plus grande différence
ainsi rapportée ; code de retour non nul si les fichiers diffèrent.
"""
import argparse
import sys


def read(path):
    with open(path) as f:
        return [[float(v) for v in line.split()] for line in f if line.split()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("file")
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=0.0)
    args = parser.parse_args()
    a = read(args.reference)
    b = read(args.file)
    if len(a) != len(b) or not a or any(len(x) != len(y) for x, y in zip(a, b)):
        print("formes différentes ou fichier vide")
        return 1
    ncolumns = max(len(x) for x in a)
    scale = [0.0] * ncolumns
    for x in a:
        for j, v in enumerate(x):
            scale[j] = max(scale[j], abs(v))
    worst = 0.0
    ok = True
    n = 0
    for x, y in zip(a, b):
        for j, (u, v) in enumerate(zip(x, y)):
            n += 1
            d = abs(u - v)
            if d > args.atol + args.rtol * scale[j]:
                ok = False
            if scale[j] != 0:
                worst = max(worst, d / scale[j])
    print(f"{n} valeurs, plus grande différence rapportée au maximum de sa "
          f"colonne : {worst:.3e}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
