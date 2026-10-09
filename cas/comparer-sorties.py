#!/usr/bin/env python3
"""Comparaison de deux fichiers de résultats numériques.

usage : cas/comparer-sorties.py référence fichier [--rtol r] [--atol a]

Compare un à un les nombres des deux fichiers (lus comme des suites de
nombres séparés par des blancs) : la différence de chaque valeur doit être
inférieure à atol + rtol * |référence|. Affiche la plus grande différence
relative ; code de retour non nul si les fichiers diffèrent.
"""
import argparse
import sys


def read(path):
    with open(path) as f:
        return [float(v) for v in f.read().split()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("file")
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=0.0)
    args = parser.parse_args()
    a = read(args.reference)
    b = read(args.file)
    if len(a) != len(b) or not a:
        print(f"nombres de valeurs différents ou nuls : {len(a)}, {len(b)}")
        return 1
    worst = 0.0
    ok = True
    for x, y in zip(a, b):
        d = abs(x - y)
        if d > args.atol + args.rtol * abs(x):
            ok = False
        if x != 0:
            worst = max(worst, d / abs(x))
    print(f"{len(a)} valeurs, plus grande différence relative : {worst:.3e}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
