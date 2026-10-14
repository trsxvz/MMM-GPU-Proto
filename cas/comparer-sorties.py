#!/usr/bin/env python3
"""Comparaison de deux fichiers de résultats numériques.

usage : cas/comparer-sorties.py référence fichier [--rtol r] [--atol a]
        [--groupes g[:r]...] [--csv]

Compare un à un les nombres des deux fichiers, rangés en lignes de
colonnes (par exemple : temps, composantes du gradient et des contraintes
macroscopiques). La différence de chaque valeur doit être inférieure à
atol + rtol * m, m étant la plus grande valeur absolue de son groupe de
colonnes dans la référence : une composante proche de zéro (une
contrainte imposée nulle) n'est pas comparée à elle-même. Par défaut,
chaque colonne forme un groupe ; avec --groupes, des colonnes de même
nature partagent une échelle, par exemple `--groupes 0 1-3 4-6` pour le
temps, le gradient et les contraintes : une contrainte transverse nulle en
théorie (bruit numérique de 1e2 Pa sous un chargement de 1e7 Pa) est
alors comparée à l'échelle des contraintes. Un groupe peut avoir sa propre
tolérance relative (`1-3:1e-6`), à la place de --rtol. Une valeur non
finie (NaN, infini) dans l'un des fichiers est une différence. Avec --csv,
les valeurs sont séparées par des virgules et la première ligne (les noms
des colonnes) est ignorée. Affiche la plus grande différence ainsi
rapportée ; code de retour non nul si les fichiers diffèrent.
"""
import argparse
import math
import sys


def read(path, csv=False):
    with open(path) as f:
        lines = f.readlines()
    if csv:
        lines = [line.replace(",", " ") for line in lines[1:]]
    return [[float(v) for v in line.split()] for line in lines if line.split()]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    parser.add_argument("file")
    parser.add_argument("--rtol", type=float, default=1e-10)
    parser.add_argument("--atol", type=float, default=0.0)
    parser.add_argument("--groupes", nargs="*", default=[],
                        help="groupes de colonnes partageant une échelle "
                        "(indices à partir de 0, intervalles i-j), suivis "
                        "éventuellement de leur tolérance relative (i-j:r)")
    parser.add_argument("--csv", action="store_true",
                        help="valeurs séparées par des virgules, première "
                        "ligne ignorée")
    args = parser.parse_args()
    a = read(args.reference, args.csv)
    b = read(args.file, args.csv)
    if len(a) != len(b) or not a or any(len(x) != len(y) for x, y in zip(a, b)):
        print("formes différentes ou fichier vide")
        return 1
    ncolumns = max(len(x) for x in a)
    scale = [0.0] * ncolumns
    rtol = [args.rtol] * ncolumns
    for x in a:
        for j, v in enumerate(x):
            scale[j] = max(scale[j], abs(v))
    for g in args.groupes:
        interval, _, tolerance = g.partition(":")
        first, _, last = interval.partition("-")
        columns = range(int(first), int(last or first) + 1)
        if any(j >= ncolumns for j in columns):
            print(f"groupe {g} : colonne absente")
            return 1
        m = max(scale[j] for j in columns)
        for j in columns:
            scale[j] = m
            if tolerance:
                rtol[j] = float(tolerance)
    worst = 0.0
    ok = True
    n = 0
    for x, y in zip(a, b):
        for j, (u, v) in enumerate(zip(x, y)):
            n += 1
            d = abs(u - v)
            if (not (math.isfinite(u) and math.isfinite(v))) or \
                    d > args.atol + rtol[j] * scale[j]:
                ok = False
            if scale[j] != 0:
                worst = max(worst, d / scale[j])
    print(f"{n} valeurs, plus grande différence rapportée au maximum de son "
          f"groupe : {worst:.3e}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
