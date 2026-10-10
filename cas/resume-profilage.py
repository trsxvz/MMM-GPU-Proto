#!/usr/bin/env python3
"""Résumé par phase des tableaux de profilage de mfem-mgis.

usage : cas/resume-profilage.py fichier...

Pour chaque fichier (sortie d'un cas lancé avec le profilage), somme le
temps de chaque section (le temps maximal sur les processus en parallèle),
toutes positions dans l'arbre confondues (les sections de même nom ne sont
pas imbriquées entre elles), et l'affiche par phase.
"""
import re
import sys

# phases : nom affiché, sections de mfem-mgis
PHASES = [
    ("total", ["root"]),
    ("intégration", ["NS::processNewUnknownsEstimate", "FPS::integrate"]),
    ("résidu", ["NS::computeResidual"]),
    ("tangente", ["NS::getJacobian"]),
    ("préconditionneur", ["NS::setLinearSolverOperator"]),
    ("résolution", ["MFEM::Mult(r,c)"]),
    ("accélération", ["FPS::accelerate"]),
    ("prédiction", ["NLEPIB::computePrediction"]),
    ("fin de pas", ["NLEPIB::update"]),
]
# tableau séquentiel : appels, temps, part ; parallèle : appels, min, moyenne,
# max, part, déséquilibre
LINE = re.compile(r"^\s*\|[\s|]*(?:-->|>)\s*(\S.*?)\s*\|\s*(\d+)\s*\|"
                  r"((?:\s*[\d.]+%?\s*\|)+)\s*$")


def summarize(path):
    times = {}
    calls = {}
    with open(path, errors="replace") as f:
        for line in f:
            m = LINE.match(line)
            if m is None:
                continue
            name = m.group(1)
            values = [v.strip() for v in m.group(3).split("|") if v.strip()]
            # temps maximal en parallèle (4e colonne), temps en séquentiel
            t = float(values[2] if len(values) > 3 else values[0])
            times[name] = times.get(name, 0.0) + t
            calls[name] = calls.get(name, 0) + int(m.group(2))
    return times, calls


def main():
    header = ["cas"] + [p for p, _ in PHASES]
    print(" | ".join(header))
    print(" | ".join(["---"] * len(header)))
    for path in sys.argv[1:]:
        times, calls = summarize(path)
        if not times:
            continue
        row = [re.sub(r"^.*/|\.txt$", "", path)]
        for _, sections in PHASES:
            t = sum(times.get(s, 0.0) for s in sections)
            n = sum(calls.get(s, 0) for s in sections)
            row.append(f"{t:.2f} ({n})" if n else "-")
        print(" | ".join(row))


if __name__ == "__main__":
    main()
