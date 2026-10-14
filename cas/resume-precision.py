#!/usr/bin/env python3
"""Résumé des mesures float contre double de cas/mesures-precision-poutre.sh.

Pour chaque calcul (précision et critère de Newton) : temps minimal et
médian des répétitions, itérations de Newton, norme de la solution et son
écart relatif à celle du double avec le même critère, puis le rapport des
temps double / float.

usage : cas/resume-precision.py <dossier de résultats>
"""
import re
import statistics
import sys
from pathlib import Path

LINE = re.compile(r"^dofs (\d+) \| nls (\S+) \| ls (\S+) \| \|\|u\|\| (\S+) \| "
                  r"iterations per step:([\d ]+)\| total (\d+)")
TIME = re.compile(r"^time ([0-9.]+) s")
NAME = re.compile(r"^(double|float)-(.+)-(\d+)\.txt$")


def read(path):
    """lit une sortie de la poutre : inconnues, norme, itérations, temps"""
    dofs = norm = total = time = None
    for line in path.read_text(errors="replace").splitlines():
        m = LINE.match(line)
        if m:
            dofs, norm, total = int(m.group(1)), float(m.group(4)), int(m.group(6))
        m = TIME.match(line)
        if m:
            time = float(m.group(1))
    return dofs, norm, total, time


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    directory = Path(sys.argv[1])
    runs = {}
    for path in sorted(directory.glob("*.txt")):
        m = NAME.match(path.name)
        if not m:
            continue
        dofs, norm, total, time = read(path)
        if time is None:
            print(f"{path.name} : pas de temps (calcul en échec ?)")
            continue
        runs.setdefault((m.group(2), m.group(1)), []).append(
            (time, norm, total, dofs))
    criteria = sorted({c for c, _ in runs})
    print(f"{'calcul':38s} {'n':>2s} {'min (s)':>9s} {'médiane':>9s} "
          f"{'itér.':>6s} {'||u||':>14s} {'écart':>9s}")
    for criterion in criteria:
        reference = runs.get((criterion, "double"))
        for precision in ("double", "float"):
            values = runs.get((criterion, precision))
            if not values:
                continue
            times = [v[0] for v in values]
            norm, total = values[0][1], values[0][2]
            gap = ""
            if reference and precision == "float":
                gap = f"{abs(norm - reference[0][1]) / reference[0][1]:.2e}"
            print(f"{precision + '-' + criterion:38s} {len(values):2d} "
                  f"{min(times):9.3f} {statistics.median(times):9.3f} "
                  f"{total:6d} {norm:14.10g} {gap:>9s}")
        if (criterion, "double") in runs and (criterion, "float") in runs:
            d = min(v[0] for v in runs[(criterion, "double")])
            f = min(v[0] for v in runs[(criterion, "float")])
            print(f"{'  double / float (temps minimaux)':38s} {d / f:9.2f}")


if __name__ == "__main__":
    main()
