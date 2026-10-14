# Prints the hashes (/<hash>) of the packages of the active spack environment,
# concretized, except the ones whose names are given and the externals: all
# the nodes of the environment are printed, roots included, so that
# `spack install` builds every package but the given ones.
#
# Used by the job `cuda` of the continuous integration to build hypre and
# MFEM alone with a single job. `spack install --only dependencies` can't be
# used for that: spack 1.2 then installs neither the roots of the
# environment, even the ones the other packages depend on (openmpi, cmake),
# nor some of their dependencies.
# usage: spack -e <environment> python scripts/spack-specs-except.py name...
import sys

import spack.environment as ev

excluded = set(sys.argv[1:])
nodes = {}
for root in ev.active_environment().concrete_roots():
    for node in root.traverse():
        if (node.name not in excluded) and (not node.external):
            nodes[node.dag_hash()] = node
print(" ".join("/" + h for h in sorted(nodes)))
