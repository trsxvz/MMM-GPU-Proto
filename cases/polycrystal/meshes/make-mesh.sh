#!/bin/bash
# Periodic RVE of N grains for the polycrystal of mm-opera-hpc: Laguerre
# tessellation of Merope (germs by saturated RSA, as mesh/1000crystals.py of
# mm-opera-hpc), mean volume of a grain 1000 (cube of side (1000 N)^(1/3)),
# mesh size 2 (about 6 elements across a grain), meshed by gmsh (msh 2.2),
# orientations by randomVectorGeneration of mm-opera-hpc (algorithm
# `angles`, seed N of numpy), then the short edges collapsed by
# clean_mesh.py (see README.md).
#
# usage: make-mesh.sh N [seed]   (in the directory of the meshes)
# Variables: MEROPE_DIR (sources of Merope, built in INSTALL-DIR), PYTHON
# (python with numpy, numba, mkl, gmsh; default python3), GMSH (default gmsh),
# MM_OPERA_HPC (sources of mm-opera-hpc, default the submodule).
set -eu
N=$1
seed=${2:-1}
here=$(cd "$(dirname "$0")" && pwd)
MM_OPERA_HPC=${MM_OPERA_HPC:-$here/../../../mm-opera-hpc}
PYTHON=${PYTHON:-python3}
GMSH=${GMSH:-gmsh}
export PYTHONPATH=$MEROPE_DIR/INSTALL-DIR/lib:$MEROPE_DIR/tools/python
mkdir -p "v$N"
cd "v$N"
L=$("$PYTHON" -c "print(($N*1000)**(1/3))")
sed -e "s/L = \[100, 100, 100\]/L = [$L, $L, $L]/" \
  -e "s/nbSpheres = 1000 /nbSpheres = $N/" -e 's/MeshSize = 1$/MeshSize = 2/' \
  -e "s/randomSeed = 1$/randomSeed = $seed/" \
  -e "s/1000crystals.geo/${N}crystals.geo/" \
  "$MM_OPERA_HPC/polycrystal/mesh/1000crystals.py" > "${N}crystals.py"
"$PYTHON" "${N}crystals.py" > geo.log 2>&1
"$GMSH" -3 -nt 8 -format msh22 "${N}crystals.geo" -o "${N}crystals-raw.msh" \
  > gmsh.log 2>&1
if grep -qE '^Error|No elements in volume' gmsh.log; then
  echo "gmsh errors (see v$N/gmsh.log): try another seed" >&2
  exit 1
fi
"$PYTHON" -c "
import sys
sys.path.insert(0, '$MM_OPERA_HPC/polycrystal/randomVectorGeneration')
import numpy as np
np.random.seed($N)
import generate_random_vectors as g
g.main(nbVec=2 * $N, algoKey='angles')" > /dev/null
mv vectors.txt "vectors_${N}crystals.txt"
"$PYTHON" -I "$here/clean_mesh.py" "${N}crystals-raw.msh" 0.3 \
  "${N}crystals.msh" | tee clean.log
