# Check that the `generic-parallel` interface refuses a material property
# whose name contains `_generic_parallel`, the suffix of the generated
# functions, whereas the `generic` interface accepts it.
#
# usage: sh ../check-refused-name.sh mfront configuration_file
#
# The sources are generated in the `refused-name` directory.
set -e
mfront=$1
configuration=$2
case "$mfront" in
  /*) ;;
  */*) mfront=$PWD/$mfront ;;
esac
case "$configuration" in
  /*) ;;
  *) configuration=$PWD/$configuration ;;
esac
directory=$(cd "$(dirname "$0")" && pwd)
file=$directory/InvalidNameTest_generic_parallel.mfront
mkdir -p refused-name
cd refused-name
"$mfront" --interface=generic "$file" > generic.log 2>&1
if "$mfront" --interface=generic-parallel \
  --configuration-file="$configuration" "$file" > generic-parallel.log 2>&1; then
  echo "the name of the material property has not been refused" >&2
  exit 1
fi
grep -q "must not contain '_generic_parallel'" generic-parallel.log
