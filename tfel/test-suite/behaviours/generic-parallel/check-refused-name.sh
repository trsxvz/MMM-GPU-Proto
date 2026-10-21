# Check that the `generic-parallel` interface refuses a behaviour and a
# model whose names contain `_generic_parallel`, the suffix of the
# generated functions, whereas the `generic` interface accepts them.
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
# the configuration file is looked for in the current directory, then in
# the configuration files installed with TFEL, as MFront does
find_configuration() {
  if [ -f "$1" ]; then
    echo "$PWD/$1"
    return
  fi
  prefix=$(cd "$(dirname "$(command -v "$mfront")")/.." && pwd)
  for d in "$prefix"/share/tfel*/mfront/configuration-files; do
    if [ -f "$d/$1" ]; then
      echo "$d/$1"
      return
    fi
  done
  echo "configuration file '$1' not found" >&2
  exit 1
}
configuration=$(find_configuration "$configuration")
directory=$(cd "$(dirname "$0")" && pwd)
mkdir -p refused-name
cd refused-name
for name in InvalidNameTest_generic_parallel \
  InvalidNameModelTest_generic_parallel; do
  file=$directory/$name.mfront
  "$mfront" --interface=generic "$file" > "$name-generic.log" 2>&1
  if "$mfront" --interface=generic-parallel \
    --configuration-file="$configuration" "$file" \
    > "$name-generic-parallel.log" 2>&1; then
    echo "the name '$name' has not been refused" >&2
    exit 1
  fi
  grep -q "must not contain '_generic_parallel'" "$name-generic-parallel.log"
done
