#!/bin/sh
# Checks the `execution_policy` option of the raja backend, which is copied
# in the generated code (generation of the sources only): a valid policy
# must be found in the generated sources, with the identifier of the
# execution policy of the backend exported between parentheses, an invalid
# one must be refused by mfront. tfel-check only tolerates the failure of a
# command, it does not require it: this script returns an error if an
# invalid policy is accepted.
#
# usage: sh check-execution-policies.sh mfront
mfront=${1:?usage: check-execution-policies.sh mfront}
source=src/ImplicitNorton-generic-parallel-raja.cxx
status=0

# generates the sources of the test behaviour with the given execution
# policy
generate() {
  {
    echo 'interfaces_options: {'
    echo '  generic-parallel: {'
    printf '    backend: {raja: {execution_policy: "%s"}}\n' "$1"
    echo '  }'
    echo '}'
  } > execution-policy.json
  rm -rf include src
  "$mfront" --interface=generic-parallel \
    --configuration-file=execution-policy.json \
    ../ImplicitNorton.mfront > /dev/null 2>&1
}

for policy in 'RAJA::cuda_exec<128, true>' \
  'RAJA::seq_exec' \
  '::RAJA::omp_parallel_for_exec' \
  'mfront::generic_parallel::raja::DefaultHostExecutionPolicy' \
  'RAJA::cuda_exec_explicit<256, 2, false>'; do
  if ! generate "$policy"; then
    echo "valid execution policy refused: '$policy'"
    status=1
  elif ! grep -qF "$policy" "$source"; then
    echo "valid execution policy not found in $source: '$policy'"
    status=1
  elif ! grep -qF '_generic_parallel_execution_policy, (' "$source"; then
    # the type of the execution policy of the backend contains the policy of
    # RAJA, whose commas must not split the arguments of the macro exporting
    # the identifier of the policy
    echo "identifier of the execution policy not protected in" \
      "$source: '$policy'"
    status=1
  fi
done

for policy in 'RAJA::cuda_exec<128' \
  'RAJA::cuda_exec<128>>' \
  'int); system(' \
  'A<B)>' \
  'RAJA::cuda_exec<128);abort();//>' \
  'RAJA::seq_exec, RAJA::seq_exec' \
  'RAJA:seq_exec' \
  'RAJA::cuda_exec<-1>' \
  ''; do
  if generate "$policy"; then
    echo "invalid execution policy accepted: '$policy'"
    status=1
  fi
done

rm -f execution-policy.json
rm -rf include src
exit $status
