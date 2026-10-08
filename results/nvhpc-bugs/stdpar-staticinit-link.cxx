// nvc++ -stdpar=gpu: some constant data of the device code are exported
// under names which are only unique in a source file
// (`__staticinit_<counter>_<function>`). Two source files giving the same
// counter to the same inline function can't be linked together:
// "nvlink error : Multiple definition of '__staticinit_...'". Two copies of
// this file, which only differ by the name of the exported function, are
// enough. Without relocatable device code (`-gpu=nordc`), the device code
// of each file is linked on its own and the library is built.
// for n in 1 2; do
//   nvc++ -O2 -stdpar=gpu -gpu=mem:managed -fPIC -DNAME=f$n \
//     -c stdpar-staticinit-link.cxx -o f$n.o
// done
// nvc++ -shared -stdpar=gpu -gpu=mem:managed f1.o f2.o -o libf.so
#include <algorithm>
#include <execution>

extern "C" void NAME(double* const values, const unsigned n) {
  std::for_each_n(std::execution::par_unseq, values, n,
                  [](double& v) { v += 1; });
}
