// Loading a library compiled by nvc++ from a program which does not load
// the math library fails: `libnvcpumath.so` (NVHPC 25.11) uses the `sincos`
// function of `libm.so.6` (an IFUNC symbol) but does not depend on it.
// usage: g++ dlopen-without-libm.cxx -o dlopen-without-libm -ldl
//        ./dlopen-without-libm lib.so (lib.so: any library built by nvc++)
// Expected: "lib.so: loaded". Obtained: "Relink `libnvcpumath.so' with
// `libm.so.6' for IFUNC symbol `sincos'" and the program stops. Programs
// using the math library (most C++ programs) are not affected.
#include <cstdio>
#include <dlfcn.h>
int main(int argc, char** argv) {
  for (int i = 1; i < argc; ++i) {
    auto* l = dlopen(argv[i], RTLD_NOW | RTLD_LOCAL);
    std::printf("%s: %s\n", argv[i], l == nullptr ? dlerror() : "loaded");
  }
  return 0;
}
