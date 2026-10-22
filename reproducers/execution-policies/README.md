# Execution policies across libraries

Checks the design of the execution policies of the `generic-parallel`
interfaces on each platform, with stand-ins for the types of CUDA, SYCL and
Kokkos (`stand-ins/StandIns.hxx`), so that it builds without these libraries.

## Design checked

- The runtime library (`MFrontGenericParallelRuntime`) defines the base class
  `AbstractBackendExecutionPolicy`: virtual destructor and pure virtual
  `getIdentifier()`. On Windows, the class is exported by `dllexport`; on the
  other systems, all the symbols are hidden except this class, as with the
  flags given by `tfel-config`.
- The policies of the backends are defined in headers only
  (`include/MFront/GenericParallel`): `cuda` (stream), `stlpar` (common
  options only), `sycl` (holds a `sycl::queue`, which holds a shared pointer),
  `kokkos` (class template on the execution space, the stand-ins of the
  execution spaces having no attribute of visibility, as in Kokkos).
- The caller (stand-in for MGIS or a user code) is linked to the runtime
  library, creates the policies with the `make...` functions and designated
  initializers, and gives them through a `C` structure
  (`MFront/GenericParallel/BatchData.h`, opaque pointer in `C`).
- The generated library is loaded at run time (`dlopen` with `RTLD_NOW`, as
  MGIS, or `LoadLibrary`). Its functions compare the identifier of the policy
  with the one of their backend, then use `static_cast`, and copy what they
  need. It is not linked to the runtime library (on Linux and FreeBSD with
  `--no-undefined`, on macOS and Windows unresolved symbols are errors by
  default): it needs no symbol of this library.

## Tests

| Test | Caller | Generated library |
|---|---|---|
| `linked` | program linked to the caller | not linked to the runtime library |
| `linked-generated-linked` | program linked to the caller | linked to the runtime library |
| `loaded` | program loading the caller at run time, as Python loads a module | not linked to the runtime library |
| `c` | `C` program, null policies | not linked to the runtime library |
| `cross` | program linked to the caller, `cuda` and `stlpar` only | built by another toolchain (`EXECUTION_POLICIES_OTHER_GENERATED_LIBRARY`) |

For each backend: default policy for a null pointer (refused for `sycl`, whose
queue is required), identifier and members read, policies of the other
backends and of another execution space refused (`-5`), asynchronous call
refused without statistics (`-5`), copy of the queue destroyed by the
generated function. The tests `linked` and `loaded` also print, for
information, the answers of `dynamic_cast`, which the design does not use
(library `dynamic-cast`, linked to the runtime library).

`run.sh [build directory] [cmake options...]` builds and runs the tests for
both orders of declaration of `getIdentifier` and of the destructor
(`EXECUTION_POLICIES_IDENTIFIER_FIRST`) and for hidden and default
visibility (`EXECUTION_POLICIES_HIDDEN_VISIBILITY`). The workflow
`.github/workflows/execution-policies.yml` runs it on Linux (gcc and
libstdc++, clang and libc++), macOS, Windows (Visual Studio and MinGW) and
FreeBSD, then, on Windows, runs the test `cross` from a program built by
Visual Studio with the generated library built by MinGW, and the reverse.
