The `@SharedData` keyword declares a data shared by all the integration
points, built on the host by a `C++` expression. This expression is
evaluated once: the data is a constant member of the behaviour, which is a
reference to a copy of the result of the expression.

The expression is evaluated by a static member function of the behaviour
class, on the host. It can't use the material properties, the state
variables or the parameters of the behaviour.

This keyword is followed by the type of the data, its name, the `=`
character and the expression, ended by a semicolon.

With the `generic-parallel` interface, the data is copied in the memory of
the device at the first call of a function of the generated library. Its
type must then be trivially copyable. The following restrictions apply:

- shared data are not supported in behaviour variables
  (`@BehaviourVariable`, `@BehaviourVariableFactory`).
- the copies are never freed (one copy per device with the `cuda` and
  `hip` backends, one copy per device and context with the `sycl`
  backend, one copy on the heap with the `stlpar` and `openmp` backends),
  except with the `kokkos` backend, which frees its copy when Kokkos is
  finalized. With the `stlpar` and `openmp` backends, the heap must be
  accessible by the device when the computations are offloaded (unified
  or managed memory, as with `nvc++`).

## Example

~~~~{.cpp}
@Includes {
#include "TFEL/Material/PolyCrystalsSlidingSystems.hxx"
#include "TFEL/Material/MetallicCFCSlidingSystems.hxx"
}

@Members {
  using PCSlidingSystems =
      PolyCrystalsSlidingSystems<Np, MetallicCFCSlidingSystems<real>, real>;
}

@SharedData PCSlidingSystems gs =
    PCSlidingSystems::getPolyCrystalsSlidingSystems("PolyCrystalsAngles.txt");
~~~~

In this example, `PCSlidingSystems` is a type alias declared in the
`@Members` block and `Np` is the number of grains, declared by the
`@IntegerConstant` keyword.
