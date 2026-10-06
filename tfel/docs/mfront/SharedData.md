The `@SharedData` keyword declares a data shared by all the integration
points, built on the host by a `C++` expression. This expression is
evaluated once: the data is a constant member of the behaviour, which is a
reference to a copy of the result of the expression.

This keyword is followed by the type of the data, its name, the `=`
character and the expression, ended by a semicolon.

With the `generic-parallel` interface, the data is copied in the memory of
the device at the first call of a function of the generated library. Its
type must then be trivially copyable.

## Example

~~~~{.cpp}
@SharedData PCSlidingSystems gs =
    PCSlidingSystems::getPolyCrystalsSlidingSystems("PolyCrystalsAngles.txt");
~~~~
