The `@Profiling` keyword is followed by a boolean. If true, several
high resolutions clocks will be introduced in the generated code to
profile performance bottlenecks. The total time spend in various
portions of the generated code will be stored and displayed when the
calling process exits. Profiling is disabled in device code, and in the
whole library with the compilers translating host and device code in a
single pass (`nvc++ -stdpar=gpu`, `nvc++ -mp=gpu`, generic compiler of
AdaptiveCpp).

## Example

~~~~{.cpp}
@Profiling true;
~~~~
