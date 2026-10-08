// A class as generated for a model with a string parameter in device mode.
// MODEL_STATIC_MEMBER: static constexpr member (variant A); otherwise, local
// constant declared at the beginning of the function (variant B).
#include <string_view>
#include "constant-string.hxx"
#ifndef MODEL_STRING
#define MODEL_STRING std::string_view
#endif
#if defined(__CUDACC__) || defined(__HIPCC__)
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif
template <typename real>
struct Model {
#ifdef MODEL_STATIC_MEMBER
  static constexpr MODEL_STRING law = "linear";
#endif
  HOST_DEVICE real execute(const real x) const {
#ifndef MODEL_STATIC_MEMBER
    [[maybe_unused]] constexpr auto law = MODEL_STRING{"linear"};
#endif
    // run time comparisons and size, as in user code
    if (law == "quadratic") {
      return x * x;
    }
    return law.size() * x + (law == "linear" ? 1 : 0);
  }
};
