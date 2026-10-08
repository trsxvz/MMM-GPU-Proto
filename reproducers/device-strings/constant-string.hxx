#include <cstddef>
#include <string_view>
#if defined(__CUDACC__) || defined(__HIPCC__)
#define CS_HOST_DEVICE __host__ __device__
#else
#define CS_HOST_DEVICE
#endif
struct ConstantString {
  template <std::size_t N>
  CS_HOST_DEVICE constexpr ConstantString(const char (&s)[N]) noexcept
      : value(s), length(N - 1) {}
  CS_HOST_DEVICE constexpr std::size_t size() const noexcept { return this->length; }
  CS_HOST_DEVICE constexpr const char* data() const noexcept { return this->value; }
  CS_HOST_DEVICE constexpr bool operator==(const ConstantString& o) const noexcept {
    if (this->length != o.length) {
      return false;
    }
    for (std::size_t i = 0; i != this->length; ++i) {
      if (this->value[i] != o.value[i]) {
        return false;
      }
    }
    return true;
  }
  operator std::string_view() const noexcept { return {this->value, this->length}; }
 private:
  const char* value;
  std::size_t length;
};
