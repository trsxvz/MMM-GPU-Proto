#pragma once

// constexpr: also callable from the code of the behaviours generated for a
// GPU (nvcc with --expt-relaxed-constexpr)
constexpr double computeHardeningModulus(const double T_actuelle) noexcept {
  double H = 1500e6;

  if (T_actuelle > 573.15) {
    H -= 1.5e6 * (T_actuelle - 573.15);
  }

  if (H < 100e6) {
    H = 100e6;
  }

  return H;
}