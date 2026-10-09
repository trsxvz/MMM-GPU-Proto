/*!
 * \file   src/MaterialArrays.cxx
 * \brief  This file implements the `MaterialArrays` class
 * \date   09/10/2026
 */

#include <limits>
#include <algorithm>
#include "MGIS/Raise.hxx"
#include "MGIS/MemoryLocation.hxx"
#include "MGIS/FieldView.hxx"
#include "MGIS/Behaviour/Behaviour.hxx"
#include "MFEMMGIS/MaterialArrays.hxx"

namespace mfem_mgis {

  /*!
   * \brief allocates an array in the memory of MFEM, set to zero on the host
   * \param[out] v: array
   * \param[in] s: size of the array
   */
  static void allocateMaterialArray(mfem::Vector& v, const size_type s) {
    if (s > static_cast<size_type>(std::numeric_limits<int>::max())) {
      mgis::raise("MaterialArrays: array too large for MFEM");
    }
    v.SetSize(static_cast<int>(s));
    v.UseDevice(true);
    if (s != 0) {
      std::fill_n(v.HostWrite(), s, real{0});
    }
  }  // end of allocateMaterialArray

  /*!
   * \return a view of the copy of the array on the host, made valid, the
   * copy on the device being invalidated
   * \param[in] v: array
   */
  static std::span<real> getHostView(mfem::Vector& v) {
    if (v.Size() == 0) {
      return {};
    }
    return {v.HostReadWrite(), static_cast<std::size_t>(v.Size())};
  }  // end of getHostView

  MaterialArrays::MaterialArrays(const Behaviour& b,
                                 const size_type n,
                                 const mgis::FieldLayout l)
      : array_layout(l) {
    using mgis::behaviour::getArraySize;
    const auto gs = getArraySize(b.gradients, b.hypothesis);
    const auto ts = getArraySize(b.thermodynamic_forces, b.hypothesis);
    const auto is = getArraySize(b.isvs, b.hypothesis);
    for (std::size_t i = 0; i != 2; ++i) {
      allocateMaterialArray(this->gradients[i], n * gs);
      allocateMaterialArray(this->thermodynamic_forces[i], n * ts);
      allocateMaterialArray(this->internal_state_variables[i], n * is);
      if (b.computesStoredEnergy) {
        allocateMaterialArray(this->stored_energies[i], n);
      }
      if (b.computesDissipatedEnergy) {
        allocateMaterialArray(this->dissipated_energies[i], n);
      }
    }
    allocateMaterialArray(this->tangent_operator_blocks,
                          n * mgis::behaviour::getTangentOperatorArraySize(b));
    // the deformation gradients are initialized to the identity, as MGIS
    // does for the arrays it allocates (the first three components are the
    // diagonal ones)
    if ((b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) &&
        (b.kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY) &&
        (n != 0)) {
      const auto pbp = l == mgis::FieldLayout::POINT_BY_POINT;
      const auto ps = pbp ? gs : size_type{1};
      const auto cs = pbp ? size_type{1} : n;
      for (auto& g : this->gradients) {
        auto* const F = g.HostReadWrite();
        for (size_type i = 0; i != n; ++i) {
          for (size_type c = 0; c != 3; ++c) {
            F[i * ps + c * cs] = real{1};
          }
        }
      }
    }
  }  // end of MaterialArrays

  mgis::behaviour::MaterialDataManagerInitializer
  MaterialArrays::getInitializer() {
    auto i = mgis::behaviour::MaterialDataManagerInitializer{};
    auto state = [this](mgis::behaviour::MaterialStateManagerInitializer& si,
                        const std::size_t s) {
      si.gradients = getHostView(this->gradients[s]);
      si.thermodynamic_forces = getHostView(this->thermodynamic_forces[s]);
      si.internal_state_variables =
          getHostView(this->internal_state_variables[s]);
      si.stored_energies = getHostView(this->stored_energies[s]);
      si.dissipated_energies = getHostView(this->dissipated_energies[s]);
    };
    state(i.s0, 0);
    state(i.s1, 1);
    i.K = getHostView(this->tangent_operator_blocks);
    i.location = mgis::MemoryLocation::HOST;
    i.layout = this->array_layout;
    return i;
  }  // end of getInitializer

  void MaterialArrays::synchronizeOnHost() {
    for (std::size_t i = 0; i != 2; ++i) {
      static_cast<void>(getHostView(this->gradients[i]));
      static_cast<void>(getHostView(this->thermodynamic_forces[i]));
      static_cast<void>(getHostView(this->internal_state_variables[i]));
      static_cast<void>(getHostView(this->stored_energies[i]));
      static_cast<void>(getHostView(this->dissipated_energies[i]));
    }
    static_cast<void>(getHostView(this->tangent_operator_blocks));
  }  // end of synchronizeOnHost

  void MaterialArrays::synchronizeOnHostForReading() const {
    auto read = [](const mfem::Vector& v) {
      if (v.Size() != 0) {
        static_cast<void>(v.HostRead());
      }
    };
    for (std::size_t i = 0; i != 2; ++i) {
      read(this->gradients[i]);
      read(this->thermodynamic_forces[i]);
      read(this->internal_state_variables[i]);
      read(this->stored_energies[i]);
      read(this->dissipated_energies[i]);
    }
    read(this->tangent_operator_blocks);
  }  // end of synchronizeOnHostForReading

  MaterialArrays::~MaterialArrays() = default;

}  // end of namespace mfem_mgis
