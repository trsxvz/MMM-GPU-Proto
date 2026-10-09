/*!
 * \file   src/MaterialArrays.cxx
 * \brief  This file implements the `MaterialArrays` class
 * \date   09/10/2026
 */

#include <limits>
#include <cstddef>
#include <algorithm>
#include "mfem/general/device.hpp"
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
  static void allocateMaterialArray(mfem::Vector& v, const std::size_t s) {
    if (s > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
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

  //! \return the index of a state
  static std::size_t getStateIndex(const TimeStepStage s) {
    return s == TimeStepStage::BEGINNING_OF_TIME_STEP ? 0 : 1;
  }  // end of getStateIndex

  //! \return the index of an array of a state
  static std::size_t getArrayIndex(const MaterialArrays::StateArray a) {
    return static_cast<std::size_t>(a);
  }  // end of getArrayIndex

  MaterialArrays::MaterialArrays(const Behaviour& b,
                                 const size_type n,
                                 const mgis::FieldLayout l)
      : array_layout(l) {
    using mgis::behaviour::getArraySize;
    // sizes computed with the type of the sizes of MGIS, checked against the
    // range of the sizes of MFEM by `allocateMaterialArray`
    const auto np = static_cast<std::size_t>(n);
    const auto gs =
        static_cast<std::size_t>(getArraySize(b.gradients, b.hypothesis));
    const auto ts = static_cast<std::size_t>(
        getArraySize(b.thermodynamic_forces, b.hypothesis));
    const auto is =
        static_cast<std::size_t>(getArraySize(b.isvs, b.hypothesis));
    for (std::size_t i = 0; i != 2; ++i) {
      auto& arrays = this->state_arrays;
      allocateMaterialArray(arrays[getArrayIndex(StateArray::GRADIENTS)][i],
                            np * gs);
      allocateMaterialArray(
          arrays[getArrayIndex(StateArray::THERMODYNAMIC_FORCES)][i], np * ts);
      allocateMaterialArray(
          arrays[getArrayIndex(StateArray::INTERNAL_STATE_VARIABLES)][i],
          np * is);
      if (b.computesStoredEnergy) {
        allocateMaterialArray(
            arrays[getArrayIndex(StateArray::STORED_ENERGIES)][i], np);
      }
      if (b.computesDissipatedEnergy) {
        allocateMaterialArray(
            arrays[getArrayIndex(StateArray::DISSIPATED_ENERGIES)][i], np);
      }
    }
    allocateMaterialArray(
        this->tangent_operator_blocks,
        np * static_cast<std::size_t>(
                 mgis::behaviour::getTangentOperatorArraySize(b)));
    // the deformation gradients are initialized to the identity, as MGIS
    // does for the arrays it allocates (the first three components are the
    // diagonal ones)
    if ((b.btype == Behaviour::STANDARDFINITESTRAINBEHAVIOUR) &&
        (b.kinematic == Behaviour::FINITESTRAINKINEMATIC_F_CAUCHY) &&
        (np != 0)) {
      const auto pbp = l == mgis::FieldLayout::POINT_BY_POINT;
      const auto ps = pbp ? gs : std::size_t{1};
      const auto cs = pbp ? std::size_t{1} : np;
      for (auto& g : this->state_arrays[getArrayIndex(StateArray::GRADIENTS)]) {
        auto* const F = g.HostReadWrite();
        for (std::size_t i = 0; i != np; ++i) {
          for (std::size_t c = 0; c != 3; ++c) {
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
                        const TimeStepStage s) {
      si.gradients =
          getHostView(this->getMutableArray(StateArray::GRADIENTS, s));
      si.thermodynamic_forces = getHostView(
          this->getMutableArray(StateArray::THERMODYNAMIC_FORCES, s));
      si.internal_state_variables = getHostView(
          this->getMutableArray(StateArray::INTERNAL_STATE_VARIABLES, s));
      si.stored_energies =
          getHostView(this->getMutableArray(StateArray::STORED_ENERGIES, s));
      si.dissipated_energies = getHostView(
          this->getMutableArray(StateArray::DISSIPATED_ENERGIES, s));
    };
    state(i.s0, TimeStepStage::BEGINNING_OF_TIME_STEP);
    state(i.s1, TimeStepStage::END_OF_TIME_STEP);
    i.K = getHostView(this->tangent_operator_blocks);
    i.location = mgis::MemoryLocation::HOST;
    i.layout = this->array_layout;
    return i;
  }  // end of getInitializer

  void MaterialArrays::synchronizeOnHost() {
    for (auto& arrays : this->state_arrays) {
      for (auto& v : arrays) {
        static_cast<void>(getHostView(v));
      }
    }
    static_cast<void>(getHostView(this->tangent_operator_blocks));
  }  // end of synchronizeOnHost

  void MaterialArrays::synchronizeOnHostForReading() const {
    auto read = [](const mfem::Vector& v) {
      if (v.Size() != 0) {
        static_cast<void>(v.HostRead());
      }
    };
    for (const auto& arrays : this->state_arrays) {
      for (const auto& v : arrays) {
        read(v);
      }
    }
    read(this->tangent_operator_blocks);
  }  // end of synchronizeOnHostForReading

  void MaterialArrays::synchronizeOnHost(const StateArray a,
                                         const TimeStepStage s) {
    static_cast<void>(getHostView(this->getMutableArray(a, s)));
  }  // end of synchronizeOnHost

  void MaterialArrays::synchronizeOnHostForReading(
      const StateArray a, const TimeStepStage s) const {
    const auto& v = this->getArray(a, s);
    if (v.Size() != 0) {
      static_cast<void>(v.HostRead());
    }
  }  // end of synchronizeOnHostForReading

  void MaterialArrays::synchronizeTangentOperatorBlocksOnHostForReading()
      const {
    if (this->tangent_operator_blocks.Size() != 0) {
      static_cast<void>(this->tangent_operator_blocks.HostRead());
    }
  }  // end of synchronizeTangentOperatorBlocksOnHostForReading

  const mfem::Vector& MaterialArrays::getArray(const StateArray a,
                                               const TimeStepStage s) const {
    return this->state_arrays[getArrayIndex(a)][getStateIndex(s)];
  }  // end of getArray

  mfem::Vector& MaterialArrays::getMutableArray(const StateArray a,
                                                const TimeStepStage s) {
    return this->state_arrays[getArrayIndex(a)][getStateIndex(s)];
  }  // end of getMutableArray

  const mfem::Vector& MaterialArrays::getTangentOperatorBlocks() const {
    return this->tangent_operator_blocks;
  }  // end of getTangentOperatorBlocks

  real* MaterialArrays::readWrite(const StateArray a,
                                  const TimeStepStage s,
                                  const bool d) {
    return this->getMutableArray(a, s).ReadWrite(d);
  }  // end of readWrite

  real* MaterialArrays::write(const StateArray a,
                              const TimeStepStage s,
                              const bool d) {
    return this->getMutableArray(a, s).Write(d);
  }  // end of write

  real* MaterialArrays::readWriteTangentOperatorBlocks(const bool d) {
    return this->tangent_operator_blocks.ReadWrite(d);
  }  // end of readWriteTangentOperatorBlocks

  real* MaterialArrays::writeTangentOperatorBlocks(const bool d) {
    return this->tangent_operator_blocks.Write(d);
  }  // end of writeTangentOperatorBlocks

  MaterialArrays::~MaterialArrays() = default;

  mgis::MemoryLocation getDeviceMemoryLocation() noexcept {
    return mfem::Device::Allows(mfem::Backend::DEBUG_DEVICE)
               ? mgis::MemoryLocation::HOST
               : mgis::MemoryLocation::DEVICE;
  }  // end of getDeviceMemoryLocation

}  // end of namespace mfem_mgis
