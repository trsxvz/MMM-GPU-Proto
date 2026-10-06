/*!
 * \file   mfront/include/MFront/GenericParallelBehaviour/SYCL.hxx
 * \brief  This file implements the integration of a batch of integration
 * points with SYCL.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELBEHAVIOUR_SYCL_HXX
#define LIB_MFRONT_GENERICPARALLELBEHAVIOUR_SYCL_HXX

#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include <sycl/sycl.hpp>
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::sycl {

  //! \brief policy for the transfer of the fields of a batch
  enum struct DataTransferPolicy {
    /*!
     * \brief the fields which are not accessible by the device associated
     * with the queue (host memory not allocated by SYCL) are copied in
     * device memory, and the outputs are copied back
     */
    AUTOMATIC,
    //! \brief all the fields are assumed to be accessible by the device
    NONE
  };

  namespace internals {

    /*!
     * \brief copies in device memory of the fields which are not accessible
     * by the device associated with a queue. The outputs must be copied back
     * by calling the `copyOutputs` method.
     */
    struct DeviceCopies {
      /*!
       * \brief constructor
       * \param[in] q: queue
       * \param[in] n: number of points of the fields
       */
      DeviceCopies(::sycl::queue& q, const size_type n)
          : number_of_points(n), queue(q) {}
      DeviceCopies(DeviceCopies&&) = delete;
      DeviceCopies(const DeviceCopies&) = delete;
      DeviceCopies& operator=(DeviceCopies&&) = delete;
      DeviceCopies& operator=(const DeviceCopies&) = delete;
      /*!
       * \brief copy a field in device memory if it is not accessible by the
       * device, and make the field refer to the copy. All the values covered
       * by the strides of the field are copied, so that the values of an
       * output which are not written are preserved.
       * \param[in,out] f: field
       * \param[in] size: number of components
       * \param[in] output: if true, the field is an output
       */
      template <typename FieldType>
      void transfer(FieldType& f, const size_type size, const bool output) {
        const auto n = this->number_of_points;
        if ((f.values == nullptr) || (size == 0) || (n == 0) ||
            (this->isAccessible(f.values))) {
          return;
        }
        const auto extent =
            (n - 1) * f.point_stride + (size - 1) * f.component_stride + 1;
        this->transfer(f.values, extent, output);
      }
      /*!
       * \brief copy an array in device memory if it is not accessible by the
       * device, and make the pointer refer to the copy.
       * \param[in,out] p: pointer to the array
       * \param[in] size: number of values
       * \param[in] output: if true, the array is an output
       */
      template <typename ValueType>
      void transfer(ValueType*& p, const size_type size, const bool output) {
        if ((p == nullptr) || (size == 0) || (this->isAccessible(p))) {
          return;
        }
        using value_type = std::remove_const_t<ValueType>;
        p = static_cast<value_type*>(
            this->copy(p, size * sizeof(value_type), output));
      }
      //! \brief copy the outputs back to their original locations
      void copyOutputs() {
        for (const auto& c : this->copies) {
          if (c.output) {
            this->queue.memcpy(c.host, c.device, c.size);
          }
        }
        this->queue.wait_and_throw();
      }
      //! \brief destructor
      ~DeviceCopies() {
        for (const auto& c : this->copies) {
          ::sycl::free(c.device, this->queue);
        }
      }

     private:
      //! \brief description of a copy in device memory
      struct Copy {
        //! \brief original location
        void* host;
        //! \brief copy
        void* device;
        //! \brief size in bytes
        std::size_t size;
        //! \brief if true, the copy must be copied back
        bool output;
      };
      //! \return if the given pointer is accessible by the device
      bool isAccessible(const void* const p) const {
        if (this->queue.get_device().is_cpu()) {
          return true;
        }
        return ::sycl::get_pointer_type(p, this->queue.get_context()) !=
               ::sycl::usm::alloc::unknown;
      }
      /*!
       * \brief copy an array in device memory
       * \return the copy
       * \param[in] p: array
       * \param[in] size: size in bytes
       * \param[in] output: if true, the copy must be copied back
       */
      void* copy(const void* const p,
                 const std::size_t size,
                 const bool output) {
        auto* const c = ::sycl::malloc_device(size, this->queue);
        if (c == nullptr) {
          throw std::runtime_error(
              "memory allocation failed on the SYCL device");
        }
        this->copies.push_back(Copy{const_cast<void*>(p), c, size, output});
        this->queue.memcpy(c, p, size);
        return c;
      }
      //! \brief number of points of the fields
      const size_type number_of_points;
      //! \brief copies
      std::vector<Copy> copies;
      //! \brief queue
      ::sycl::queue& queue;
    };

    /*!
     * \brief batch data whose fields are accessible by the device
     * associated with a queue. The fields which are not accessible are
     * copied in device memory. The outputs must be copied back by calling
     * the `copyOutputs` method.
     * \tparam PointTraits: description of the data of an integration point
     */
    template <typename PointTraits>
    struct DeviceBatchData : DeviceCopies {
      /*!
       * \brief constructor
       * \param[in] q: queue
       * \param[in] d: batch data
       */
      DeviceBatchData(::sycl::queue& q, const BatchData& d)
          : DeviceCopies(q, d.number_of_points), data(d) {
        constexpr auto sizes = PointDataSizes{PointTraits::value};
        auto& s0 = this->data.s0;
        auto& s1 = this->data.s1;
        this->transfer(s0.gradients, sizes.gradients, false);
        this->transfer(s1.gradients, sizes.gradients, false);
        this->transfer(s0.thermodynamic_forces, sizes.thermodynamic_forces,
                       false);
        this->transfer(s1.thermodynamic_forces, sizes.thermodynamic_forces,
                       true);
        this->transfer(s0.mass_density, 1, false);
        this->transfer(s1.mass_density, 1, false);
        this->transfer(s0.material_properties, sizes.material_properties,
                       false);
        this->transfer(s1.material_properties, sizes.material_properties,
                       false);
        this->transfer(s0.internal_state_variables,
                       sizes.internal_state_variables, false);
        this->transfer(s1.internal_state_variables,
                       sizes.internal_state_variables, true);
        this->transfer(s0.stored_energy, 1, false);
        this->transfer(s1.stored_energy, 1, true);
        this->transfer(s0.dissipated_energy, 1, false);
        this->transfer(s1.dissipated_energy, 1, true);
        this->transfer(s0.external_state_variables,
                       sizes.external_state_variables, false);
        this->transfer(s1.external_state_variables,
                       sizes.external_state_variables, false);
        constexpr const auto& mps = PointTraits::material_properties_sizes;
        constexpr const auto& esvs =
            PointTraits::external_state_variables_sizes;
        this->transferVariables(s0.material_properties_by_variable, mps,
                                this->material_properties0);
        this->transferVariables(s1.material_properties_by_variable, mps,
                                this->material_properties1);
        this->transferVariables(s0.external_state_variables_by_variable, esvs,
                                this->external_state_variables0);
        this->transferVariables(s1.external_state_variables_by_variable, esvs,
                                this->external_state_variables1);
        this->transfer(this->data.tangent_operators, sizes.tangent_operator,
                       true);
        this->transfer(this->data.rdts, 1, true);
        this->transfer(this->data.speed_of_sound, 1, true);
        this->transfer(this->data.statuses, this->data.number_of_points, true);
        q.wait_and_throw();
      }
      //! \brief batch data whose fields are accessible by the device
      BatchData data;

     private:
      /*!
       * \brief copy in device memory the fields of the variables of a group
       * given variable by variable which are not accessible by the device
       * \param[in,out] fields: fields of the variables, replaced by the
       * fields referring to the copies (this array is only read on the host)
       * \param[in] sizes: number of components of each variable
       * \param[out] copies: storage of the fields referring to the copies
       */
      template <std::size_t N>
      void transferVariables(const ConstField*& fields,
                             const std::array<size_type, N>& sizes,
                             std::vector<ConstField>& copies) {
        if (fields == nullptr) {
          return;
        }
        copies.assign(fields, fields + N);
        for (std::size_t v = 0; v != N; ++v) {
          this->transfer(copies[v], sizes[v], false);
        }
        fields = copies.data();
      }
      //! \brief fields of the material properties, variable by variable
      std::vector<ConstField> material_properties0, material_properties1;
      //! \brief fields of the external state variables, variable by variable
      std::vector<ConstField> external_state_variables0,
          external_state_variables1;
    };

    //! \brief a simple alias
    using AtomicCounter =
        ::sycl::atomic_ref<unsigned long long,
                           ::sycl::memory_order::relaxed,
                           ::sycl::memory_scope::device,
                           ::sycl::access::address_space::global_space>;

    /*!
     * \brief kernel applying a function at one point and updating the
     * statistics of the batch.
     *
     * The statistics are updated by atomic operations, which are only
     * performed if the point failed, gave an unreliable result or decreased
     * the time step scaling factor. Parallel reductions are not used since
     * their synchronisation barriers make the compilation of the large
     * kernels of some behaviours fail with AdaptiveCpp on CPUs.
     *
     * \tparam PointFunction: function object applied at one point
     */
    template <typename PointFunction>
    struct Kernel {
      void operator()(const ::sycl::id<1> i) const {
        const auto r = this->f(i[0]);
        if (r.number_of_failures != 0) {
          AtomicCounter{this->statistics->number_of_failures}.fetch_add(1);
        }
        if (r.number_of_unreliable_results != 0) {
          AtomicCounter{this->statistics->number_of_unreliable_results}
              .fetch_add(1);
        }
        const auto rdt = toOrderedKey(r.rdt);
        if (rdt < this->initial_rdt) {
          AtomicCounter{this->statistics->rdt}.fetch_min(rdt);
        }
      }
      //! \brief function object applied at one point
      PointFunction f;
      //! \brief statistics, allocated on the device
      AtomicBatchStatistics* statistics;
      //! \brief ordered key of the initial time step scaling factor
      unsigned long long initial_rdt;
    };

  }  // end of namespace internals

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch. This function waits for the completion of the
   * computations.
   *
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] q: queue
   * \param[in] rdt0: initial value of the time step scaling factor
   * \param[in] f: function
   * \param[in] device_stack_size: minimal stack size, in bytes, of the
   * threads of the device when the computations are performed on an NVIDIA
   * GPU through the CUDA runtime (see
   * `mfront::generic_parallel::setMinimalCUDAStackSize`). The stack size is
   * not modified if this value is null.
   * \return the status
   */
  template <typename PointFunction>
  int apply(OutputStatus& s,
            const size_type n,
            ::sycl::queue& q,
            const real rdt0,
            const PointFunction& f,
            const std::size_t device_stack_size = 0) {
    if (n == 0) {
      return s.status;
    }
    if (!mfront::generic_parallel::setMinimalCUDAStackSize(device_stack_size)) {
      return report(s, -3, "unable to set the stack size of the threads");
    }
    auto statistics = AtomicBatchStatistics::initial(rdt0);
    const auto deleter = [&q](AtomicBatchStatistics* const ptr) {
      ::sycl::free(ptr, q);
    };
    const auto device_statistics =
        std::unique_ptr<AtomicBatchStatistics, decltype(deleter)>(
            ::sycl::malloc_device<AtomicBatchStatistics>(1, q), deleter);
    if (device_statistics == nullptr) {
      return report(s, -3, "memory allocation failed on the SYCL device");
    }
    q.copy(&statistics, device_statistics.get(), 1).wait_and_throw();
    q.parallel_for(::sycl::range<1>{n},
                   internals::Kernel<PointFunction>{f, device_statistics.get(),
                                                    statistics.rdt})
        .wait_and_throw();
    q.copy(device_statistics.get(), &statistics, 1).wait_and_throw();
    return report(s, statistics.get());
  }  // end of apply

  /*!
   * \brief integrate the behaviour at the points of a batch.
   *
   * The `queue` member of the batch data must point to a `sycl::queue`. With
   * the `NONE` data transfer policy, the fields must be accessible by the
   * device associated with this queue (shared, device or host USM
   * allocations, or any memory if the device is the host). This function
   * waits for the completion of the computations.
   *
   * \tparam Behaviour: class implementing the behaviour
   * \tparam PointTraits: description of the data of an integration point
   * \tparam ParametersType: type of the structure holding the parameters
   * \param[out] s: output status
   * \param[in] d: batch data
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] transfers: policy for the transfer of the fields
   * \param[in] device_stack_size: minimal stack size, in bytes, of the
   * threads of the device when the computations are performed on an NVIDIA
   * GPU through the CUDA runtime (see
   * `mfront::generic_parallel::setMinimalCUDAStackSize`). The stack size is
   * not modified if this value is null.
   * \return the status
   */
  template <typename Behaviour, typename PointTraits, typename ParametersType>
  int integrate(
      OutputStatus& s,
      const BatchData& d,
      const ParametersType& parameters,
      const tfel::material::OutOfBoundsPolicy p,
      const DataTransferPolicy transfers = DataTransferPolicy::AUTOMATIC,
      const std::size_t device_stack_size = 0) {
    if (!check<PointTraits>(s, d)) {
      return s.status;
    }
    if (d.queue == nullptr) {
      return report(s, -2, "no SYCL queue given");
    }
    if (d.number_of_points == 0) {
      return s.status;
    }
    auto& q = *static_cast<::sycl::queue*>(d.queue);
    auto transferred =
        std::unique_ptr<internals::DeviceBatchData<PointTraits>>{};
    if (transfers == DataTransferPolicy::AUTOMATIC) {
      transferred =
          std::make_unique<internals::DeviceBatchData<PointTraits>>(q, d);
    }
    const auto& dd = transferred == nullptr ? d : transferred->data;
    const auto r = dispatch<Behaviour, PointTraits>(
        dd, parameters, p, [&](const auto& integrator) {
          return apply(s, dd.number_of_points, q, dd.rdt, integrator,
                       device_stack_size);
        });
    if (transferred != nullptr) {
      transferred->copyOutputs();
    }
    return r;
  }  // end of integrate

  /*!
   * \brief execute an initialize function at the points of a batch.
   *
   * \tparam Behaviour: class implementing the behaviour
   * \tparam PointTraits: description of the data of an integration point
   * \tparam InitializeFunction: function object calling the method
   * implementing the initialize function
   * \tparam number_of_inputs: number of values of the inputs
   * \param[out] s: output status
   * \param[in] d: batch data
   * \param[in] inputs: values of the inputs of the initialize function
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] transfers: policy for the transfer of the fields
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename InitializeFunction,
            size_type number_of_inputs,
            typename ParametersType>
  int executeInitializeFunction(
      OutputStatus& s,
      const BatchData& d,
      const ConstField& inputs,
      const ParametersType& parameters,
      const tfel::material::OutOfBoundsPolicy p,
      const DataTransferPolicy transfers = DataTransferPolicy::AUTOMATIC,
      const std::size_t device_stack_size = 0) {
    if ((!check<PointTraits>(s, d)) ||
        (!checkValues(s, inputs, number_of_inputs, "invalid inputs"))) {
      return s.status;
    }
    if (d.queue == nullptr) {
      return report(s, -2, "no SYCL queue given");
    }
    auto& q = *static_cast<::sycl::queue*>(d.queue);
    auto transferred =
        std::unique_ptr<internals::DeviceBatchData<PointTraits>>{};
    auto dinputs = inputs;
    if (transfers == DataTransferPolicy::AUTOMATIC) {
      transferred =
          std::make_unique<internals::DeviceBatchData<PointTraits>>(q, d);
      transferred->transfer(dinputs, number_of_inputs, false);
      q.wait_and_throw();
    }
    const auto& dd = transferred == nullptr ? d : transferred->data;
    const auto r =
        dispatchInitializeFunction<Behaviour, PointTraits, InitializeFunction,
                                   number_of_inputs>(
            dd, dinputs, parameters, p, [&](const auto& f) {
              return apply(s, dd.number_of_points, q, dd.rdt, f,
                           device_stack_size);
            });
    if (transferred != nullptr) {
      transferred->copyOutputs();
    }
    return r;
  }  // end of executeInitializeFunction

  /*!
   * \brief execute a post-processing at the points of a batch.
   *
   * \tparam Behaviour: class implementing the behaviour
   * \tparam PointTraits: description of the data of an integration point
   * \tparam PostProcessing: function object calling the method
   * implementing the post-processing
   * \tparam use_initial_state: boolean stating if the post-processing uses
   * the initial state of the material
   * \tparam number_of_outputs: number of values of the outputs
   * \param[out] s: output status
   * \param[out] outputs: values of the outputs of the post-processing
   * \param[in] d: batch data
   * \param[in] parameters: values of the parameters
   * \param[in] p: out of bounds policy
   * \param[in] transfers: policy for the transfer of the fields
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename Behaviour,
            typename PointTraits,
            typename PostProcessing,
            bool use_initial_state,
            size_type number_of_outputs,
            typename ParametersType>
  int executePostProcessing(
      OutputStatus& s,
      const Field& outputs,
      const BatchData& d,
      const ParametersType& parameters,
      const tfel::material::OutOfBoundsPolicy p,
      const DataTransferPolicy transfers = DataTransferPolicy::AUTOMATIC,
      const std::size_t device_stack_size = 0) {
    if ((!check<PointTraits>(s, d)) ||
        (!checkValues(s, outputs, number_of_outputs, "invalid outputs"))) {
      return s.status;
    }
    if (d.queue == nullptr) {
      return report(s, -2, "no SYCL queue given");
    }
    auto& q = *static_cast<::sycl::queue*>(d.queue);
    auto transferred =
        std::unique_ptr<internals::DeviceBatchData<PointTraits>>{};
    auto doutputs = outputs;
    if (transfers == DataTransferPolicy::AUTOMATIC) {
      transferred =
          std::make_unique<internals::DeviceBatchData<PointTraits>>(q, d);
      transferred->transfer(doutputs, number_of_outputs, true);
      q.wait_and_throw();
    }
    const auto& dd = transferred == nullptr ? d : transferred->data;
    const auto r =
        dispatchPostProcessing<Behaviour, PointTraits, PostProcessing,
                               use_initial_state, number_of_outputs>(
            doutputs, dd, parameters, p, [&](const auto& f) {
              return apply(s, dd.number_of_points, q, dd.rdt, f,
                           device_stack_size);
            });
    if (transferred != nullptr) {
      transferred->copyOutputs();
    }
    return r;
  }  // end of executePostProcessing

  /*!
   * \brief rotate the values of the points of a batch of an orthotropic
   * behaviour.
   *
   * \tparam PointRotation: function object rotating the values of a point
   * (see `RotationExecutor`)
   * \param[out] s: output status
   * \param[in] d: rotation data
   * \param[in] transfers: policy for the transfer of the fields
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   * \return the status
   */
  template <typename PointRotation>
  int rotate(OutputStatus& s,
             const RotationData& d,
             const DataTransferPolicy transfers = DataTransferPolicy::AUTOMATIC,
             const std::size_t device_stack_size = 0) {
    if (!checkRotationData(s, d)) {
      return s.status;
    }
    if (d.queue == nullptr) {
      return report(s, -2, "no SYCL queue given");
    }
    if (d.number_of_points == 0) {
      return s.status;
    }
    constexpr auto n = static_cast<size_type>(PointRotation::size);
    auto& q = *static_cast<::sycl::queue*>(d.queue);
    auto copies = std::optional<internals::DeviceCopies>{};
    auto dd = d;
    if (transfers == DataTransferPolicy::AUTOMATIC) {
      copies.emplace(q, d.number_of_points);
      copies->transfer(dd.source, n, false);
      copies->transfer(dd.rotation_matrices, 9, false);
      copies->transfer(dd.destination, n, true);
      q.wait_and_throw();
    }
    const auto r =
        apply(s, dd.number_of_points, q, 1, RotationExecutor<PointRotation>{dd},
              device_stack_size);
    if (copies.has_value()) {
      copies->copyOutputs();
    }
    return r;
  }  // end of rotate

}  // end of namespace mfront::gpb::sycl

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_SYCL_HXX */
