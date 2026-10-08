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

#include <cstddef>
#include <cstring>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#include <sycl/sycl.hpp>
#include "MFront/GenericParallel/DeviceStackSize.hxx"
#include "MFront/GenericParallel/SYCLDeviceCopies.hxx"
#include "MFront/GenericParallel/StatusBlockPool.hxx"
#include "MFront/GenericParallelBehaviour/Integrate.hxx"

namespace mfront::gpb::sycl {

  //! \brief policy for the transfer of the fields of a batch
  using mfront::generic_parallel::sycl::DataTransferPolicy;

  namespace internals {

    //! \brief copies of the fields which are not accessible by the device
    using mfront::generic_parallel::sycl::DeviceCopies;

    /*!
     * \brief batch data whose fields are accessible by the device
     * associated with a queue. With the `AUTOMATIC` data transfer policy,
     * the fields which are not accessible are copied in device memory. The
     * outputs must be copied back by calling the `copyOutputs` method.
     * \tparam PointTraits: description of the data of an integration point
     */
    template <typename PointTraits>
    struct DeviceBatchData : DeviceCopies {
      /*!
       * \brief constructor
       * \param[in] q: queue
       * \param[in] d: batch data
       * \param[in] p: policy for the transfer of the fields
       */
      DeviceBatchData(::sycl::queue& q,
                      const BatchData& d,
                      const DataTransferPolicy p)
          : DeviceCopies(q, d.number_of_points, p), data(d) {
        // the sizes of the thermodynamic forces and of the tangent operator
        // of the finite strain behaviours depend on the computation options
        using Integrator = typename PointTraits::Integrator;
        const auto sizes = Integrator::template getPointDataSizes<PointTraits>(
            getComputationOptions(d));
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
        // the fields of a group given variable by variable are ignored
        if (s0.material_properties_by_variable == nullptr) {
          this->transfer(s0.material_properties, sizes.material_properties,
                         false);
        }
        if (s1.material_properties_by_variable == nullptr) {
          this->transfer(s1.material_properties, sizes.material_properties,
                         false);
        }
        this->transfer(s0.internal_state_variables,
                       sizes.internal_state_variables, false);
        this->transfer(s1.internal_state_variables,
                       sizes.internal_state_variables, true);
        this->transfer(s0.stored_energy, 1, false);
        this->transfer(s1.stored_energy, 1, true);
        this->transfer(s0.dissipated_energy, 1, false);
        this->transfer(s1.dissipated_energy, 1, true);
        if (s0.external_state_variables_by_variable == nullptr) {
          this->transfer(s0.external_state_variables,
                         sizes.external_state_variables, false);
        }
        if (s1.external_state_variables_by_variable == nullptr) {
          this->transfer(s1.external_state_variables,
                         sizes.external_state_variables, false);
        }
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
        this->transfer(this->data.execution.statistics, 1, true);
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
       * \param[out] copies: storage of the fields referring to the copies,
       * only filled if a field is copied
       */
      template <std::size_t N>
      void transferVariables(const ConstField*& fields,
                             const std::array<size_type, N>& sizes,
                             std::vector<ConstField>& copies) {
        if ((fields == nullptr) || (!this->isEnabled())) {
          return;
        }
        auto accessible = true;
        for (std::size_t v = 0; (v != N) && (accessible); ++v) {
          accessible = (fields[v].values == nullptr) || (sizes[v] == 0) ||
                       (this->isAccessible(fields[v].values));
        }
        if (accessible) {
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
          AtomicCounter{this->statistics->first_failure_index}.fetch_min(
              static_cast<unsigned long long>(i[0]));
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
      //! \brief statistics, accessible by the device
      AtomicBatchStatistics* statistics;
      //! \brief ordered key of the initial time step scaling factor
      unsigned long long initial_rdt;
    };

    /*!
     * \return the pool of blocks of shared memory holding the statistics
     * of the synchronous computations, shared by all the functions of the
     * library. The blocks are identified by the device and the context of
     * the queue.
     */
    inline auto& getStatisticsBlocks() noexcept {
      using Pool = mfront::generic_parallel::StatusBlockPool<
          AtomicBatchStatistics, std::pair<::sycl::device, ::sycl::context>>;
      static auto pool = Pool{};
      return pool;
    }  // end of getStatisticsBlocks

  }  // end of namespace internals

  /*!
   * \brief apply a function at the points of a batch and report the
   * statistics of the batch.
   *
   * If the execution options give statistics, this function returns once
   * the computations are submitted to the queue. Otherwise, it waits for
   * their completion.
   *
   * \tparam PointFunction: function object taking the index of a point and
   * returning the statistics associated with this point
   * \param[in,out] s: output status, which must have been initialized
   * \param[in] n: number of points
   * \param[in] o: execution options. The queue must point to a
   * `sycl::queue`.
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
            const ExecutionOptions& o,
            const real rdt0,
            const PointFunction& f,
            const std::size_t device_stack_size = 0) {
    if (n == 0) {
      return s.status;
    }
    if (o.queue == nullptr) {
      return report(s, -2, "no SYCL queue given");
    }
    if (!mfront::generic_parallel::setMinimalCUDAStackSize(
            device_stack_size, o.preserve_stack_size != 0)) {
      return report(s, -3, mfront::generic_parallel::stack_size_error_message);
    }
    auto& q = *static_cast<::sycl::queue*>(o.queue);
    // the kernel must be executed after the computations already submitted
    // to the queue, which may compute the inputs. This is guaranteed by an
    // in-order queue. With an out-of-order queue, the computations are
    // synchronous and the queue is synchronized before the kernel.
    if (!q.is_in_order()) {
      if (o.statistics != nullptr) {
        return report(s, -2,
                      "asynchronous computations require an in-order SYCL "
                      "queue");
      }
      q.wait_and_throw();
    }
    const auto range = ::sycl::range<1>{n};
    const auto rdt_key = toOrderedKey(rdt0);
    using Kernel = internals::Kernel<PointFunction>;
    if (o.statistics != nullptr) {
      q.parallel_for(range, Kernel{f, o.statistics, rdt_key});
      return s.status;
    }
    // synchronous computations: the statistics are held by a block of
    // shared memory, initialized and read by the host. The device only
    // accesses this block if an event occurs (see internals::Kernel): no
    // data is transferred otherwise.
    const auto block = internals::getStatisticsBlocks().acquire(
        {q.get_device(), q.get_context()},
        [&q] { return ::sycl::malloc_shared<AtomicBatchStatistics>(1, q); });
    auto* const statistics = block.get();
    if (statistics == nullptr) {
      return report(s, -3, "unable to allocate the statistics of the batch");
    }
    mfront_gpb_initializeBatchStatistics(statistics, rdt0);
    q.parallel_for(range, Kernel{f, statistics, rdt_key}).wait_and_throw();
    return report(s, getBatchStatistics(*statistics));
  }  // end of apply

  /*!
   * \return the primitive of the backend, as expected by
   * `mfront::gpb::integrateBatch` and the similar functions of
   * `MFront/GenericParallelBehaviour/Integrate.hxx`
   * \param[in,out] s: output status
   * \param[in] device_stack_size: minimal stack size of the threads of the
   * device (see `apply`)
   */
  inline auto makeApply(OutputStatus& s,
                        const std::size_t device_stack_size = 0) {
    return [&s, device_stack_size](const size_type n, const ExecutionOptions& o,
                                   const real rdt0, const auto& f) {
      return sycl::apply(s, n, o, rdt0, f, device_stack_size);
    };
  }  // end of makeApply

  /*!
   * \return a copy of a data in the memory of the device associated with a
   * queue
   * \param[in] v: data
   * \param[in] q: queue
   */
  template <typename T>
  const T* copyToDevice(const T& v, ::sycl::queue& q) {
    static_assert(std::is_trivially_copyable_v<T>,
                  "the shared data must be trivially copyable");
    auto* const p = ::sycl::malloc_device<T>(1, q);
    if (p == nullptr) {
      throw std::runtime_error("memory allocation failed on the SYCL device");
    }
    q.copy(&v, p, 1).wait_and_throw();
    return p;
  }  // end of copyToDevice

  /*!
   * \return a copy of a data shared by all the integration points (see the
   * `@SharedData` keyword) in the memory of the device associated with the
   * queue given by the execution options. A copy is made at the first call
   * for each device and context, and is never freed.
   * \return the copy, or a null pointer if no queue is given (the call
   * reports this invalid argument before any computation)
   * \tparam get_shared_data: function returning the data
   * \param[in] o: execution options
   */
  template <auto get_shared_data>
  const auto* getSharedDataCopy(const ExecutionOptions& o) {
    using T = std::remove_cvref_t<decltype(get_shared_data())>;
    using Device = std::pair<::sycl::device, ::sycl::context>;
    static auto copies = SharedDataCopies<T, Device>{};
    if (o.queue == nullptr) {
      return static_cast<const T*>(nullptr);
    }
    auto& q = *static_cast<::sycl::queue*>(o.queue);
    return copies.get(Device{q.get_device(), q.get_context()},
                      [&q] { return copyToDevice(get_shared_data(), q); });
  }  // end of getSharedDataCopy

  namespace internals {

    /*!
     * \brief call a function with batch data whose fields are accessible by
     * the device associated with the queue of the batch (see
     * `DeviceBatchData`), and copy the outputs back. If some fields were
     * copied, this function waits for the completion of the computations.
     * \param[out] s: output status
     * \param[in] d: batch data
     * \param[in] transfers: policy for the transfer of the fields
     * \param[in] call: function called with the batch data and the copies,
     * which allow to transfer additional fields
     * \return the status
     */
    template <typename PointTraits, typename Function>
    int callWithDeviceData(OutputStatus& s,
                           const BatchData& d,
                           const DataTransferPolicy transfers,
                           Function&& call) {
      if (d.execution.queue == nullptr) {
        initialize(s, d.rdt);
        return report(s, -2, "no SYCL queue given");
      }
      auto& q = *static_cast<::sycl::queue*>(d.execution.queue);
      auto dd = DeviceBatchData<PointTraits>(q, d, transfers);
      const auto r = call(dd.data, static_cast<DeviceCopies&>(dd));
      dd.copyOutputs();
      return r;
    }  // end of callWithDeviceData

  }  // end of namespace internals

  /*!
   * \brief integrate the behaviour at the points of a batch (see
   * `mfront::gpb::integrateBatch`).
   *
   * The queue given by the execution options must point to a `sycl::queue`.
   * With the `NONE` data transfer policy, the fields must be accessible by
   * the device associated with this queue (shared, device or host USM
   * allocations, or any memory if the device is the host). The computations
   * are asynchronous if the execution options give statistics and if no
   * field is copied.
   *
   * \param[in] transfers: policy for the transfer of the fields
   * \param[in] apply: primitive of the backend (see `makeApply`)
   */
  template <typename Behaviour,
            typename PointTraits,
            typename ParametersType,
            typename Apply>
  int integrateBatch(OutputStatus& s,
                     const BatchData& d,
                     const ParametersType& parameters,
                     const tfel::material::OutOfBoundsPolicy p,
                     const DataTransferPolicy transfers,
                     const Apply& apply) {
    return internals::callWithDeviceData<PointTraits>(
        s, d, transfers, [&](const BatchData& dd, internals::DeviceCopies&) {
          return mfront::gpb::integrateBatch<Behaviour, PointTraits>(
              s, dd, parameters, p, apply);
        });
  }  // end of integrateBatch

  /*!
   * \brief execute an initialize function at the points of a batch (see
   * `mfront::gpb::executeInitializeFunctionOnBatch` and `integrateBatch`)
   */
  template <typename Behaviour,
            typename PointTraits,
            typename InitializeFunction,
            size_type number_of_inputs,
            typename ParametersType,
            typename Apply>
  int executeInitializeFunctionOnBatch(
      OutputStatus& s,
      const BatchData& d,
      const ConstField& inputs,
      const ParametersType& parameters,
      const tfel::material::OutOfBoundsPolicy p,
      const DataTransferPolicy transfers,
      const Apply& apply) {
    return internals::callWithDeviceData<PointTraits>(
        s, d, transfers,
        [&](const BatchData& dd, internals::DeviceCopies& copies) {
          auto dinputs = inputs;
          copies.transfer(dinputs, number_of_inputs, false);
          return mfront::gpb::executeInitializeFunctionOnBatch<
              Behaviour, PointTraits, InitializeFunction, number_of_inputs>(
              s, dd, dinputs, parameters, p, apply);
        });
  }  // end of executeInitializeFunctionOnBatch

  /*!
   * \brief execute a post-processing at the points of a batch (see
   * `mfront::gpb::executePostProcessingOnBatch` and `integrateBatch`)
   */
  template <typename Behaviour,
            typename PointTraits,
            typename PostProcessing,
            bool use_initial_state,
            size_type number_of_outputs,
            typename ParametersType,
            typename Apply>
  int executePostProcessingOnBatch(OutputStatus& s,
                                   const Field& outputs,
                                   const BatchData& d,
                                   const ParametersType& parameters,
                                   const tfel::material::OutOfBoundsPolicy p,
                                   const DataTransferPolicy transfers,
                                   const Apply& apply) {
    return internals::callWithDeviceData<PointTraits>(
        s, d, transfers,
        [&](const BatchData& dd, internals::DeviceCopies& copies) {
          auto doutputs = outputs;
          copies.transfer(doutputs, number_of_outputs, true);
          return mfront::gpb::executePostProcessingOnBatch<
              Behaviour, PointTraits, PostProcessing, use_initial_state,
              number_of_outputs>(s, doutputs, dd, parameters, p, apply);
        });
  }  // end of executePostProcessingOnBatch

  /*!
   * \brief rotate the values of the points of a batch of an orthotropic
   * behaviour (see `mfront::gpb::rotateBatch` and `integrateBatch`)
   */
  template <typename PointRotation, typename Apply>
  int rotateBatch(OutputStatus& s,
                  const RotationData& d,
                  const DataTransferPolicy transfers,
                  const Apply& apply) {
    if (d.execution.queue == nullptr) {
      initialize(s, 1);
      return report(s, -2, "no SYCL queue given");
    }
    constexpr auto n = static_cast<size_type>(PointRotation::size);
    auto copies =
        internals::DeviceCopies(*static_cast<::sycl::queue*>(d.execution.queue),
                                d.number_of_points, transfers);
    auto dd = d;
    copies.transfer(dd.source, n, false);
    copies.transfer(dd.rotation_matrices, 9, false);
    copies.transfer(dd.destination, n, true);
    copies.transfer(dd.execution.statistics, 1, true);
    const auto r = mfront::gpb::rotateBatch<PointRotation>(s, dd, apply);
    copies.copyOutputs();
    return r;
  }  // end of rotateBatch

}  // end of namespace mfront::gpb::sycl

#endif /* LIB_MFRONT_GENERICPARALLELBEHAVIOUR_SYCL_HXX */
