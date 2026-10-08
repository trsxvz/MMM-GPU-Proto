/*!
 * \file MFront/GenericParallel/MaterialProperty/BackendBase.hxx
 *
 * \brief  This file declares the `BackendBase` class
 * \author Thomas Helfer
 * \date   07/04/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_BACKENDBASE_HXX
#define LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_BACKENDBASE_HXX

#include <string>
#include <vector>
#include <functional>
#include <string_view>
#include "MFront/MFrontConfig.hxx"
#include "MFront/GenericParallel/MaterialProperty/AbstractBackend.hxx"

namespace mfront {

  // forward declaration
  struct VariableDescription;

}  // end of namespace mfront

namespace mfront::generic_parallel::material_property {

  //! \brief this class provides useful methods for most backends
  struct MFRONT_VISIBILITY_EXPORT BackendBase : AbstractBackend {
    //! \brief default constructor
    BackendBase();
    //
    void writeCDeclarations(std::ostream&,
                            const GenericParallelMaterialPropertyInterface&,
                            const MaterialPropertyDescription&) const override;
    [[nodiscard]] std::size_t getNumberOfStatuses(
        const MaterialPropertyDescription&) const override;
    void writeStatusesReportFunction(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&) const override;
    // do nothing implementation
    void writeProxyObjects(std::ostream&,
                           const GenericParallelMaterialPropertyInterface&,
                           const MaterialPropertyDescription&) const override;
    //! \brief destructor
    ~BackendBase() noexcept override;

   protected:
    /*!
     * \return true if the kernels report the violations of the bounds and
     * the invalid results in an array of statuses (runtime checks enabled)
     * \param[in] mpd: material property description
     */
    [[nodiscard]] static bool isBoundsCheckingRequired(
        const MaterialPropertyDescription&);
    /*!
     * \brief write the beginning of a C function computing the material
     * property on a set of points: signature, error reporting
     * (`mfront_report`), initialization of the output status, checks of the
     * arguments and of the strides. The function returns if there is no
     * point. If strides are treated, the strides of the arguments are copied
     * in the `mfront_strides` array, and `mfront_areAllArgumentsUniform`
     * states if all the arguments are uniform. The execution options are
     * copied in `mfront_options` (default options if the pointer given by
     * the calling code is null), and `mfront_asynchronous` states if the
     * statuses are given by the calling code (asynchronous evaluation).
     *
     * \param[in] os: output stream
     * \param[in] i: interface
     * \param[in] mpd: material property description
     * \param[in] treatStrides: if true, the function takes the strides of
     * the output and of the arguments. Otherwise, its name ends with `2`.
     */
    void writeCFunctionBeginning(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&,
        const bool) const;
    /*!
     * \brief write the body of a kernel, which evaluates the material
     * property at the point `mfront_idx`. The kernel reads the arguments in
     * `mfront_arg0`, `mfront_arg1`, etc. (with the strides `mfront_strides`
     * if strides are handled) and writes the result in `mfront_output`
     * (with the stride `mfront_output_stride` if strides are handled). The
     * parameters are read in `mfront_parameters`, unless they are treated
     * as static variables, and the violations of the bounds and the invalid
     * results are written in `mfront_bounds_statuses`.
     *
     * \param[in] os: output stream
     * \param[in] i: interface
     * \param[in] mpd: material property description
     * \param[in] handleStrides: if true, the arguments and the output are
     * strided
     */
    void writeKernelBody(std::ostream&,
                         const GenericParallelMaterialPropertyInterface&,
                         const MaterialPropertyDescription&,
                         const bool) const;
    /*!
     * \brief write the selection of the kernel called at each point.
     *
     * If strides are treated, a single point is evaluated if all the
     * arguments are uniform and the output stride is null, the kernel
     * without strides is used if all the strides are equal to one, and the
     * kernel handling strides otherwise.
     *
     * \param[in] os: output stream
     * \param[in] i: interface
     * \param[in] mpd: material property description
     * \param[in] treatStrides: boolean stating if strides are treated
     * \param[in] write_call: function writing the evaluation of the given
     * number of points (`1` or `mfront_npoints`) by the kernel handling
     * strides if its second argument is true, or by the kernel without
     * strides otherwise
     */
    void writeKernelSelection(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&,
        const bool,
        const std::function<void(std::string_view, bool)>&) const;
    /*!
     * \brief write the analysis of the results, once the kernels are
     * completed (see `writeStatusesReportFunction`), if the runtime checks
     * are enabled and if the evaluation is synchronous
     *
     * \param[in] os: output stream
     * \param[in] i: interface
     * \param[in] mpd: material property description
     */
    void writeResultsAnalysis(std::ostream&,
                              const GenericParallelMaterialPropertyInterface&,
                              const MaterialPropertyDescription&) const;
    //
    /*!
     * \brief write the check of the physical bounds of a variable in a kernel
     * \param[in] os: output stream
     * \param[in] v: variable
     * \param[in] i: index of the variable (the output is after the inputs)
     */
    virtual void writePhysicalBounds(std::ostream&,
                                     const VariableDescription&,
                                     const size_t) const;
    /*!
     * \brief write the check of the standard bounds of a variable in a kernel
     * \param[in] os: output stream
     * \param[in] prefix: prefix of the out of bounds policies
     * \param[in] v: variable
     * \param[in] i: index of the variable (the output is after the inputs)
     */
    virtual void writeBounds(std::ostream&,
                             const std::string&,
                             const VariableDescription&,
                             const size_t) const;
    /*!
     * \brief write the code analysing the array `mfront_bounds_statuses`,
     * in the function reporting the result (see
     * `writeStatusesReportFunction`)
     *
     * \param[in] os: output stream
     * \param[in] i: interface
     * \param[in] mpd: material property description
     */
    virtual void writeBoundsStatusesAnalysis(
        std::ostream&,
        const GenericParallelMaterialPropertyInterface&,
        const MaterialPropertyDescription&) const;
    /*!
     * \return the size of the array of statuses filled by the kernels. For
     * each argument and for the output, the array contains a status for the
     * physical bounds and a status for the standard bounds. The last status
     * states if the result is invalid (`nan` or infinite) at one point
     * at least.
     * \param[in] mpd: material property description
     */
    [[nodiscard]] static std::size_t getStatusesArraySize(
        const MaterialPropertyDescription&);
    /*!
     * \return the index of the status stating that the result is invalid at
     * one point at least.
     * \param[in] mpd: material property description
     */
    [[nodiscard]] static std::size_t getInvalidResultStatusIndex(
        const MaterialPropertyDescription&);
    /*!
     * \brief write the declaration of the array of statuses filled by the
     * kernels (`mfront_bounds_statuses`).
     *
     * For an asynchronous evaluation, this array is given by the execution
     * options and initialized by the calling code. Otherwise, it is a block
     * of a pool of blocks accessible by the host and by the device (see
     * `mfront::generic_parallel::StatusBlockPool`), initialized to zero: no
     * memory is allocated per call. The device only accesses the array if a
     * bound is violated or if the result is invalid.
     *
     * \param[in] os: output stream
     * \param[in] mpd: material property description
     * \param[in] device_type: type identifying a device
     * \param[in] device: expression identifying the current device
     * \param[in] allocation: body of a function allocating a block of
     * `mfront_n` integers, returning a null pointer on failure
     */
    static void writeBoundsStatusesDeclaration(
        std::ostream&,
        const MaterialPropertyDescription&,
        const std::string_view,
        const std::string_view,
        const std::string_view);
    /*!
     * \brief write the code setting the minimal stack size of the threads
     * of the device (see `MFront/GenericParallel/DeviceStackSize.hxx`), or
     * checking it if the execution options require to preserve it
     * \param[in] os: output stream
     * \param[in] stack_size: minimal stack size, in bytes. Nothing is
     * written if null
     */
    static void writeDeviceStackSizeSetting(std::ostream&, const std::size_t);
    /*!
     * \brief write the declaration of the `mfront_report` function, writing
     * an error message in the output status `mfront_output_status`
     * \param[in] os: output stream
     */
    static void writeErrorReporting(std::ostream&);
    /*!
     * \brief write the code marking the result as invalid if it is not
     * finite. This code is written in the kernel, after the evaluation of
     * the material property.
     * \param[in] os: output stream
     * \param[in] mpd: material property description
     */
    virtual void writeInvalidResultCheck(
        std::ostream&, const MaterialPropertyDescription&) const;
    /*!
     * \brief return if the backend shall handle the data transfer from a host
     * to a device
     *
     * If this method returns true, the argument of the C function
     * implementations containing pointers to values of the arguments
     * is named `mfront_args_host`. The backend shall
     * then define a pointer named `mfront_args` to the beginning of an
     * internal array of pointers. The output argument is then named
     * `mfront_output_host`, and the backend shall define a pointer named
     * `mfront_output`.
     *
     * If this method returns false, those arguments are directly named
     * `mfront_args` and `mfront_output`.
     *
     * \note by default, no data transfer is handled
     */
    virtual bool handlesDataTransfer() const;
  };

}  // end of namespace mfront::generic_parallel::material_property

#endif /* LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_BACKENDBASE_HXX */
