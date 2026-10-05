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
    // do nothing implementation
    void writeProxyObjects(std::ostream&,
                           const GenericParallelMaterialPropertyInterface&,
                           const MaterialPropertyDescription&) const override;
    //! \brief destructor
    ~BackendBase() noexcept override;

   protected:
    /*!
     * \brief structure decribing the extra arguments of the C functions
     */
    struct ExtraArgumentOfCFunctions {
      //! \brief type of the argument
      std::string type;
      //! \brief name of the argument. Shall start with `mfront_`
      std::string name;
      //! \brief description of the argument.
      std::string description;
      //! \brief boolean stating if the argument is passed as a pointer
      bool is_pointer;
      //! \brief boolean stating if the argument is mutable
      bool is_mutable;
    };
    //! \return the list of extra argument to the C functions
    [[nodiscard]] virtual std::vector<ExtraArgumentOfCFunctions>
    getExtraArgumentsOfCFunctions() const;
    //
    virtual void writePhysicalBounds(std::ostream&,
                                     const VariableDescription&,
                                     const size_t,
                                     const bool) const;
    virtual void writeBounds(std::ostream&,
                             const std::string&,
                             const VariableDescription&,
                             const size_t,
                             const bool) const;
    /*!
     * \brief write the code analysying the array `mfront_bounds_statuses`
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
     * If this method return true, the argument of the C function
     * implementations containing pointers to values of the state variables
     * is named `mfront_args_host`. The backend shall
     * then defined a pointer named `mfront_args` to the beginning of an
     * internal array of pointers.
     *
     * If this method returns false, this argument is directly named
     * `mfront_args`.
     *
     * \note by default, no data transfer is handled
     */
    virtual bool handlesDataTransfer() const;
  };

}  // end of namespace mfront::generic_parallel::material_property

#endif /* LIB_MFRONT_GENERICPARALLEL_MATERIALPROPERTY_BACKENDBASE_HXX */
