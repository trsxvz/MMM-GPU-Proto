/*!
 * \file   mfront/include/MFront/GenericParallel/Behaviour/AbstractBackend.hxx
 * \brief  This file declares the base class of the backends of the
 * `generic-parallel` behaviour interface
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_ABSTRACTBACKEND_HXX
#define LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_ABSTRACTBACKEND_HXX

#include <iosfwd>
#include <string>
#include <string_view>
#include "MFront/MFrontConfig.hxx"

namespace mfront::generic_parallel::behaviour {

  /*!
   * \brief base class of the backends of the `generic-parallel` behaviour
   * interface.
   *
   * The integration of a batch is implemented by function templates in the
   * `MFront/GenericParallelBehaviour` headers. A backend selects those
   * headers and writes the call to the right function.
   */
  struct MFRONT_VISIBILITY_EXPORT AbstractBackend {
    /*!
     * \return the name of the backend, which is used as a suffix for the
     * names of the generated library and files
     */
    [[nodiscard]] virtual std::string getName() const = 0;
    //! \return the extension of the generated source file
    [[nodiscard]] virtual std::string getSourceFileExtension() const = 0;
    /*!
     * \brief write the headers required by the backend
     * \param[out] os: output stream
     */
    virtual void writeSpecificIncludes(std::ostream&) const = 0;
    /*!
     * \brief write the statement returning the result of a function of the
     * backend, defined in the `MFront/GenericParallelBehaviour` headers
     * (`integrate`, `executeInitializeFunction` or `executePostProcessing`).
     * The template arguments and the arguments specific to the backend are
     * appended to the given ones.
     *
     * \param[out] os: output stream
     * \param[in] f: name of the function
     * \param[in] template_arguments: template arguments
     * \param[in] arguments: arguments
     */
    virtual void writeCall(std::ostream&,
                           std::string_view,
                           std::string_view,
                           std::string_view) const = 0;
    //! \brief destructor
    virtual ~AbstractBackend() noexcept;
  };  // end of AbstractBackend

}  // end of namespace mfront::generic_parallel::behaviour

#endif /* LIB_MFRONT_GENERICPARALLEL_BEHAVIOUR_ABSTRACTBACKEND_HXX */
