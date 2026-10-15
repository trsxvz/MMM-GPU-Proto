/*!
 * \file   include/MFEMMGIS/PreconditionerFactory.hxx
 * \brief  This file declares the `PreconditionerFactory` class
 * \date   15/10/2026
 */

#ifndef LIB_MFEM_MGIS_PRECONDITIONERFACTORY_HXX
#define LIB_MFEM_MGIS_PRECONDITIONERFACTORY_HXX

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/LinearSolverDescription.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \brief factory of the preconditioners of the linear solvers, selected by
   * the `Preconditioner` parameter of a linear solver (a name, `Name`, and
   * options, `Options`). Each preconditioner is registered with its
   * description (see `LinearSolverDescription`).
   * \tparam parallel: boolean stating if parallel preconditioners are
   * considered
   */
  template <bool parallel>
  struct MFEM_MGIS_EXPORT PreconditionerFactory {
    /*!
     * \brief generator of a preconditioner: the options are the ones given
     * by the `Options` parameter of the preconditioner (empty if not given)
     */
    using Generator = std::function<std::unique_ptr<LinearSolverPreconditioner>(
        Context&, FiniteElementSpace<parallel>&, const Parameters&)>;
    //! \return the unique instance of the class
    static PreconditionerFactory& getFactory();
    /*!
     * \brief register a new preconditioner
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the preconditioner
     * \param[in] d: description of the preconditioner
     * \param[in] g: generator of the preconditioner
     * \return true on success
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           LinearSolverDescription d,
                           Generator g) noexcept;
    /*!
     * \return the requested preconditioner, or a null pointer on error
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the preconditioner
     * \param[in] fespace: finite element space
     * \param[in] opts: options of the preconditioner
     */
    [[nodiscard]] std::unique_ptr<LinearSolverPreconditioner> generate(
        Context& ctx,
        std::string_view n,
        FiniteElementSpace<parallel>& fespace,
        const Parameters& opts) const;
    /*!
     * \return the description of a preconditioner, or a null pointer if no
     * preconditioner has this name
     * \param[in] n: name of the preconditioner
     */
    [[nodiscard]] const LinearSolverDescription* getDescription(
        std::string_view n) const noexcept;
    //! \return the names of the registered preconditioners
    [[nodiscard]] std::vector<std::string> getNames() const;

   private:
    //! \brief a registered preconditioner
    struct Entry {
      //! \brief description
      LinearSolverDescription description;
      //! \brief generator
      Generator generator;
    };
    //! \brief default constructor, which registers the preconditioners
    PreconditionerFactory();
    //! \brief destructor
    ~PreconditionerFactory();
    //! \brief registered preconditioners
    std::map<std::string, Entry, std::less<>> entries;
  };  // end of struct PreconditionerFactory

#ifdef MFEM_USE_MPI
  extern template struct PreconditionerFactory<true>;
#endif /* MFEM_USE_MPI */
  extern template struct PreconditionerFactory<false>;

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_PRECONDITIONERFACTORY_HXX */
