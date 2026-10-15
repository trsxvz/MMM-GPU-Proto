/*!
 * \file   include/MFEMMGIS/LinearSolverFactory.hxx
 * \brief  This file declares the `LinearSolverFactory` class
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEM_MGIS_LINEARSOLVERFACTORY_HXX
#define LIB_MFEM_MGIS_LINEARSOLVERFACTORY_HXX

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <string_view>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/LinearSolverHandler.hxx"
#include "MFEMMGIS/LinearSolverDescription.hxx"

namespace mfem_mgis {

  // forward declaration
  struct Parameters;

  /*!
   * \brief an abstract factory for linear solvers. Each linear solver is
   * registered with its description (see `LinearSolverDescription`). The
   * preconditioners are registered in the `PreconditionerFactory`.
   * \tparam parallel: boolean stating if parallel linear solvers are
   * considered
   *
   * \note a linear solver which is neither an iterative solver of MFEM nor
   * a Krylov solver of hypre reports the outcome of its resolutions by
   * implementing the `LinearSolverConvergenceReport` interface (see the
   * `hasConverged` function in `SolverUtilities.hxx`)
   */
  template <bool parallel>
  struct LinearSolverFactory;

#ifdef MFEM_USE_MPI

  //! \brief specialisation in parallel
  template <>
  struct MFEM_MGIS_EXPORT LinearSolverFactory<true> {
    //! \brief a simple alias
    using Generator = std::function<LinearSolverHandler(
        Context&, FiniteElementSpace<true>&, const Parameters&)>;
    //! \return the unique instance of the class
    static LinearSolverFactory& getFactory();
    /*!
     * \brief register a new linear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] d: description of the linear solver
     * \param[in] g: generator of the linear solver
     * \return true on success
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           LinearSolverDescription d,
                           Generator g) noexcept;
    /*!
     * \brief register a new linear solver, whose description states that
     * it computes where the device of MFEM computes and accepts any
     * preconditioner (the `Preconditioner` parameter being left to its
     * generator)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] g: generator of the linear solver
     * \return true on success
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           Generator g) noexcept;
    /*!
     * \return the description of a linear solver, or a null pointer if no
     * linear solver has this name
     * \param[in] n: name of the linear solver
     */
    [[nodiscard]] const LinearSolverDescription* getDescription(
        std::string_view n) const noexcept;
    //! \return the names of the registered linear solvers
    [[nodiscard]] std::vector<std::string> getNames() const;
    /*!
     * \return the requested linear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] fespace: finite element space
     * \param[in] params: parameters passed to the linear solver
     */
    LinearSolverHandler generate(Context& ctx,
                                 std::string_view n,
                                 FiniteElementSpace<true>& fespace,
                                 const Parameters& params) const;

   private:
    //! \brief default constructor
    LinearSolverFactory();
    //! \brief destructor
    ~LinearSolverFactory();
    //! \brief a registered linear solver
    struct Entry {
      //! \brief description
      LinearSolverDescription description;
      //! \brief generator
      Generator generator;
    };
    //! \brief registered linear solvers
    std::map<std::string, Entry, std::less<>> entries;
  };  // end of struct LinearSolverFactory

#endif /* MFEM_USE_MPI */

  //! \brief specialisation in sequential
  template <>
  struct MFEM_MGIS_EXPORT LinearSolverFactory<false> {
    //! \brief a simple alias
    using Generator = std::function<LinearSolverHandler(
        Context&, FiniteElementSpace<false>&, const Parameters&)>;
    //! \return the unique instance of the class
    static LinearSolverFactory& getFactory();
    /*!
     * \brief register a new linear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] d: description of the linear solver
     * \param[in] g: generator of the linear solver
     * \return true on success
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           LinearSolverDescription d,
                           Generator g) noexcept;
    /*!
     * \brief register a new linear solver, whose description states that
     * it computes where the device of MFEM computes and accepts any
     * preconditioner (the `Preconditioner` parameter being left to its
     * generator)
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] g: generator of the linear solver
     * \return true on success
     */
    [[nodiscard]] bool add(Context& ctx,
                           std::string_view n,
                           Generator g) noexcept;
    /*!
     * \return the description of a linear solver, or a null pointer if no
     * linear solver has this name
     * \param[in] n: name of the linear solver
     */
    [[nodiscard]] const LinearSolverDescription* getDescription(
        std::string_view n) const noexcept;
    //! \return the names of the registered linear solvers
    [[nodiscard]] std::vector<std::string> getNames() const;
    /*!
     * \return the requested linear solver
     * \param[in, out] ctx: execution context
     * \param[in] n: name of the linear solver
     * \param[in] fespace: finite element space
     * \param[in] params: parameters passed to the linear solver
     */
    LinearSolverHandler generate(Context& ctx,
                                 std::string_view n,
                                 FiniteElementSpace<false>& fespace,
                                 const Parameters& params) const;

   private:
    //! \brief default constructor
    LinearSolverFactory();
    //! \brief destructor
    ~LinearSolverFactory();
    //! \brief a registered linear solver
    struct Entry {
      //! \brief description
      LinearSolverDescription description;
      //! \brief generator
      Generator generator;
    };
    //! \brief registered linear solvers
    std::map<std::string, Entry, std::less<>> entries;
  };  // end of struct LinearSolverFactory

}  // end of namespace mfem_mgis

#endif /* LIB_MFEM_MGIS_LINEARSOLVERFACTORY_HXX */
