/*!
 * \file   include/MFEMMGIS/ParaviewExportResults.hxx
 * \brief  This file declares the `ParaviewExportResults` class
 * \author Thomas Helfer
 * \date   24/03/2021
 */

#ifndef LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX
#define LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX

#include <string>
#include <utility>
#include <vector>
#include <functional>
#include <string_view>
#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/TimeStepStage.hxx"
#include "MFEMMGIS/PostProcessing/NonLinearEvolutionProblemPostProcessingBase.hxx"

namespace mfem_mgis {

  /*!
   * \brief a post-processing to export the results to paraview: the
   * unknowns and, optionally, other fields given by the calling code (see
   * `addField`), in the same files
   */
  template <bool parallel>
  struct ParaviewExportResults final
      : public NonLinearEvolutionProblemPostProcessingBase<parallel> {
    /*!
     * \brief function giving a field exported with the unknowns (see
     * `addField`). It is called at each export, before the files are
     * written, with the time step stage of the exported unknowns: the
     * beginning of the time step for the export at the initial time, the
     * end of the time step otherwise. It returns the field, defined on the
     * mesh of the export, which must stay alive until the files are
     * written, or a null pointer on failure, the error being registered in
     * the execution context.
     */
    using FieldFunction =
        std::function<GridFunction<parallel>*(Context&, const TimeStepStage)>;
    /*!
     * \brief constructor
     * \param[in, out] ctx: execution context
     * \param[in] pb: non linear problem
     * \param[in] params: parameters passed to the post-processing
     */
    ParaviewExportResults(mgis::Context& ctx,
                          NonLinearEvolutionProblemImplementation<parallel>& pb,
                          const Parameters& params);
    /*!
     * \brief adds a field exported with the unknowns, in the same files:
     * the mesh is written once per export for all the fields
     * \param[in] n: name of the field in the files
     * \param[in] f: function giving the field at each export
     */
    void addField(std::string_view n, FieldFunction f);
    /*!
     * \brief execute the post-processing at the initial time
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem
     * \param[in] t: initial time
     * \return true on success
     */
    [[nodiscard]] bool executeInitialPostProcessing(
        mgis::Context& ctx,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const real t) noexcept override;
    /*!
     * \brief execute the post-processing
     * \param[in, out] ctx: execution context
     * \param[in] p: non linear evolution problem, unused
     * \param[in] t: time at the beginning of the time step
     * \param[in] dt: time increment
     * \return true on success
     */
    [[nodiscard]] bool execute(
        mgis::Context& ctx,
        NonLinearEvolutionProblemImplementation<parallel>& p,
        const real t,
        const real dt) noexcept override;
    //! \brief destructor
    ~ParaviewExportResults() override;

   private:
    /*!
     * \brief writes the files of an export
     * \return true on success
     * \param[in, out] ctx: execution context
     * \param[in] u: exported unknowns (true degrees of freedom), copied in
     * `result`. A reference to the unknowns (`mfem::GridFunction::MakeTRef`)
     * would be an alias whose validity flags are the ones of the unknowns
     * when it is made: with or without prolongation, MFEM would copy stale
     * values of the host over the unknowns computed on a device.
     * \param[in] t: time
     * \param[in] s: time step stage of the unknowns
     */
    [[nodiscard]] bool exportResults(Context& ctx,
                                     const mfem::Vector& u,
                                     const real t,
                                     const TimeStepStage s) noexcept;
    //! \brief paraview exporter
    mfem::ParaViewDataCollection exporter;
    //! \brief exported grid function
    mfem_mgis::GridFunction<parallel> result;
    //! \brief other exported fields: names and functions giving them
    std::vector<std::pair<std::string, FieldFunction>> fields;
    //! \brief number of records
    size_type cycle;
    /*!
     * \brief boolean stating if the results shall be exported at the initial
     * time of the simulation
     */
    const bool shallExecuteInitialPostProcessing;
    //! \brief submesh defined when exporting data for domain or boundary
    //! attributes
    std::shared_ptr<mfem_mgis::SubMesh<parallel>> submesh;
    //! \brief fespace defined when exporting data for domain or boundary
    //! attributes
    std::shared_ptr<mfem_mgis::FiniteElementSpace<parallel>> fes_sm;
    //! \brief result_sm is used to transfer data from the result gridfunction
    //! on submesh
    std::shared_ptr<mfem_mgis::GridFunction<parallel>> result_sm;
  };  // end of struct ParaviewExportResults

}  // end of namespace mfem_mgis

#include "MFEMMGIS/ParaviewExportResults.ixx"

#endif /* LIB_MFEMMGIS_PARAVIEWEXPORTRESULTS_HXX */
