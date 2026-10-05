/*!
 * \file   mfront/include/MFront/GenericParallelModelInterface.hxx
 * \brief  This file declares the `generic-parallel` interface for models.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICPARALLELMODELINTERFACE_HXX
#define LIB_MFRONT_GENERICPARALLELMODELINTERFACE_HXX

#include <memory>
#include "MFront/AbstractModelInterface.hxx"
#include "MFront/GenericParallel/Behaviour/AbstractBackend.hxx"

namespace mfront {

  /*!
   * \brief the `generic-parallel` interface for models.
   *
   * As with the `generic` interface for models, models are treated as
   * behaviours without gradients: the outputs of the model are the internal
   * state variables and its inputs are the external state variables. The
   * generated functions execute the model over a batch of points, with the
   * backends and the batch data of the `generic-parallel` interface for
   * behaviours.
   */
  struct GenericParallelModelInterface final : public AbstractModelInterface {
    //! \return the name of the interface
    [[nodiscard]] static std::string getName();
    //! \brief constructor
    GenericParallelModelInterface();
    //
    void setOptions(const DataMap&) override;
    [[nodiscard]] std::pair<bool, tokens_iterator> treatKeyword(
        const std::string&,
        const std::vector<std::string>&,
        tokens_iterator,
        const tokens_iterator) override;
    void declareReservedNames(std::set<std::string>&) override;
    void writeOutputFiles(const FileDescription&,
                          const ModelDescription&) override;
    void getTargetsDescription(TargetsDescription&,
                               const ModelDescription&) override;
    //! \brief destructor
    ~GenericParallelModelInterface() override;

   private:
    //! \return the name of the generated library
    [[nodiscard]] std::string getLibraryName(const ModelDescription&) const;
    //! \return the name of the header file
    [[nodiscard]] std::string getHeaderFileName(const std::string&) const;
    //! \return the name of the source file
    [[nodiscard]] std::string getSourceFileName(const std::string&) const;
    //! \brief write the header file
    void writeHeaderFile(const FileDescription&, const ModelDescription&);
    //! \brief write the source file
    void writeSourceFile(const FileDescription&, const ModelDescription&);
    //! \brief backend
    std::unique_ptr<generic_parallel::behaviour::AbstractBackend> backend;
  };  // end of struct GenericParallelModelInterface

}  // end of namespace mfront

#endif /* LIB_MFRONT_GENERICPARALLELMODELINTERFACE_HXX */
