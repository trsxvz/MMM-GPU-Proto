/*!
 * \file   mfront/include/MFront/GenericModelInterface.hxx
 * \brief
 * \author Thomas Helfer
 * \date   28/04/2022
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICMODELINTERFACE_HXX
#define LIB_MFRONT_GENERICMODELINTERFACE_HXX

#include "MFront/AbstractModelInterface.hxx"

namespace mfront {

  /*!
   * \brief options of the generation of the class implementing a model
   */
  struct GenericModelClassOptions {
    /*!
     * \brief if true, the class can be used in device code: its constructor
     * takes the behaviour data and a structure holding the values of the
     * real parameters, its methods are host and device functions and, in
     * device code, the material laws are evaluated by their device
     * implementations (the functions of the `mfront` interface are used on
     * the host). String parameters are not supported in this case (an
     * exception is thrown).
     */
    bool device = false;
    //! \brief type of the structure holding the values of the parameters
    std::string parameters_type = "";
  };

  /*!
   * \brief write, in a namespace opened by the caller (`mfront::gm`), the
   * class implementing a model and, if the model has real parameters which
   * are not treated as static variables, the class holding their values
   * (see `getGenericModelParametersInitializerName`). Each function of the
   * model is implemented by a method returning false if a bounds check
   * fails, and the `execute` method executes all the functions of the model.
   * \param[out] os: output stream
   * \param[in] fd: file description
   * \param[in] md: model description
   * \param[in] o: options
   */
  MFRONT_VISIBILITY_EXPORT void writeGenericModelClass(
      std::ostream&,
      const FileDescription&,
      const ModelDescription&,
      const GenericModelClassOptions&);
  //! \return true if the given parameter is a real parameter
  MFRONT_VISIBILITY_EXPORT bool isGenericModelRealParameter(
      const VariableDescription&);
  /*!
   * \return true if the values of the real parameters of the given model are
   * held by the class whose name is returned by
   * `getGenericModelParametersInitializerName`
   */
  MFRONT_VISIBILITY_EXPORT bool hasGenericModelParametersInitializer(
      const ModelDescription&);
  //! \return the name of the class holding the values of the parameters
  MFRONT_VISIBILITY_EXPORT std::string getGenericModelParametersInitializerName(
      const ModelDescription&);

  /*!
   * \brief a generic interface for point-wise models implemented using the
   * `Model` domain specific language.
   *
   * This interface is based on the on-the-fly conversion of the
   * model description in a behaviour description and using the interface for
   * behaviours.
   */
  struct GenericModelInterface final : public AbstractModelInterface {
    //! \return the name of the interface
    [[nodiscard]] static std::string getName();
    //! \brief constructor
    GenericModelInterface();
    //
    void setOptions(const DataMap &) override;
    //! \return the name of the generated library
    [[nodiscard]] std::string getLibraryName(const ModelDescription &) const;
    [[nodiscard]] std::pair<bool, tokens_iterator> treatKeyword(
        const std::string &,
        const std::vector<std::string> &,
        tokens_iterator,
        const tokens_iterator) override;
    void declareReservedNames(std::set<std::string> &) override;
    void writeOutputFiles(const FileDescription &,
                          const ModelDescription &) override;
    void getTargetsDescription(TargetsDescription &,
                               const ModelDescription &) override;
    //! \brief destructor
    ~GenericModelInterface() override;

   protected:
    /*!
     * \brief generate the header file
     * \param fd: file description
     * \param md: model description
     */
    virtual void writeHeaderFile(const FileDescription &,
                                 const ModelDescription &);
    /*!
     * \brief generate the source file
     * \param fd: file description
     * \param md: model description
     */
    virtual void writeSourceFile(const FileDescription &,
                                 const ModelDescription &);

  };  // end of struct GenericModelInterface

}  // end of namespace mfront

#endif /* LIB_MFRONT_GENERICMODELINTERFACE_HXX */
