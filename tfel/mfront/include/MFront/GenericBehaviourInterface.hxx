/*!
 * \file   mfront/include/GenericBehaviourInterface.hxx
 * \brief
 * \author Thomas Helfer
 * \date 30/06/2018
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#ifndef LIB_MFRONT_GENERICBEHAVIOURINTERFACE_HXX
#define LIB_MFRONT_GENERICBEHAVIOURINTERFACE_HXX

#include <set>
#include <optional>
#include "MFront/BehaviourInterfaceBase.hxx"

namespace mfront {

  /*!
   * \brief a generic behaviour interface
   * See https://github.com/thelfer/MFrontGenericInterfaceSupport
   */
  struct MFRONT_VISIBILITY_EXPORT GenericBehaviourInterface
      : public BehaviourInterfaceBase {
    //! \return the interface name
    [[nodiscard]] static std::string getName();
    /*!
     * \return the current API version
     *
     * This API version changes everytime the API of the generic interface
     * changes.
     */
    [[nodiscard]] static unsigned short getAPIVersion();
    //
    [[nodiscard]] std::string getInterfaceName() const override;
    [[nodiscard]] std::pair<bool, tokens_iterator> treatKeyword(
        BehaviourDescription&,
        const std::string&,
        const std::vector<std::string>&,
        tokens_iterator,
        const tokens_iterator) override;
    [[nodiscard]] std::set<Hypothesis> getModellingHypothesesToBeTreated(
        const BehaviourDescription&) const override;
    void writeInterfaceSpecificIncludes(
        std::ostream&, const BehaviourDescription&) const override;
    void endTreatment(const BehaviourDescription&,
                      const FileDescription&) const override;
    void getTargetsDescription(TargetsDescription&,
                               const BehaviourDescription&) override;
    void writeBehaviourConstructorHeader(std::ostream&,
                                         const BehaviourDescription&,
                                         const Hypothesis,
                                         const std::string&) const override;
    void writeBehaviourConstructorBody(std::ostream&,
                                       const BehaviourDescription&,
                                       const Hypothesis) const override;
    void writeBehaviourDataConstructor(
        std::ostream&,
        const Hypothesis,
        const BehaviourDescription&) const override;
    void writeIntegrationDataConstructor(
        std::ostream&,
        const Hypothesis,
        const BehaviourDescription&) const override;
    void writeBehaviourDataMainVariablesSetters(
        std::ostream&, const BehaviourDescription&) const override;
    void writeIntegrationDataMainVariablesSetters(
        std::ostream&, const BehaviourDescription&) const override;
    void exportMechanicalData(std::ostream&,
                              const Hypothesis,
                              const BehaviourDescription&) const override;
    void writeBehaviourInitializeFunctions(std::ostream&,
                                           const BehaviourDescription&,
                                           const Hypothesis) const override;
    void writeBehaviourPostProcessings(std::ostream&,
                                       const BehaviourDescription&,
                                       const Hypothesis) const override;
    [[nodiscard]] std::string getFunctionNameBasis(
        const std::string&) const override;
    //
    [[nodiscard]] virtual std::string getLibraryName(
        const BehaviourDescription&) const;
    /*!
     * \return the name of the function generated for the given hypothesis
     * \param[in] base name for the function to be generated
     */
    [[nodiscard]] virtual std::string getFunctionNameForHypothesis(
        const std::string&, const Hypothesis) const;
    //! \brief destructor
    ~GenericBehaviourInterface() override;

   protected:
    //! \brief options of the generation of the functions of this interface
    struct GenericFunctionsOptions {
      //! \brief name of the header declaring the functions
      std::string header_file;
      //! \brief header guard
      std::string header_guard;
      //! \brief name of the source file
      std::string source_file;
      /*!
       * \brief if true, the source file also defines the symbols describing
       * the behaviour and the functions changing the parameters and the out
       * of bounds policy
       */
      bool symbols = true;
      /*!
       * \brief if true, the function giving the out of bounds policy has an
       * external linkage, so that the policy is shared with the other
       * sources of the library defining this function in the same way
       */
      bool shared_out_of_bounds_policy = false;
    };
    /*!
     * \brief write the header and the source file of the functions of this
     * interface (integration, initialize functions, post-processings and
     * rotations of a point)
     * \param[in] bd: behaviour description
     * \param[in] fd: file description
     * \param[in] o: options
     */
    void writeGenericFunctions(const BehaviourDescription&,
                               const FileDescription&,
                               const GenericFunctionsOptions&) const;
    /*!
     * \brief write the generation of an MTest file for the given hypothesis
     * \param[out] out: output stream
     * \param[in]  bd: behaviour description
     * \param[in]  h: modelling hypothesis
     */
    virtual void generateMTestFile(std::ostream&,
                                   const BehaviourDescription&,
                                   const Hypothesis) const;
    /*!
     * \return the arguments of the constructors of the behaviour data and
     * integration data classes written by this interface.
     *
     * The data of the integration point are given by an object of type
     * `mfront::gb::BehaviourData` named `mgb_d`. Derived interfaces may add
     * arguments to distinguish their constructors from the ones of this
     * interface.
     */
    [[nodiscard]] virtual std::string getDataConstructorsArguments() const;
    /*!
     * \return the values passed by the behaviour constructor to the
     * constructors of the behaviour data and integration data classes
     */
    [[nodiscard]] virtual std::string getDataConstructorsArgumentsValues()
        const;
    /*!
     * \return the template declaration of the behaviour constructor written
     * by this interface, or an empty string if the constructor is not a
     * template
     */
    [[nodiscard]] virtual std::string
    getBehaviourConstructorTemplateDeclaration() const;
    //! \return the arguments of the behaviour constructor
    [[nodiscard]] virtual std::string getBehaviourConstructorArguments() const;
    /*!
     * \brief write the specialisations of the
     * `mfront::gb::GenericBehaviourTraits` class used by the integration
     * drivers
     * \param[out] out: output stream
     * \param[in] bd: behaviour description
     * \param[in] mhs: modelling hypotheses
     */
    void writeGenericBehaviourTraits(std::ostream&,
                                     const BehaviourDescription&,
                                     const std::set<Hypothesis>&) const;
    //! \brief list of selected hypotheses
    std::optional<std::set<GenericBehaviourInterface::Hypothesis>>
        selectedHypotheses;
  };  // end of struct GenericBehaviourInterface

}  // end of namespace mfront

#endif /* LIB_MFRONT_GENERICBEHAVIOURINTERFACE_HXX */
