/*!
 * \file   mfront/src/GenericParallelBehaviourInterface.cxx
 * \brief  This file implements the `generic-parallel` behaviour interface
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <set>
#include <vector>
#include <string_view>
#include <fstream>
#include <sstream>
#include <functional>
#include <algorithm>
#include "TFEL/Raise.hxx"
#include "TFEL/Config/GetInstallPath.hxx"
#include "TFEL/System/System.hxx"
#include "TFEL/Utilities/Data.hxx"
#include "MFront/MFrontLogStream.hxx"
#include "MFront/MFrontUtilities.hxx"
#include "MFront/CodeGeneratorUtilities.hxx"
#include "MFront/FileDescription.hxx"
#include "MFront/LibraryDescription.hxx"
#include "MFront/TargetsDescription.hxx"
#include "MFront/BehaviourDescription.hxx"
#include "MFront/BehaviourVariableDescription.hxx"
#include "MFront/BehaviourRotationFunctions.hxx"
#include "MFront/GenericBehaviourSymbolsGenerator.hxx"
#include "MFront/GenericParallel/Behaviour/AbstractBackend.hxx"
#include "MFront/GenericParallel/Behaviour/BackendFactory.hxx"
#include "MFront/GenericParallel/Behaviour/ParallelSTLBackend.hxx"
#include "MFront/GenericParallelBehaviourInterface.hxx"

namespace mfront {

  //! \brief sizes of the data of an integration point
  struct GenericParallelBehaviourPointDataSizes {
    int gradients = 0;
    int thermodynamic_forces = 0;
    int material_properties = 0;
    int internal_state_variables = 0;
    int external_state_variables = 0;
    int tangent_operator = 0;
  };

  /*!
   * \return the number of values of the variables used by a code block (the
   * inputs of an initialize function or the outputs of a post-processing)
   * \param[in] c: code block
   * \param[in] a: name of the attribute giving the variables
   * \param[in] h: modelling hypothesis
   */
  static int getNumberOfValues(
      const CodeBlock& c,
      const char* const a,
      const tfel::material::ModellingHypothesis::Hypothesis h) {
    auto s = SupportedTypes::TypeSize{};
    if (c.attributes.count(a) != 0) {
      for (const auto& v :
           c.attributes.at(a).get<std::vector<VariableDescription>>()) {
        s += SupportedTypes::getTypeSize(v.type, v.arraySize);
      }
    }
    return s.getValueForModellingHypothesis(h);
  }  // end of getNumberOfValues

  //! \return true if the given behaviour is a finite strain behaviour
  static bool isFiniteStrainBehaviour(const BehaviourDescription& bd) {
    const auto type = bd.getBehaviourType();
    return (type == BehaviourDescription::STANDARDFINITESTRAINBEHAVIOUR) ||
           ((type == BehaviourDescription::STANDARDSTRAINBASEDBEHAVIOUR) &&
            (bd.isStrainMeasureDefined()) &&
            (bd.getStrainMeasure() != BehaviourDescription::LINEARISED));
  }  // end of isFiniteStrainBehaviour

  /*!
   * \return the integration strategy of the given behaviour
   * \param[in] bd: behaviour description
   */
  static std::string getIntegrator(const BehaviourDescription& bd) {
    if (bd.getBehaviourType() ==
        BehaviourDescription::STANDARDFINITESTRAINBEHAVIOUR) {
      return "mfront::gpb::FiniteStrainIntegrator";
    }
    if (isFiniteStrainBehaviour(bd)) {
      const auto ms = bd.getStrainMeasure();
      if (ms == BehaviourDescription::GREENLAGRANGE) {
        return "mfront::gpb::GreenLagrangeStrainIntegrator";
      }
      if (ms != BehaviourDescription::HENCKY) {
        tfel::raise(
            "GenericParallelBehaviourInterface: unsupported strain measure");
      }
      return "mfront::gpb::LogarithmicStrainIntegrator";
    }
    return "mfront::gpb::StandardIntegrator";
  }  // end of getIntegrator

  /*!
   * \brief write the declaration of an array giving the number of
   * components of the variables of a group
   * \param[out] os: output stream
   * \param[in] name: name of the array
   * \param[in] sizes: number of components of each variable
   */
  static void writeVariablesSizes(std::ostream& os,
                                  const std::string_view name,
                                  const std::vector<int>& sizes) {
    os << "static constexpr std::array<mfront::gpb::size_type, "
       << sizes.size() << "> " << name << " = {";
    for (std::size_t i = 0; i != sizes.size(); ++i) {
      os << (i == 0 ? "" : ", ") << sizes[i];
    }
    os << "};\n";
  }  // end of writeVariablesSizes

  /*!
   * \return the number of components of each variable of a list, the
   * elements of an array being distinct variables
   * \param[in] variables: variables
   * \param[in] h: modelling hypothesis
   */
  template <typename VariablesContainer>
  static std::vector<int> getVariablesSizes(
      const VariablesContainer& variables,
      const tfel::material::ModellingHypothesis::Hypothesis h) {
    auto sizes = std::vector<int>{};
    for (const auto& v : variables) {
      const auto s = SupportedTypes::getTypeSize(v.type, 1)
                         .getValueForModellingHypothesis(h);
      sizes.insert(sizes.end(), v.arraySize, s);
    }
    return sizes;
  }  // end of getVariablesSizes

  static GenericParallelBehaviourPointDataSizes getPointDataSizes(
      const GenericParallelBehaviourInterface& i,
      const BehaviourDescription& bd,
      const tfel::material::ModellingHypothesis::Hypothesis h) {
    auto size = [h](const SupportedTypes::TypeSize& s) {
      return s.getValueForModellingHypothesis(h);
    };
    auto variables_size = [](const VariableDescriptionContainer& c) {
      auto s = SupportedTypes::TypeSize{};
      for (const auto& v : c) {
        s += SupportedTypes::getTypeSize(v.type, v.arraySize);
      }
      return s;
    };
    const auto& d = bd.getBehaviourData(h);
    auto r = GenericParallelBehaviourPointDataSizes{};
    const auto mvsizes = bd.getMainVariablesSize();
    r.gradients = size(mvsizes.first);
    r.thermodynamic_forces = size(mvsizes.second);
    // the material properties read by the behaviour start with the elastic
    // properties and the thermal expansion coefficients if the behaviour
    // requires the stiffness tensor or the thermal expansion coefficient
    // tensor. Those properties are not exported in the list of the material
    // properties (they are added by the calling code, as in `MGIS`), but they
    // are part of the material properties of a point.
    const auto mprops = i.buildMaterialPropertiesList(bd, h);
    if (!mprops.first.empty()) {
      const auto& last = mprops.first.back();
      auto s = last.offset;
      s += SupportedTypes::getTypeSize(last.type, last.arraySize);
      r.material_properties = size(s);
    }
    r.internal_state_variables =
        size(variables_size(d.getPersistentVariables()));
    r.external_state_variables =
        size(variables_size(d.getExternalStateVariables()));
    if (isFiniteStrainBehaviour(bd)) {
      const auto T = static_cast<int>(tfel::material::getTensorSize(h));
      r.gradients = T;
      r.thermodynamic_forces = T;
      r.tangent_operator = T * T;
      return r;
    }
    for (const auto& [v1, v2] : bd.getTangentOperatorBlocks()) {
      r.tangent_operator +=
          size(SupportedTypes::getTypeSize(v1.type, v1.arraySize)) *
          size(SupportedTypes::getTypeSize(v2.type, v2.arraySize));
    }
    return r;
  }  // end of getPointDataSizes

  /*!
   * \return the name of the structure holding the values of the parameters
   * for the given hypothesis, or `mfront::gpb::NoParameters` if this
   * structure is not needed.
   * \param[in] bd: behaviour description
   * \param[in] h: modelling hypothesis
   */
  static std::string getParametersStructureName(
      const BehaviourDescription& bd,
      const tfel::material::ModellingHypothesis::Hypothesis h) {
    using tfel::material::ModellingHypothesis;
    if ((areParametersTreatedAsStaticVariables(bd)) ||
        (bd.getBehaviourData(h).getParameters().empty())) {
      return "mfront::gpb::NoParameters";
    }
    return bd.getClassName() + ModellingHypothesis::toString(h) + "Parameters";
  }  // end of getParametersStructureName

  /*!
   * \brief write the structure holding the values of the parameters for the
   * given hypothesis. The names and the types of the members of this
   * structure are the ones of the parameters initializers.
   * \param[out] os: output stream
   * \param[in] bd: behaviour description
   * \param[in] h: modelling hypothesis
   */
  static void writeParametersStructure(
      std::ostream& os,
      const BehaviourDescription& bd,
      const tfel::material::ModellingHypothesis::Hypothesis h) {
    const auto n = getParametersStructureName(bd, h);
    if (n == "mfront::gpb::NoParameters") {
      return;
    }
    os << "//! \\brief values of the parameters\n"
       << "struct " << n << "{\n";
    for (const auto& p : bd.getBehaviourData(h).getParameters()) {
      if (p.type == "int") {
        os << "int " << p.name << ";\n";
      } else if (p.type == "ushort") {
        os << "unsigned short " << p.name << ";\n";
      } else if (p.arraySize == 1) {
        os << "double " << p.name << ";\n";
      } else {
        os << "tfel::math::fsarray<" << p.arraySize << ", double> " << p.name
           << ";\n";
      }
    }
    os << "};\n\n";
  }  // end of writeParametersStructure

  /*!
   * \brief write the checks of the parameters of the behaviour variables
   * (and of their own behaviour variables).
   *
   * In device code, behaviour variables are built by their default
   * constructor, which uses the default values of their parameters: their
   * parameters initializers, which may read a parameters file, are only
   * available on the host. The generated code reports an error if the value
   * of a parameter given by an initializer differs from its default value.
   * \param[out] os: output stream
   * \param[in] bd: behaviour description
   * \param[in] h: modelling hypothesis
   */
  static void writeBehaviourVariablesParametersChecks(
      std::ostream& os,
      const BehaviourDescription& bd,
      const tfel::material::ModellingHypothesis::Hypothesis h) {
    using tfel::material::ModellingHypothesis;
    constexpr auto uh = ModellingHypothesis::UNDEFINEDHYPOTHESIS;
    const auto& d = bd.getBehaviourData(h);
    auto check = [&os, h](const BehaviourVariableDescription& v,
                          const std::string& kind) {
      const auto& ib = v.behaviour;
      if (!areParametersTreatedAsStaticVariables(ib)) {
        const auto ns = ib.getAttribute<std::string>(
            BehaviourDescription::internalNamespace, "");
        const auto prefix = ns.empty() ? std::string{"tfel::material::"}
                                       : "tfel::material::" + ns + "::";
        // default values are written as in the code generated for the
        // behaviour variable
        const auto precision = os.precision(14);
        for (const auto& p : ib.getBehaviourData(h).getParameters()) {
          const auto initializer =
              prefix + ib.getClassName() +
              (ib.hasParameter(uh, p.name) ? ""
                                           : ModellingHypothesis::toString(h)) +
              "ParametersInitializer::get()." + p.name;
          auto conditions = std::vector<std::string>{};
          if (p.type == "int") {
            conditions.push_back(
                initializer + " != " +
                std::to_string(ib.getIntegerParameterDefaultValue(h, p.name)));
          } else if (p.type == "ushort") {
            conditions.push_back(
                initializer + " != " +
                std::to_string(
                    ib.getUnsignedShortParameterDefaultValue(h, p.name)));
          } else {
            for (unsigned short i = 0; i != p.arraySize; ++i) {
              auto c = std::ostringstream{};
              c.copyfmt(os);
              c << initializer
                << (p.arraySize == 1 ? "" : "[" + std::to_string(i) + "]")
                << " != double("
                << (p.arraySize == 1
                        ? ib.getFloattingPointParameterDefaultValue(h, p.name)
                        : ib.getFloattingPointParameterDefaultValue(h, p.name,
                                                                    i))
                << ")";
              conditions.push_back(c.str());
            }
          }
          for (const auto& c : conditions) {
            os << "if (" << c << ") {\n"
               << "return mfront::gpb::report(mfront_status, -2, \"the "
               << "parameter '" << p.name << "' of the " << kind << " '"
               << v.name << "' differs from its default value, which is "
               << "used in device code\");\n"
               << "}\n";
          }
        }
        os.precision(precision);
      }
      writeBehaviourVariablesParametersChecks(os, ib, h);
    };
    for (const auto& v : d.getBehaviourVariables()) {
      check(v, "behaviour variable");
    }
    for (const auto& v : d.getBehaviourVariableFactories()) {
      check(v, "behaviour variable factory");
    }
  }  // end of writeBehaviourVariablesParametersChecks

  /*!
   * \brief write the initialisation of the values of the parameters from
   * the parameters initializers, which are only available on the host.
   * \param[out] os: output stream
   * \param[in] bd: behaviour description
   * \param[in] h: modelling hypothesis
   */
  static void writeParametersInitialisation(
      std::ostream& os,
      const BehaviourDescription& bd,
      const tfel::material::ModellingHypothesis::Hypothesis h) {
    using tfel::material::ModellingHypothesis;
    constexpr auto uh = ModellingHypothesis::UNDEFINEDHYPOTHESIS;
    const auto n = getParametersStructureName(bd, h);
    os << "auto mfront_parameters = " << n << "{};\n";
    if (n == "mfront::gpb::NoParameters") {
      return;
    }
    for (const auto& p : bd.getBehaviourData(h).getParameters()) {
      const auto initializer = [&bd, &p, h] {
        if (bd.hasParameter(uh, p.name)) {
          return bd.getClassName() + "ParametersInitializer";
        }
        return bd.getClassName() + ModellingHypothesis::toString(h) +
               "ParametersInitializer";
      }();
      os << "mfront_parameters." << p.name << " = " << initializer << "::get()."
         << p.name << ";\n";
    }
  }  // end of writeParametersInitialisation

  /*!
   * \return the name of the function object rotating the values of a point
   * \param[in] bd: behaviour description
   * \param[in] h: modelling hypothesis
   * \param[in] f: rotation function
   */
  static std::string getPointRotationName(
      const BehaviourDescription& bd,
      const tfel::material::ModellingHypothesis::Hypothesis h,
      const BehaviourRotationFunction& f) {
    return bd.getClassName() +
           tfel::material::ModellingHypothesis::toString(h) + "_" + f.name;
  }  // end of getPointRotationName

  /*!
   * \brief write the function objects rotating the values of a point of an
   * orthotropic behaviour
   * \param[out] os: output stream
   * \param[in] bd: behaviour description
   * \param[in] h: modelling hypothesis
   */
  static void writePointRotations(
      std::ostream& os,
      const BehaviourDescription& bd,
      const tfel::material::ModellingHypothesis::Hypothesis h) {
    for (const auto& f : getBehaviourRotationFunctions(bd, h)) {
      if (!f.error.empty()) {
        continue;
      }
      os << "//! \\brief rotate " << f.description << " at one point\n"
         << "struct " << getPointRotationName(bd, h, f) << " {\n"
         << "//! \\brief number of values of a point\n"
         << "static constexpr mfront::gpb::size_type size = " << f.size << ";\n"
         << "TFEL_HOST_DEVICE void operator()("
         << "mfront::gb::real* const dest, "
         << "const mfront::gb::real* const src, "
         << "const mfront::gb::real* const rv) const noexcept {\n";
      writePointRotation(os, f);
      os << "}\n"
         << "};\n\n";
    }
  }  // end of writePointRotations

  //! \return the name of the function object calling an initialize function
  static std::string getInitializeFunctionCallerName(
      const BehaviourDescription& bd, const std::string& f) {
    return bd.getClassName() + "InitializeFunction_" + f;
  }  // end of getInitializeFunctionCallerName

  //! \return the name of the function object calling a post-processing
  static std::string getPostProcessingCallerName(const BehaviourDescription& bd,
                                                 const std::string& f) {
    return bd.getClassName() + "PostProcessing_" + f;
  }  // end of getPostProcessingCallerName

  /*!
   * \brief write the function objects calling the initialize functions and
   * the post-processings of the behaviour by their names, since calls
   * through pointers to members are not supported in device code by some
   * compilers.
   * \param[out] os: output stream
   * \param[in] bd: behaviour description
   * \param[in] mhs: modelling hypotheses
   */
  static void writeMethodCallers(
      std::ostream& os,
      const BehaviourDescription& bd,
      const std::set<tfel::material::ModellingHypothesis::Hypothesis>& mhs) {
    auto initialize_functions = std::set<std::string>{};
    auto postprocessings = std::set<std::string>{};
    for (const auto h : mhs) {
      const auto& d = bd.getBehaviourData(h);
      for (const auto& f : d.getInitializeFunctions()) {
        initialize_functions.insert(f.first);
      }
      for (const auto& f : d.getPostProcessings()) {
        postprocessings.insert(f.first);
      }
    }
    for (const auto& f : initialize_functions) {
      os << "//! \\brief function object calling the initialize function '" << f
         << "'\n"
         << "struct " << getInitializeFunctionCallerName(bd, f) << " {\n"
         << "template <typename Behaviour>\n"
         << "TFEL_HOST_DEVICE void operator()(Behaviour& b, "
         << "const mfront::gb::real* const values) const {\n"
         << "b.execute" << f << "InitializeFunction(values);\n"
         << "}\n"
         << "};\n\n";
    }
    for (const auto& f : postprocessings) {
      os << "//! \\brief function object calling the post-processing '" << f
         << "'\n"
         << "struct " << getPostProcessingCallerName(bd, f) << " {\n"
         << "template <typename Behaviour>\n"
         << "TFEL_HOST_DEVICE void operator()(Behaviour& b, "
         << "mfront::gb::real* const outputs, "
         << "const typename Behaviour::BehaviourData& initial_state) const {\n"
         << "b.execute" << f << "PostProcessing(outputs, initial_state);\n"
         << "}\n"
         << "};\n\n";
    }
  }  // end of writeMethodCallers

  std::string GenericParallelBehaviourInterface::getName() {
    return "generic-parallel";
  }  // end of getName

  GenericParallelBehaviourInterface::GenericParallelBehaviourInterface()
      : backend(std::make_unique<
                generic_parallel::behaviour::ParallelSTLBackend>()) {
  }  // end of GenericParallelBehaviourInterface

  std::string GenericParallelBehaviourInterface::getInterfaceName() const {
    return "GenericParallel";
  }  // end of getInterfaceName

  void GenericParallelBehaviourInterface::setOptions(const DataMap& opts) {
    using generic_parallel::behaviour::BackendFactory;
    using tfel::utilities::DataStructure;
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addDataTypeValidator<DataStructure>("backend");
    validator.validate(opts);
    if (!contains(opts, "backend")) {
      return;
    }
    const auto ds = DataStructure::convert(get(opts, "backend"));
    if (getVerboseMode() >= VERBOSE_DEBUG) {
      getLogStream() << "GenericParallelBehaviourInterface::setOptions: "
                        "creating backend '"
                     << ds.name << "'\n";
    }
    this->backend = BackendFactory::get().generate(ds.name, ds.data);
  }  // end of setOptions

  std::pair<bool, GenericParallelBehaviourInterface::tokens_iterator>
  GenericParallelBehaviourInterface::treatKeyword(
      BehaviourDescription& bd,
      const std::string& k,
      const std::vector<std::string>& i,
      tokens_iterator current,
      const tokens_iterator end) {
    const auto is_hypotheses_selection =
        (k == "@SelectedModellingHypothesis") ||
        (k == "@SelectedModellingHypotheses");
    if (!i.empty()) {
      if (std::find(i.begin(), i.end(), getName()) == i.end()) {
        return {false, current};
      }
      if (!is_hypotheses_selection) {
        tfel::raise(
            "GenericParallelBehaviourInterface::treatKeyword: "
            "unsupported keyword '" +
            k + "'");
      }
    }
    if (is_hypotheses_selection) {
      return GenericBehaviourInterface::treatKeyword(bd, k, {}, current, end);
    }
    return {false, current};
  }  // end of treatKeyword

  std::string
  GenericParallelBehaviourInterface::getParametersSourceInBehaviourConstructor(
      const BehaviourDescription& bd, const Hypothesis) const {
    if (areParametersTreatedAsStaticVariables(bd)) {
      return "";
    }
    return "mgpb_parameters";
  }  // end of getParametersSourceInBehaviourConstructor

  bool
  GenericParallelBehaviourInterface::requiresDeviceCompatibleBehaviourClass()
      const {
    return true;
  }  // end of requiresDeviceCompatibleBehaviourClass

  std::string GenericParallelBehaviourInterface::
      getBehaviourConstructorTemplateDeclaration() const {
    return "template <typename ParametersType>";
  }  // end of getBehaviourConstructorTemplateDeclaration

  std::string
  GenericParallelBehaviourInterface::getBehaviourConstructorArguments() const {
    // as in the `generic` interface, the names of the arguments start with
    // the `mgb_` prefix (`mgpb_` for the arguments specific to this
    // interface), whereas the local variables of the generated functions
    // start with the `mfront_` prefix
    return "const mfront::gb::BehaviourData& mgb_d, "
           "[[maybe_unused]] const ParametersType& mgpb_parameters";
  }  // end of getBehaviourConstructorArguments

  std::string GenericParallelBehaviourInterface::getLibraryName(
      const BehaviourDescription& bd) const {
    const auto base = [&bd]() -> std::string {
      if (!bd.getLibrary().empty()) {
        return bd.getLibrary();
      }
      if (!bd.getMaterialName().empty()) {
        return bd.getMaterialName();
      }
      return "Behaviour";
    }();
    return this->getInterfaceName() + base + '-' + this->backend->getName();
  }  // end of getLibraryName

  std::string GenericParallelBehaviourInterface::getHeaderFileName(
      const std::string& name) const {
    return name + '-' + getName() + '-' + this->backend->getName() + ".hxx";
  }  // end of getHeaderFileName

  std::string GenericParallelBehaviourInterface::getSourceFileName(
      const std::string& name) const {
    return name + '-' + getName() + '-' + this->backend->getName() + '.' +
           this->backend->getSourceFileExtension();
  }  // end of getSourceFileName

  void GenericParallelBehaviourInterface::getTargetsDescription(
      TargetsDescription& d, const BehaviourDescription& bd) {
    const auto lib = this->getLibraryName(bd);
    const auto name = bd.getLibrary() + bd.getClassName();
    const auto tfel_config = tfel::getTFELConfigExecutableName();
    auto& l = d.getLibrary(lib);
    insert_if(l.sources, this->getSourceFileName(name));
    d.headers.push_back("MFront/GenericParallelBehaviour/" +
                        this->getHeaderFileName(name));
    insert_if(l.link_directories,
              "$(shell " + tfel_config + " --library-path)");
    insert_if(l.link_libraries, "$(shell " + tfel_config +
                                    " --library-dependency "
                                    "--material --mfront-profiling)");
    for (const auto h : this->getModellingHypothesesToBeTreated(bd)) {
      insert_if(l.epts, this->getFunctionNameForHypothesis(name, h));
    }
  }  // end of getTargetsDescription

  void GenericParallelBehaviourInterface::endTreatment(
      const BehaviourDescription& bd, const FileDescription& fd) const {
    tfel::system::systemCall::mkdir("include/MFront");
    tfel::system::systemCall::mkdir("include/MFront/GenericParallelBehaviour");
    this->writeHeaderFile(bd, fd);
    this->writeSourceFile(bd, fd);
  }  // end of endTreatment

  void GenericParallelBehaviourInterface::writeHeaderFile(
      const BehaviourDescription& bd, const FileDescription& fd) const {
    const auto name = bd.getLibrary() + bd.getClassName();
    const auto header = this->getHeaderFileName(name);
    const auto fn = "include/MFront/GenericParallelBehaviour/" + header;
    std::ofstream os(fn);
    tfel::raise_if(!os,
                   "GenericParallelBehaviourInterface::writeHeaderFile: "
                   "could not open file '" +
                       fn + "'");
    os.exceptions(std::ios::badbit | std::ios::failbit);
    os << "/*!\n"
       << "* \\file   " << header << '\n'
       << "* \\brief  This file declares the " << getName()
       << " interface for the " << bd.getClassName() << " behaviour (backend "
       << this->backend->getName() << ")\n"
       << "* \\author " << fd.authorName << '\n'
       << "* \\date   " << fd.date << '\n'
       << "*/\n\n";
    const auto hg = this->getHeaderGuard(bd) + "_" +
                    makeUpperCase(this->backend->getName());
    os << "#ifndef " << hg << "\n"
       << "#define " << hg << "\n\n"
       << "#include\"TFEL/Config/TFELConfig.hxx\"\n"
       << "#include\"MFront/GenericParallelBehaviour/BatchData.h\"\n\n";
    this->writeVisibilityDefines(os);
    os << "#ifdef __cplusplus\n"
       << "extern \"C\"{\n"
       << "#endif /* __cplusplus */\n\n";
    this->writeSetOutOfBoundsPolicyFunctionDeclaration(os, name);
    this->writeSetParametersFunctionsDeclarations(os, bd, name);
    for (const auto h : this->getModellingHypothesesToBeTreated(bd)) {
      os << "/*!\n"
         << " * \\brief integrate the behaviour over a batch of integration "
            "points\n"
         << " * \\param[out] s: output status\n"
         << " * \\param[in] d: batch data\n"
         << " * \\return the status, also stored in the output status\n"
         << " */\n"
         << "MFRONT_SHAREDOBJ int "
         << this->getFunctionNameForHypothesis(name, h)
         << "(mfront_gpb_OutputStatus* const,\n"
         << "const mfront_gpb_BatchData* const);\n\n";
      const auto& d = bd.getBehaviourData(h);
      for (const auto& ifct : d.getInitializeFunctions()) {
        os << "/*!\n"
           << " * \\brief execute the initialize function '" << ifct.first
           << "' over a batch of integration points\n"
           << " * \\param[out] s: output status\n"
           << " * \\param[in] d: batch data\n"
           << " * \\param[in] inputs: inputs of the initialize function, which "
           << "may be null if the initialize function has no input\n"
           << " * \\return the status, also stored in the output status\n"
           << " */\n"
           << "MFRONT_SHAREDOBJ int "
           << this->getFunctionNameForHypothesis(name, h)
           << "_InitializeFunction_" << ifct.first
           << "(mfront_gpb_OutputStatus* const,\n"
           << "const mfront_gpb_BatchData* const,\n"
           << "const mfront_gpb_ConstField* const);\n\n";
      }
      for (const auto& pp : d.getPostProcessings()) {
        os << "/*!\n"
           << " * \\brief execute the post-processing '" << pp.first
           << "' over a batch of integration points\n"
           << " * \\param[out] s: output status\n"
           << " * \\param[out] outputs: outputs of the post-processing\n"
           << " * \\param[in] d: batch data\n"
           << " * \\return the status, also stored in the output status\n"
           << " */\n"
           << "MFRONT_SHAREDOBJ int "
           << this->getFunctionNameForHypothesis(name, h) << "_PostProcessing_"
           << pp.first << "(mfront_gpb_OutputStatus* const,\n"
           << "const mfront_gpb_Field* const,\n"
           << "const mfront_gpb_BatchData* const);\n\n";
      }
      for (const auto& f : getBehaviourRotationFunctions(bd, h)) {
        os << "/*!\n"
           << " * \\brief rotate " << f.description
           << " over a batch of integration points\n"
           << " * \\param[out] s: output status\n"
           << " * \\param[in] d: rotation data\n"
           << " * \\return the status, also stored in the output status\n"
           << " */\n"
           << "MFRONT_SHAREDOBJ int "
           << this->getFunctionNameForHypothesis(name, h) << "_" << f.name
           << "(mfront_gpb_OutputStatus* const,\n"
           << "const mfront_gpb_RotationData* const);\n\n";
      }
    }
    os << "#ifdef __cplusplus\n"
       << "}\n"
       << "#endif /* __cplusplus */\n\n"
       << "#endif /* " << hg << " */\n";
  }  // end of writeHeaderFile

  void GenericParallelBehaviourInterface::writeSourceFile(
      const BehaviourDescription& bd, const FileDescription& fd) const {
    using tfel::material::ModellingHypothesis;
    const auto mhs = this->getModellingHypothesesToBeTreated(bd);
    const auto name = bd.getLibrary() + bd.getClassName();
    const auto src = this->getSourceFileName(name);
    const auto fn = "src/" + src;
    std::ofstream os(fn);
    tfel::raise_if(!os,
                   "GenericParallelBehaviourInterface::writeSourceFile: "
                   "could not open file '" +
                       fn + "'");
    os.exceptions(std::ios::badbit | std::ios::failbit);
    os << "/*!\n"
       << "* \\file   " << src << '\n'
       << "* \\brief  This file implements the " << getName()
       << " interface for the " << bd.getClassName() << " behaviour (backend "
       << this->backend->getName() << ")\n"
       << "* \\author " << fd.authorName << '\n'
       << "* \\date   " << fd.date << '\n'
       << "*/\n\n";
    this->writeVisibilityDefines(os);
    os << "#include\"TFEL/Material/OutOfBoundsPolicy.hxx\"\n"
       << "#include\"TFEL/Material/" << bd.getClassName() << ".hxx\"\n"
       << "#include\"MFront/GenericBehaviour/GenericBehaviourTraits.hxx\"\n"
       << "#include\"MFront/GenericParallelBehaviour/Integrate.hxx\"\n";
    if (isFiniteStrainBehaviour(bd)) {
      os << "#include\"MFront/GenericParallelBehaviour/"
         << "FiniteStrainIntegrators.hxx\"\n";
    }
    this->backend->writeSpecificIncludes(os);
    os << "#include\"MFront/GenericParallelBehaviour/"
       << this->getHeaderFileName(name) << "\"\n\n";
    this->writeGetOutOfBoundsPolicyFunctionImplementation(os, bd, name);
    this->writeGenericBehaviourTraits(os, bd, mhs);
    os << "namespace {\n\n";
    writeMethodCallers(os, bd, mhs);
    for (const auto h : mhs) {
      const auto sizes = getPointDataSizes(*this, bd, h);
      os << "//! \\brief description of the data of an integration point\n"
         << "struct " << bd.getClassName() << ModellingHypothesis::toString(h)
         << "PointTraits {\n"
         << "//! \\brief maximal sizes of the data of an integration point\n"
         << "static constexpr auto value = mfront::gpb::PointDataSizes{"
         << ".gradients = " << sizes.gradients << ", "
         << ".thermodynamic_forces = " << sizes.thermodynamic_forces << ", "
         << ".material_properties = " << sizes.material_properties << ", "
         << ".internal_state_variables = " << sizes.internal_state_variables
         << ", "
         << ".external_state_variables = " << sizes.external_state_variables
         << ", "
         << ".tangent_operator = " << sizes.tangent_operator << "};\n";
      // the material properties of a point (including the elastic
      // properties added when the behaviour requires the stiffness tensor)
      // and the external state variables, variable by variable
      os << "//! \\brief number of components of each material property\n";
      writeVariablesSizes(
          os, "material_properties_sizes",
          getVariablesSizes(this->buildMaterialPropertiesList(bd, h).first,
                            h));
      os << "//! \\brief number of components of each external state "
            "variable\n";
      writeVariablesSizes(
          os, "external_state_variables_sizes",
          getVariablesSizes(
              bd.getBehaviourData(h).getExternalStateVariables(), h));
      os << "//! \\brief integration strategy\n"
         << "using Integrator = " << getIntegrator(bd) << ";\n"
         << "};\n\n";
      writeParametersStructure(os, bd, h);
      writePointRotations(os, bd, h);
    }
    os << "} // end of anonymous namespace\n\n";
    os << "#ifdef __cplusplus\n"
       << "extern \"C\"{\n"
       << "#endif /* __cplusplus */\n\n";
    GenericBehaviourSymbolsGenerator sg;
    sg.generateGeneralSymbols(os, *this, bd, fd, mhs, name);
    this->backend->writeSymbols(os, this->getFunctionNameBasis(name));
    if (!bd.areAllMechanicalDataSpecialised(mhs)) {
      const auto uh = ModellingHypothesis::UNDEFINEDHYPOTHESIS;
      sg.generateSymbols(os, *this, bd, fd, name, uh);
    }
    for (const auto& h : mhs) {
      if (bd.hasSpecialisedMechanicalData(h)) {
        sg.generateSymbols(os, *this, bd, fd, name, h);
      }
    }
    this->writeSetOutOfBoundsPolicyFunctionImplementation(os, bd, name);
    this->writeSetParametersFunctionsImplementations(os, bd, name);
    // write a function of the generated library
    // \param[in] h: modelling hypothesis
    // \param[in] f: name of the function
    // \param[in] arguments: additional arguments of the function
    // \param[in] captures: captures of the additional arguments
    // \param[in] call: function writing the call to the backend
    auto write_function = [this, &os, &bd, &name](
                              const Hypothesis h, const std::string& f,
                              const std::string& arguments,
                              const std::string& captures,
                              const std::function<void()>& call) {
      os << "MFRONT_SHAREDOBJ int " << f
         << "(mfront_gpb_OutputStatus* const mfront_output_status,\n"
         << arguments << "){\n"
         << "return mfront::gpb::handleExceptions(mfront_output_status, "
         << "[mfront_output_status" << captures << "] {\n"
         << "if (mfront_batch_data == nullptr) {\n"
         << "mfront::gpb::initialize(*mfront_output_status, 1);\n"
         << "return mfront::gpb::report(*mfront_output_status, -2, "
         << "\"no batch data given\");\n"
         << "}\n"
         << "using namespace tfel::material;\n"
         << "using real = mfront::gb::real;\n"
         << "constexpr auto h = ModellingHypothesis::"
         << ModellingHypothesis::toUpperCaseString(h) << ";\n";
      if (bd.useQt()) {
        os << "using Behaviour = " << bd.getClassName() << "<h, real, true>;\n";
      } else {
        os << "using Behaviour = " << bd.getClassName()
           << "<h, real, false>;\n";
      }
      os << "using mfront_point_traits = " << bd.getClassName()
         << ModellingHypothesis::toString(h) << "PointTraits;\n"
         << "auto& mfront_status = *mfront_output_status;\n"
         << "const auto& mfront_data = *mfront_batch_data;\n";
      writeParametersInitialisation(os, bd, h);
      writeBehaviourVariablesParametersChecks(os, bd, h);
      os << "const auto mfront_policy = " << this->getFunctionNameBasis(name)
         << "_getOutOfBoundsPolicy();\n";
      call();
      os << "});\n"
         << "} // end of " << f << "\n\n";
    };
    for (const auto h : mhs) {
      const auto f = this->getFunctionNameForHypothesis(name, h);
      write_function(h, f,
                     "const mfront_gpb_BatchData* const mfront_batch_data",
                     ", mfront_batch_data", [this, &os] {
                       this->backend->writeCall(
                           os, "integrateBatch",
                           "Behaviour, mfront_point_traits", "mfront_status",
                           "mfront_data, mfront_parameters, mfront_policy");
                     });
      const auto& d = bd.getBehaviourData(h);
      for (const auto& [ifct, c] : d.getInitializeFunctions()) {
        const auto n = getNumberOfValues(
            c, CodeBlock::used_initialize_function_variables, h);
        write_function(h, f + "_InitializeFunction_" + ifct,
                       "const mfront_gpb_BatchData* const mfront_batch_data,\n"
                       "const mfront_gpb_ConstField* const mfront_inputs",
                       ", mfront_batch_data, mfront_inputs",
                       [this, &os, &bd, &ifct, n] {
                         // the inputs may be null if the function has no input
                         this->backend->writeCall(
                             os, "executeInitializeFunctionOnBatch",
                             "Behaviour, mfront_point_traits, " +
                                 getInitializeFunctionCallerName(bd, ifct) +
                                 ", " + std::to_string(n),
                             "mfront_status",
                             "mfront_data, "
                             "mfront::gpb::getField(mfront_inputs), "
                             "mfront_parameters, mfront_policy");
                       });
      }
      for (const auto& [pp, c] : d.getPostProcessings()) {
        const auto n =
            getNumberOfValues(c, CodeBlock::used_postprocessing_variables, h);
        write_function(
            h, f + "_PostProcessing_" + pp,
            "const mfront_gpb_Field* const mfront_outputs,\n"
            "const mfront_gpb_BatchData* const mfront_batch_data",
            ", mfront_outputs, mfront_batch_data", [this, &os, &bd, &pp, n] {
              this->backend->writeCall(
                  os, "executePostProcessingOnBatch",
                  "Behaviour, mfront_point_traits, " +
                      getPostProcessingCallerName(bd, pp) + ", true, " +
                      std::to_string(n),
                  "mfront_status",
                  "mfront::gpb::getField(mfront_outputs), "
                  "mfront_data, mfront_parameters, mfront_policy");
            });
      }
    }
    for (const auto h : mhs) {
      this->writeRotationFunctions(os, bd, name, h);
    }
    os << "#ifdef __cplusplus\n"
       << "}\n"
       << "#endif /* __cplusplus */\n\n";
  }  // end of writeSourceFile

  void GenericParallelBehaviourInterface::writeRotationFunctions(
      std::ostream& os,
      const BehaviourDescription& bd,
      const std::string& name,
      const Hypothesis h) const {
    const auto fn = this->getFunctionNameForHypothesis(name, h);
    for (const auto& f : getBehaviourRotationFunctions(bd, h)) {
      os << "MFRONT_SHAREDOBJ int " << fn << "_" << f.name
         << "(mfront_gpb_OutputStatus* const mfront_output_status,\n"
         << "const mfront_gpb_RotationData* const mfront_rotation_data){\n"
         << "return mfront::gpb::handleExceptions(mfront_output_status, "
         << "[mfront_output_status, mfront_rotation_data] {\n"
         << "auto& mfront_status = *mfront_output_status;\n"
         << "if (mfront_rotation_data == nullptr) {\n"
         << "mfront::gpb::initialize(mfront_status, 1);\n"
         << "return mfront::gpb::report(mfront_status, -2, "
         << "\"no rotation data given\");\n"
         << "}\n";
      if (f.error.empty()) {
        os << "const auto& mfront_data = *mfront_rotation_data;\n";
        this->backend->writeCall(os, "rotateBatch",
                                 getPointRotationName(bd, h, f),
                                 "mfront_status", "mfront_data");
      } else {
        os << "mfront::gpb::initialize(mfront_status, 1);\n"
           << "return mfront::gpb::report(mfront_status, -2, \"" << f.error
           << "\");\n";
      }
      os << "});\n"
         << "} // end of " << fn << "_" << f.name << "\n\n";
    }
  }  // end of writeRotationFunctions

  GenericParallelBehaviourInterface::~GenericParallelBehaviourInterface() =
      default;

}  // end of namespace mfront
