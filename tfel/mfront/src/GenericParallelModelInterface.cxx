/*!
 * \file   mfront/src/GenericParallelModelInterface.cxx
 * \brief  This file implements the `generic-parallel` interface for models.
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <fstream>
#include <string_view>
#include <utility>
#include <vector>
#include <algorithm>
#include "TFEL/Raise.hxx"
#include "TFEL/System/System.hxx"
#include "TFEL/Utilities/Data.hxx"
#include "TFEL/Material/ModellingHypothesis.hxx"
#include "MFront/CodeGeneratorUtilities.hxx"
#include "MFront/FileDescription.hxx"
#include "MFront/TargetsDescription.hxx"
#include "MFront/ModelDescription.hxx"
#include "MFront/BehaviourDescription.hxx"
#include "MFront/GenericBehaviourSymbolsGenerator.hxx"
#include "MFront/GenericParallelBehaviourInterface.hxx"
#include "MFront/GenericModelInterface.hxx"
#include "MFront/GenericParallel/Behaviour/BackendFactory.hxx"
#include "MFront/GenericParallel/Behaviour/ParallelSTLBackend.hxx"
#include "MFront/GenericParallelModelInterface.hxx"

namespace mfront {

  /*!
   * \return the name of the structure holding the values of the parameters,
   * declared in the `mfront::gm` namespace
   */
  static std::string getParametersStructureName(const ModelDescription& md) {
    if (!hasGenericModelParametersInitializer(md)) {
      return "mfront::gpb::NoParameters";
    }
    return "mfront::gm::" + md.className + "Parameters";
  }  // end of getParametersStructureName

  std::string GenericParallelModelInterface::getName() {
    return "generic-parallel";
  }  // end of getName

  GenericParallelModelInterface::GenericParallelModelInterface()
      : backend(std::make_unique<
                generic_parallel::behaviour::ParallelSTLBackend>()) {
  }  // end of GenericParallelModelInterface

  void GenericParallelModelInterface::setOptions(const DataMap& opts) {
    using generic_parallel::behaviour::BackendFactory;
    // options of the interface, the other ones being the ones of the backend
    auto backend_options = opts;
    if (contains(opts, "generic_functions")) {
      this->generic_functions = get<bool>(opts, "generic_functions");
      backend_options.erase("generic_functions");
    }
    if (auto b = BackendFactory::get().generate(backend_options);
        b != nullptr) {
      this->backend = std::move(b);
    }
  }  // end of setOptions

  std::pair<bool, GenericParallelModelInterface::tokens_iterator>
  GenericParallelModelInterface::treatKeyword(const std::string& k,
                                              const std::vector<std::string>& i,
                                              tokens_iterator current,
                                              const tokens_iterator) {
    tfel::raise_if(std::find(i.begin(), i.end(), getName()) != i.end(),
                   "GenericParallelModelInterface::treatKeyword: "
                   "unsupported keyword '" +
                       k + "'");
    return {false, current};
  }  // end of treatKeyword

  std::string GenericParallelModelInterface::getLibraryName(
      const ModelDescription& md) const {
    return this->backend->getLibraryName(md.library, md.material, "Model");
  }  // end of getLibraryName

  std::string GenericParallelModelInterface::getHeaderFileName(
      const std::string& name) const {
    return this->backend->getHeaderFileName(name);
  }  // end of getHeaderFileName

  std::string GenericParallelModelInterface::getSourceFileName(
      const std::string& name) const {
    return this->backend->getSourceFileName(name);
  }  // end of getSourceFileName

  void GenericParallelModelInterface::writeOutputFiles(
      const FileDescription& fd, const ModelDescription& md) {
    GenericParallelBehaviourInterface::checkName(
        "GenericParallelModelInterface::writeOutputFiles",
        md.library + md.className);
    this->writeHeaderFile(fd, md);
    this->writeSourceFile(fd, md);
  }  // end of writeOutputFiles

  void GenericParallelModelInterface::writeHeaderFile(
      const FileDescription& fd, const ModelDescription& md) {
    using tfel::material::ModellingHypothesis;
    tfel::system::systemCall::mkdir("include/MFront");
    tfel::system::systemCall::mkdir("include/MFront/GenericParallelModel");
    const auto name = md.library + md.className;
    const auto header = this->getHeaderFileName(name);
    std::ofstream os("include/MFront/GenericParallelModel/" + header);
    tfel::raise_if(!os,
                   "GenericParallelModelInterface::writeHeaderFile: "
                   "could not open file '" +
                       header + "'");
    os.exceptions(std::ios::badbit | std::ios::failbit);
    os << "/*!\n"
       << " * \\file   " << header << '\n'
       << " * \\brief  This file declares the " << getName()
       << " interface for the " << md.className << " model (backend "
       << this->backend->getName() << ")\n"
       << " * \\author " << fd.authorName << '\n'
       << " * \\date   " << fd.date << '\n'
       << " */\n\n";
    const auto hg = "LIB_GENERICPARALLELMODEL_" + makeUpperCase(name) + "_" +
                    makeUpperCase(this->backend->getName()) + "_HXX";
    os << "#ifndef " << hg << "\n"
       << "#define " << hg << "\n\n"
       << "#include\"TFEL/Config/TFELConfig.hxx\"\n"
       << "#include\"MFront/GenericParallelBehaviour/BatchData.h\"\n";
    if (this->generic_functions) {
      os << "#include\"MFront/GenericBehaviour/BehaviourData.h\"\n";
    }
    os << '\n';
    writeExportDirectives(os, false);
    os << "#ifdef __cplusplus\n"
       << "extern \"C\"{\n"
       << "#endif /* __cplusplus */\n\n";
    if (hasGenericModelParametersInitializer(md)) {
      os << "/*!\n"
         << " * \\brief function allowing to modify the parameters of the "
         << name << " model\n"
         << " * \\param[in] n: name of the parameter\n"
         << " * \\param[in] v: value of the parameter\n"
         << " */\n"
         << "MFRONT_SHAREDOBJ int\n"
         << name << "_setParameter(const char *const, const double);\n\n";
    }
    {
      GenericParallelBehaviourInterface bi;
      GenericParallelBehaviourInterface::writeReportFunctionDeclaration(
          os, bi.getFunctionNameBasis(name));
    }
    for (const auto& h : ModellingHypothesis::getModellingHypotheses()) {
      os << "/*!\n"
         << " * \\brief execute the " << name
         << " model over a batch of points for the "
         << ModellingHypothesis::toString(h) << " modelling hypothesis\n"
         << " * \\param[out] s: output status\n"
         << " * \\param[in] d: batch data\n"
         << " * \\return the status, also stored in the output status\n"
         << " */\n"
         << "MFRONT_SHAREDOBJ int "
         << GenericParallelBehaviourInterface::getBatchFunctionName(
                name + "_" + ModellingHypothesis::toString(h))
         << "(mfront_gpb_OutputStatus* const,\n"
         << "const mfront_gpb_BatchData* const);\n\n";
    }
    if (this->generic_functions) {
      // functions of the `generic` interface
      for (const auto& h : ModellingHypothesis::getModellingHypotheses()) {
        os << "/*!\n"
           << " * \\brief execute the " << name << " model at a point for the "
           << ModellingHypothesis::toString(h)
           << " modelling hypothesis, as the `generic` interface\n"
           << " * \\param[in,out] d: data of the point\n"
           << " */\n"
           << "MFRONT_SHAREDOBJ int " << name << "_"
           << ModellingHypothesis::toString(h)
           << "(mfront_gb_BehaviourData* const);\n\n";
      }
    }
    os << "#ifdef __cplusplus\n"
       << "}\n"
       << "#endif /* __cplusplus */\n\n"
       << "#endif /* " << hg << " */\n";
  }  // end of writeHeaderFile

  void GenericParallelModelInterface::writeSourceFile(
      const FileDescription& fd, const ModelDescription& md) {
    using tfel::material::ModellingHypothesis;
    constexpr auto uh = ModellingHypothesis::UNDEFINEDHYPOTHESIS;
    const auto name = md.library + md.className;
    const auto src = this->getSourceFileName(name);
    std::ofstream os("src/" + src);
    os.precision(14);
    tfel::raise_if(!os,
                   "GenericParallelModelInterface::writeSourceFile: "
                   "could not open file '" +
                       src + "'");
    os.exceptions(std::ios::badbit | std::ios::failbit);
    const auto use_qt = md.use_qt.has_value() ? *(md.use_qt) : false;
    os << "/*!\n"
       << " * \\file   " << src << '\n'
       << " * \\brief  This file implements the " << getName()
       << " interface for the " << md.className << " model (backend "
       << this->backend->getName() << ")\n"
       << " * \\author " << fd.authorName << '\n'
       << " * \\date   " << fd.date << '\n'
       << " */\n\n";
    this->backend->writeSpecificIncludes(os);
    os << "#include <cmath>\n"
       << "#include \"MFront/GenericParallel/MemorySpace.hxx\"\n"
       << "#include <cstring>\n"
       << "#include \"TFEL/PhysicalConstants.hxx\"\n"
       << "#include \"TFEL/Config/TFELTypes.hxx\"\n"
       << "#include \"TFEL/Math/Array/View.hxx\"\n"
       << "#include \"TFEL/Math/General/DerivativeType.hxx\"\n";
    if (use_qt) {
      os << "#include \"TFEL/Math/qt.hxx\"\n";
    }
    os << "#include \"TFEL/Material/BoundsCheck.hxx\"\n"
       << "#include \"MFront/GenericBehaviour/Integrate.hxx\"\n"
       << "#include \"MFront/GenericParallelBehaviour/ModelIntegrator.hxx\"\n"
       << "#include \"MFront/GenericParallelModel/"
       << this->getHeaderFileName(name) << "\"\n\n";
    if (!md.includes.empty()) {
      os << md.includes << "\n\n";
    }
    const auto& mhs = ModellingHypothesis::getModellingHypotheses();
    const auto hypotheses =
        std::set<ModellingHypothesis::Hypothesis>{mhs.begin(), mhs.end()};
    const auto bd = convertToBehaviourDescription(md);
    const auto sgp = GenericParallelBehaviourInterface::makeSymbolsGenerator(
        this->generic_functions);
    auto& sg = *sgp;
    // the symbols are the ones of a behaviour of the `generic-parallel`
    // interface: in particular, the interface name tells that the entry
    // points integrate a batch
    GenericParallelBehaviourInterface bi;
    sg.generateGeneralSymbols(os, bi, bd, fd, hypotheses, name);
    this->backend->writeSymbols(os, bi.getFunctionNameBasis(name));
    sg.generateSymbols(os, bi, bd, fd, name, uh);
    if (!md.sources.empty()) {
      os << md.sources << "\n\n";
    }
    const auto parameters = getParametersStructureName(md);
    // the generated classes have internal linkage: two libraries
    // generated for the same model do not share their parameters
    // (static variables of inline functions are unique symbols)
    os << "namespace mfront::gm{\n\n"
       << "namespace {\n\n";
    if (hasGenericModelParametersInitializer(md)) {
      // structure holding the values of the parameters, filled on the host
      // and passed by value to the device
      os << "//! \\brief values of the parameters of the " << md.className
         << " model\n"
         << "struct " << md.className << "Parameters {\n";
      for (const auto& p : md.parameters) {
        if (isGenericModelRealParameter(p)) {
          os << "double " << p.name << ";\n";
        }
      }
      os << "};\n\n";
    }
    writeGenericModelClass(os, fd, md,
                           GenericModelClassOptions{
                               .device = true, .parameters_type = parameters});
    os << "} // end of namespace\n\n"
       << "} // end of namespace mfront::gm\n\n";
    // description of the data of a point
    const auto point_traits = md.className + "PointTraits";
    os << "namespace {\n\n"
       << "//! \\brief description of the data of a point\n"
       << "struct " << point_traits << " {\n"
       << "//! \\brief sizes of the data of a point\n"
       << "static constexpr auto value = mfront::gpb::PointDataSizes{"
       << ".gradients = 0, .thermodynamic_forces = 0, "
       << ".material_properties = " << md.constantMaterialProperties.size()
       << ", .internal_state_variables = " << md.outputs.size()
       << ", .external_state_variables = " << md.inputs.size()
       << ", .tangent_operator = 0};\n";
    // the material properties and the inputs of a model are scalars
    auto write_scalars = [&os](const std::string_view n,
                               const std::size_t size) {
      os << "static constexpr std::array<mfront::gpb::size_type, " << size
         << "> " << n << " = {";
      for (std::size_t i = 0; i != size; ++i) {
        os << (i == 0 ? "1" : ", 1");
      }
      os << "};\n";
    };
    os << "//! \\brief number of components of each material property\n";
    write_scalars("material_properties_sizes",
                  md.constantMaterialProperties.size());
    os << "//! \\brief number of components of each input\n";
    write_scalars("external_state_variables_sizes", md.inputs.size());
    os << "//! \\brief integration strategy\n"
       << "using Integrator = mfront::gpb::ModelIntegrator;\n"
       << "};\n\n"
       << "} // end of anonymous namespace\n\n";
    os << "#ifdef __cplusplus\n"
       << "extern \"C\"{\n"
       << "#endif /* __cplusplus */\n\n";
    GenericParallelBehaviourInterface::writeReportFunction(
        os, bi.getFunctionNameBasis(name));
    const auto initializer = getGenericModelParametersInitializerName(md);
    if (hasGenericModelParametersInitializer(md)) {
      os << "MFRONT_SHAREDOBJ int\n"
         << name << "_setParameter(const char *const n, const double v){\n"
         << "auto& parameters = mfront::gm::" << initializer << "::get();\n";
      auto first = true;
      for (const auto& p : md.parameters) {
        if (!isGenericModelRealParameter(p)) {
          continue;
        }
        os << (first ? "if" : "} else if")  //
           << "(std::strcmp(\"" << p.getExternalName() << "\", n) == 0){\n"
           << "parameters." << p.name << " = "
           << "decltype(parameters." << p.name << "){v};\n"
           << "return 1;\n";
        first = false;
      }
      os << "}\n"
         << "return 0;\n"
         << "} // end of " << name << "_setParameter\n\n";
    }
    // values of the parameters, read on the host
    auto write_parameters = [&os, &md, &initializer, &parameters] {
      if (hasGenericModelParametersInitializer(md)) {
        os << "const auto& mfront_initializer = mfront::gm::" << initializer
           << "::get();\n"
           << "const auto mfront_parameters = " << parameters << "{";
        auto first = true;
        for (const auto& p : md.parameters) {
          if (isGenericModelRealParameter(p)) {
            os << (first ? "" : ", ")
               << "tfel::math::base_type_cast(mfront_initializer." << p.name
               << ")";
            first = false;
          }
        }
        os << "};\n";
      } else {
        os << "const auto mfront_parameters = " << parameters << "{};\n";
      }
    };
    for (const auto& h : ModellingHypothesis::getModellingHypotheses()) {
      const auto f = GenericParallelBehaviourInterface::getBatchFunctionName(
          name + "_" + ModellingHypothesis::toString(h));
      GenericParallelBehaviourInterface::writeBatchFunctionBeginning(
          os, f, "const mfront_gpb_BatchData* const mfront_batch_data",
          ", mfront_batch_data");
      os << "using mfront_point_traits = " << point_traits << ";\n";
      write_parameters();
      // the bounds are checked by the model with its default out of bounds
      // policy, as in the `generic` interface
      os << "const auto mfront_policy = tfel::material::"
         << getDefaultOutOfBoundsPolicyAsString(md) << ";\n";
      this->backend->writeCall(
          os, "integrateBatch",
          "mfront::gm::" + md.className + ", mfront_point_traits",
          "mfront_status", "mfront_data, mfront_parameters, mfront_policy");
      GenericParallelBehaviourInterface::writeBatchFunctionEnd(os, f);
    }
    if (this->generic_functions) {
      // functions of the `generic` interface, with the same parameters as
      // the functions of batches of points
      for (const auto& h : ModellingHypothesis::getModellingHypotheses()) {
        const auto f = name + "_" + ModellingHypothesis::toString(h);
        os << "MFRONT_SHAREDOBJ int " << f
           << "(mfront_gb_BehaviourData* const d){\n"
           << "try{\n";
        write_parameters();
        os << "const mfront::gm::" << md.className
           << " m(*d, mfront_parameters);\n"
           << "if(!m.execute(*d)){\n"
           << "mfront::gb::reportError(*d, \"model execution failed\");\n"
           << "*(d->rdt) = 0.1;\n"
           << "return -1;\n"
           << "}\n"
           << "} catch(...){\n"
           << "mfront::gb::reportFailureByException(*d);\n"
           << "*(d->rdt) = 0.1;\n"
           << "return -1;\n"
           << "}\n"
           << "return 1;\n"
           << "} // end of " << f << "\n\n";
      }
    }
    os << "#ifdef __cplusplus\n"
       << "}\n"
       << "#endif /* __cplusplus */\n\n";
  }  // end of writeSourceFile

  void GenericParallelModelInterface::getTargetsDescription(
      TargetsDescription& td, const ModelDescription& md) {
    using tfel::material::ModellingHypothesis;
    const auto name = md.library + md.className;
    auto entry_points = std::vector<std::string>{};
    for (const auto& h : ModellingHypothesis::getModellingHypotheses()) {
      entry_points.push_back(
          GenericParallelBehaviourInterface::getBatchFunctionName(
              name + "_" + ModellingHypothesis::toString(h)));
    }
    if (this->generic_functions) {
      for (const auto& h : ModellingHypothesis::getModellingHypotheses()) {
        entry_points.push_back(name + "_" + ModellingHypothesis::toString(h));
      }
    }
    GenericParallelBehaviourInterface::describeLibrary(
        td, this->getLibraryName(md), this->getSourceFileName(name),
        "MFront/GenericParallelModel/" + this->getHeaderFileName(name),
        entry_points);
  }  // end of getTargetsDescription

  GenericParallelModelInterface::~GenericParallelModelInterface() = default;

}  // end of namespace mfront
