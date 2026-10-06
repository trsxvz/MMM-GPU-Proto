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
#include <algorithm>
#include "TFEL/Raise.hxx"
#include "TFEL/Config/GetInstallPath.hxx"
#include "TFEL/System/System.hxx"
#include "TFEL/Utilities/Data.hxx"
#include "TFEL/Material/ModellingHypothesis.hxx"
#include "MFront/MFrontLogStream.hxx"
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
    using tfel::utilities::DataStructure;
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addDataTypeValidator<DataStructure>("backend");
    validator.validate(opts);
    if (!contains(opts, "backend")) {
      return;
    }
    const auto ds = DataStructure::convert(get(opts, "backend"));
    if (getVerboseMode() >= VERBOSE_DEBUG) {
      getLogStream() << "GenericParallelModelInterface::setOptions: "
                        "creating backend '"
                     << ds.name << "'\n";
    }
    this->backend = BackendFactory::get().generate(ds.name, ds.data);
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

  void GenericParallelModelInterface::declareReservedNames(
      std::set<std::string>&) {}  // end of declareReservedNames

  std::string GenericParallelModelInterface::getLibraryName(
      const ModelDescription& md) const {
    const auto base = [&md]() -> std::string {
      if (!md.library.empty()) {
        return md.library;
      }
      if (!md.material.empty()) {
        return md.material;
      }
      return "Model";
    }();
    return "GenericParallel" + base + '-' + this->backend->getName();
  }  // end of getLibraryName

  std::string GenericParallelModelInterface::getHeaderFileName(
      const std::string& name) const {
    return name + '-' + getName() + '-' + this->backend->getName() + ".hxx";
  }  // end of getHeaderFileName

  std::string GenericParallelModelInterface::getSourceFileName(
      const std::string& name) const {
    return name + '-' + getName() + '-' + this->backend->getName() + '.' +
           this->backend->getSourceFileExtension();
  }  // end of getSourceFileName

  void GenericParallelModelInterface::writeOutputFiles(
      const FileDescription& fd, const ModelDescription& md) {
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
       << "#include\"MFront/GenericParallelBehaviour/BatchData.h\"\n\n";
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
    for (const auto& h : ModellingHypothesis::getModellingHypotheses()) {
      os << "/*!\n"
         << " * \\brief execute the " << name
         << " model over a batch of points for the "
         << ModellingHypothesis::toString(h) << " modelling hypothesis\n"
         << " * \\param[out] s: output status\n"
         << " * \\param[in] d: batch data\n"
         << " * \\return the status, also stored in the output status\n"
         << " */\n"
         << "MFRONT_SHAREDOBJ int " << name << "_"
         << ModellingHypothesis::toString(h)
         << "(mfront_gpb_OutputStatus* const,\n"
         << "const mfront_gpb_BatchData* const);\n\n";
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
    const auto mhs = ModellingHypothesis::getModellingHypotheses();
    const auto hypotheses =
        std::set<ModellingHypothesis::Hypothesis>{mhs.begin(), mhs.end()};
    const auto bd = convertToBehaviourDescription(md);
    GenericBehaviourSymbolsGenerator sg;
    // the symbols are the ones of a behaviour of the `generic-parallel`
    // interface: in particular, the interface name tells that the entry
    // points integrate a batch
    GenericParallelBehaviourInterface bi;
    sg.generateGeneralSymbols(os, bi, bd, fd, hypotheses, name);
    os << "MFRONT_EXPORT_SYMBOL(const char*, "
       << bi.getFunctionNameBasis(name) << "_generic_parallel_backend, \""
       << this->backend->getName() << "\");\n";
    os << '\n';
    sg.generateSymbols(os, bi, bd, fd, name, uh);
    if (!md.sources.empty()) {
      os << md.sources << "\n\n";
    }
    const auto parameters = getParametersStructureName(md);
    os << "namespace mfront::gm{\n\n";
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
    os << "} // end of namespace mfront::gm\n\n";
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
    for (const auto& h : ModellingHypothesis::getModellingHypotheses()) {
      const auto f = name + "_" + ModellingHypothesis::toString(h);
      os << "MFRONT_SHAREDOBJ int " << f
         << "(mfront_gpb_OutputStatus* const mfront_output_status,\n"
         << "const mfront_gpb_BatchData* const mfront_batch_data){\n"
         << "return mfront::gpb::handleExceptions(mfront_output_status, "
         << "[mfront_output_status, mfront_batch_data] {\n"
         << "if (mfront_batch_data == nullptr) {\n"
         << "return mfront::gpb::report(*mfront_output_status, -2, "
         << "\"no batch data given\");\n"
         << "}\n"
         << "using mfront_point_traits = " << point_traits << ";\n"
         << "auto& mfront_status = *mfront_output_status;\n"
         << "const auto& mfront_data = *mfront_batch_data;\n";
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
      // unused: the bounds are checked with the default out of bounds policy
      // of the model, as in the `generic` interface
      os << "const auto mfront_policy = tfel::material::None;\n";
      this->backend->writeCall(
          os, "integrate",
          "mfront::gm::" + md.className + ", mfront_point_traits",
          "mfront_status, mfront_data, "
          "mfront_parameters, mfront_policy");
      os << "});\n"
         << "} // end of " << f << "\n\n";
    }
    os << "#ifdef __cplusplus\n"
       << "}\n"
       << "#endif /* __cplusplus */\n\n";
  }  // end of writeSourceFile

  void GenericParallelModelInterface::getTargetsDescription(
      TargetsDescription& td, const ModelDescription& md) {
    using tfel::material::ModellingHypothesis;
    const auto lib = this->getLibraryName(md);
    const auto name = md.library + md.className;
    const auto tfel_config = tfel::getTFELConfigExecutableName();
    auto& l = td.getLibrary(lib);
    insert_if(l.sources, this->getSourceFileName(name));
    td.headers.push_back("MFront/GenericParallelModel/" +
                         this->getHeaderFileName(name));
    insert_if(l.link_directories,
              "$(shell " + tfel_config + " --library-path)");
    insert_if(l.link_libraries, "$(shell " + tfel_config +
                                    " --library-dependency "
                                    "--material --mfront-profiling)");
    for (const auto& h : ModellingHypothesis::getModellingHypotheses()) {
      insert_if(l.epts, name + "_" + ModellingHypothesis::toString(h));
    }
  }  // end of getTargetsDescription

  GenericParallelModelInterface::~GenericParallelModelInterface() = default;

}  // end of namespace mfront
