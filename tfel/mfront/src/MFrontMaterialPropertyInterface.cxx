/*!
 * \file   mfront/src/MFrontMaterialPropertyInterface.cxx
 * \brief
 * \author Thomas Helfer
 * \date   06 mai 2008
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <sstream>
#include <stdexcept>
#include <algorithm>
#include "TFEL/Raise.hxx"
#include "TFEL/Config/GetInstallPath.hxx"
#include "MFront/CodeGeneratorUtilities.hxx"
#include "MFront/MFrontUtilities.hxx"
#include "MFront/MFrontHeader.hxx"
#include "MFront/FileDescription.hxx"
#include "MFront/TargetsDescription.hxx"
#include "MFront/MaterialPropertyDescription.hxx"
#include "MFront/MFrontMaterialPropertyInterface.hxx"

namespace mfront {

  std::string MFrontMaterialPropertyInterface::getName() { return "mfront"; }

  MFrontMaterialPropertyInterface::MFrontMaterialPropertyInterface() = default;

  std::vector<std::string>
  MFrontMaterialPropertyInterface::getSupportedFloatingPointTypes() const {
    return {"float", "double", "long double"};
  }  // end of getSupportedFloatingPointTypes

  bool MFrontMaterialPropertyInterface::
      shallGenerateOverloadedFunctionForQuantities() const {
    return true;
  }  // end of shallGenerateOverloadedFunctionForQuantities

  std::pair<bool,
            tfel::utilities::CxxTokenizer::TokensContainer::const_iterator>
  MFrontMaterialPropertyInterface::treatKeyword(
      const std::string& k,
      const std::vector<std::string>& i,
      tokens_iterator current,
      const tokens_iterator) {
    tfel::raise_if(std::find(i.begin(), i.end(), "mfront") != i.end(),
                   "MFrontMaterialPropertyInterface::treatKeyword: "
                   "unsupported key '" +
                       k + "'");
    return {false, current};
  }  // end of treatKeyword

  void MFrontMaterialPropertyInterface::getTargetsDescription(
      TargetsDescription& d, const MaterialPropertyDescription& mpd) const {
    const auto lib = "MFrontMaterialLaw";
    const auto name = this->getSrcFileName(mpd.material, mpd.className);
    const auto f = mpd.material.empty() ? mpd.className
                                        : mpd.material + "_" + mpd.className;
    const auto header = this->getHeaderFileName(mpd.material, mpd.className);
    const auto tfel_config = tfel::getTFELConfigExecutableName();
    auto& l = d.getLibrary(lib);
    if (!header.empty()) {
      insert_if(d.headers, header + ".hxx");
    }
    insert_if(l.sources, name + ".cxx");
    insert_if(l.link_libraries, "m");
    insert_if(l.epts, {f, f + "_checkBounds"});
  }  // end of getTargetsDescription

  std::string MFrontMaterialPropertyInterface::getHeaderFileName(
      const std::string& material, const std::string& className) const {
    if (material.empty()) {
      return className + "-mfront";
    }
    return material + "_" + className + "-mfront";
  }  // end of getHeaderFileName

  std::string MFrontMaterialPropertyInterface::getSrcFileName(
      const std::string& material, const std::string& className) const {
    if (material.empty()) {
      return className + "-mfront";
    }
    return material + "_" + className + "-mfront";
  }  // end of getSrcFileName

  void MFrontMaterialPropertyInterface::writeBeginHeaderNamespace(
      std::ostream& os) const {
    os << "namespace mfront{\n\n";
  }  // end of writeBeginHeaderNamespace

  void MFrontMaterialPropertyInterface::writeEndHeaderNamespace(
      std::ostream& os) const {
    os << "} // end of namespace mfront\n\n";
  }  // end of writeEndHeaderNamespace()

  void MFrontMaterialPropertyInterface::writeBeginSrcNamespace(
      std::ostream& os) const {
    os << "namespace mfront{\n\n";
  }  // end of writeBeginSrcNamespace

  void MFrontMaterialPropertyInterface::writeEndSrcNamespace(
      std::ostream& os) const {
    os << "} // end of namespace mfront\n";
  }  // end of writeEndSrcNamespace()

  std::string MFrontMaterialPropertyInterface::getFunctionName(
      const MaterialPropertyDescription& mpd) const {
    const auto material = mpd.material;
    const auto law = mpd.className;
    return material.empty() ? law : material + '_' + law;
  }  // end of getFunctionName

  bool MFrontMaterialPropertyInterface::requiresCheckBoundsFunction() const {
    return true;
  }

  std::string MFrontMaterialPropertyInterface::getCheckBoundsFunctionName(
      const MaterialPropertyDescription& mpd) const {
    return this->getFunctionName(mpd) + "_checkBounds";
  }  // end of getCheckBoundsFunctionName

  void MFrontMaterialPropertyInterface::writeFileDescriptionSymbols(
      std::ostream&,
      const MaterialPropertyDescription&,
      const FileDescription&) const {}  // end of writeFileDescriptionSymbols

  void MFrontMaterialPropertyInterface::writeValidatorSymbol(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeValidatorSymbol

  void MFrontMaterialPropertyInterface::writeBuildIdentifierSymbol(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeBuildIdentifierSymbol

  void MFrontMaterialPropertyInterface::writeEntryPointSymbol(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeEntryPointSymbol

  void MFrontMaterialPropertyInterface::writeTFELVersionSymbol(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeTFELVersionSymbol

  void MFrontMaterialPropertyInterface::writeUnitSystemSymbol(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeUnitSystemSymbol

  void MFrontMaterialPropertyInterface::writeLawSymbol(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeMaterialSymbol

  void MFrontMaterialPropertyInterface::writeMaterialSymbol(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeMaterialSymbol

  void MFrontMaterialPropertyInterface::writeInterfaceSymbol(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeInterfaceSymbol

  void MFrontMaterialPropertyInterface::writeMaterialKnowledgeTypeSymbol(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeEntryPointSymbol

  void MFrontMaterialPropertyInterface::writeVariablesNamesSymbol(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeVariablesNamesSymbol

  void MFrontMaterialPropertyInterface::writeVariablesBoundsSymbols(
      std::ostream&, const MaterialPropertyDescription&) const {
  }  // end of writeVariablesBoundsSymbols

  void MFrontMaterialPropertyInterface::writeDeviceImplementations(
      std::ostream& os,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd) const {
    const auto name = this->getFunctionName(mpd);
    os << "#include <cmath>\n"
       << "#include <limits>\n"
       << "#include \"TFEL/Config/TFELConfig.hxx\"\n"
       << "#include \"TFEL/Config/TFELTypes.hxx\"\n"
       << "#include \"TFEL/PhysicalConstants.hxx\"\n"
       << "#include \"TFEL/Math/General/IEEE754.hxx\"\n"
       << "#include \"TFEL/Math/General/DerivativeType.hxx\"\n";
    if (useQuantities(mpd)) {
      os << "#include \"TFEL/Math/qt.hxx\"\n";
    }
    if (!mpd.includes.empty()) {
      os << mpd.includes << "\n";
    }
    os << "\nnamespace mfront::device {\n\n"
       << "/*!\n"
       << " * \\brief implementation of the " << name
       << " material property which can be\n"
       << " * called in device code. This implementation does not use\n"
       << " * exceptions nor errno: an invalid value (NaN) is returned if\n"
       << " * a physical bound is not satisfied.\n"
       << " */\n";
    if (mpd.inputs.empty()) {
      os << "template <typename mfront_ValueType = double>\n"
         << "TFEL_HOST_DEVICE inline mfront_ValueType " << name << "(){\n";
      this->writeDeviceMaterialPropertyBody(os, mpd, fd, "mfront_ValueType");
    } else {
      os << "TFEL_HOST_DEVICE inline double " << name << "(";
      this->writeArgumentsList(os, mpd, "double", false);
      os << "){\n";
      this->writeDeviceMaterialPropertyBody(os, mpd, fd, "double");
    }
    os << "} // end of " << name << "\n\n"
       << "} // end of namespace mfront::device\n\n";
  }  // end of writeDeviceImplementations

  void MFrontMaterialPropertyInterface::writeDeviceMaterialPropertyBody(
      std::ostream& os,
      const MaterialPropertyDescription& mpd,
      const FileDescription& fd,
      const std::string_view t) const {
    const auto nan = "std::numeric_limits<" + std::string{t} + ">::quiet_NaN()";
    auto write_physical_bounds_check = [&os, &nan](const VariableDescription& v,
                                                   const std::string& n) {
      if (!v.hasPhysicalBounds()) {
        return;
      }
      const auto& b = v.getPhysicalBounds();
      if (b.boundsType == VariableBoundsDescription::LOWER) {
        os << "if(" << n << " < " << b.lowerBound << "){\n";
      } else if (b.boundsType == VariableBoundsDescription::UPPER) {
        os << "if(" << n << " > " << b.upperBound << "){\n";
      } else {
        os << "if((" << n << " < " << b.lowerBound << ") || (" << n << " > "
           << b.upperBound << ")){\n";
      }
      os << "return " << nan << ";\n"
         << "}\n";
    };
    writeBeginningOfMaterialPropertyBody(os, mpd, fd, t, true);
    for (const auto& p : mpd.parameters) {
      tfel::raise_if(!p.hasAttribute(VariableDescription::defaultValue),
                     "MFrontMaterialPropertyInterface::"
                     "writeDeviceMaterialPropertyBody: "
                     "internal error (can't find value of parameter '" +
                         p.name + "')");
      os << "constexpr auto " << p.name << " = " << p.type << "("
         << p.getAttribute<double>(VariableDescription::defaultValue)
         << ");\n";
    }
    this->writeInterfaceSpecificVariables(os, mpd);
    const auto use_qt = useQuantities(mpd);
    if (!areRuntimeChecksDisabled(mpd)) {
      for (const auto& i : mpd.inputs) {
        write_physical_bounds_check(i, (use_qt ? "mfront_" : "") + i.name);
      }
    }
    os << "auto " << mpd.output.name << " = " << mpd.output.type << "{};\n"
       << mpd.f.body << "\n";
    const auto output =
        use_qt ? mpd.output.name + ".getValue()" : mpd.output.name;
    if (!areRuntimeChecksDisabled(mpd)) {
      write_physical_bounds_check(mpd.output, output);
    }
    os << "return " << output << ";\n";
  }  // end of writeDeviceMaterialPropertyBody

  MFrontMaterialPropertyInterface::~MFrontMaterialPropertyInterface() = default;

}  // end of namespace mfront
