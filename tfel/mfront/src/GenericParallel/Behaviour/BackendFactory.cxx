/*!
 * \file   mfront/src/GenericParallel/Behaviour/BackendFactory.cxx
 * \brief
 * \date   05/10/2026
 * \copyright Copyright (C) 2006-2025 CEA/DEN, EDF R&D. All rights
 * reserved.
 * This project is publicly released under either the GNU GPL Licence with
 * linking exception or the CECILL-A licence. A copy of thoses licences are
 * delivered with the sources of TFEL. CEA or EDF may also distribute this
 * project under specific licensing conditions.
 */

#include <utility>
#include "TFEL/Raise.hxx"
#include "MFront/MFrontLogStream.hxx"
#include "MFront/GenericParallel/Behaviour/GPUBackend.hxx"
#include "MFront/GenericParallel/Behaviour/KokkosBackend.hxx"
#include "MFront/GenericParallel/Behaviour/OpenMPBackend.hxx"
#include "MFront/GenericParallel/Behaviour/ParallelSTLBackend.hxx"
#include "MFront/GenericParallel/Behaviour/SYCLBackend.hxx"
#include "MFront/GenericParallel/Behaviour/BackendFactory.hxx"

namespace mfront::generic_parallel::behaviour {

  BackendFactory& BackendFactory::get() noexcept {
    static BackendFactory f;
    return f;
  }  // end of get

  BackendFactory::BackendFactory() {
    this->add("cuda", [](const DataMap& opts) {
      return std::make_unique<CUDABackend>(opts);
    });
    this->add("hip", [](const DataMap& opts) {
      return std::make_unique<HIPBackend>(opts);
    });
    this->add("kokkos", [](const DataMap& opts) {
      return std::make_unique<KokkosBackend>(opts);
    });
    this->add("openmp", [](const DataMap& opts) {
      return std::make_unique<OpenMPBackend>(opts);
    });
    this->add("stlpar", [](const DataMap& opts) {
      return std::make_unique<ParallelSTLBackend>(opts);
    });
    this->add("sycl", [](const DataMap& opts) {
      return std::make_unique<SYCLBackend>(opts);
    });
  }  // end of BackendFactory

  void BackendFactory::add(std::string_view n, Generator g) {
    if (!this->generators.insert({std::string{n}, std::move(g)}).second) {
      tfel::raise("BackendFactory::add: a backend named '" + std::string{n} +
                  "' has already been registred");
    }
  }  // end of add

  std::unique_ptr<AbstractBackend> BackendFactory::generate(
      std::string_view n, const DataMap& opts) const {
    const auto p = this->generators.find(n);
    if (p == this->generators.end()) {
      auto names = std::string{};
      for (const auto& g : this->generators) {
        names += (names.empty() ? "" : ", ") + g.first;
      }
      tfel::raise("BackendFactory::generate: no backend named '" +
                  std::string{n} + "' (available backends: " + names + ")");
    }
    try {
      return p->second(opts);
    } catch (std::exception& e) {
      tfel::raise("BackendFactory::generate: error while generating backend '" +
                  std::string{n} + "'\n" + std::string{e.what()});
    }
  }  // end of generate

  std::unique_ptr<AbstractBackend> BackendFactory::generate(
      const DataMap& opts) const {
    using tfel::utilities::DataStructure;
    auto validator = tfel::utilities::DataMapValidator{};
    validator.addDataTypeValidator<DataStructure>("backend");
    validator.validate(opts);
    if (!contains(opts, "backend")) {
      return {};
    }
    const auto ds =
        DataStructure::convert(tfel::utilities::get(opts, "backend"));
    if (getVerboseMode() >= VERBOSE_DEBUG) {
      getLogStream() << "BackendFactory::generate: creating backend '"
                     << ds.name << "'\n";
    }
    return this->generate(ds.name, ds.data);
  }  // end of generate

  BackendFactory::~BackendFactory() = default;

}  // end of namespace mfront::generic_parallel::behaviour
