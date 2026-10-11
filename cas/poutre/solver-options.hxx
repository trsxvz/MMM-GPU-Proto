/*!
 * \file   cas/poutre/solver-options.hxx
 * \brief  Options of the multigrids (BoomerAMG, AmgX) given on the command
 * line of the beam cases
 */

#ifndef LIB_CAS_POUTRE_SOLVER_OPTIONS_HXX
#define LIB_CAS_POUTRE_SOLVER_OPTIONS_HXX

#include <string>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include <string_view>
#include "MFEMMGIS/Parameters.hxx"

/*!
 * \return the options given as `Name=value` separated by commas: integers,
 * reals or strings
 * \param[in] o: options
 * \param[in] what: name of the options, for the error messages
 */
inline mfem_mgis::Parameters parseNameValueOptions(const std::string& o,
                                                   const std::string& what) {
  auto params = mfem_mgis::Parameters{};
  auto start = std::string::size_type{};
  while (start < o.size()) {
    const auto end = std::min(o.find(',', start), o.size());
    const auto option = o.substr(start, end - start);
    const auto equal = option.find('=');
    if ((equal == std::string::npos) || (equal == 0) ||
        (equal + 1 == option.size())) {
      std::cerr << "invalid option of " << what << " '" << option << "'\n";
      std::exit(EXIT_FAILURE);
    }
    const auto name = option.substr(0, equal);
    const auto value = option.substr(equal + 1);
    // integer, real or string, the whole value being read
    auto n = std::size_t{};
    try {
      const auto i = std::stoi(value, &n);
      if (n == value.size()) {
        params.replaceOrInsert(name, i);
        start = end + 1;
        continue;
      }
      const auto r = std::stod(value, &n);
      if (n == value.size()) {
        params.replaceOrInsert(name, r);
        start = end + 1;
        continue;
      }
    } catch (std::exception&) {
    }
    params.replaceOrInsert(name, value);
    start = end + 1;
  }
  return params;
}  // end of parseNameValueOptions

/*!
 * \return the options of BoomerAMG given as `Name=value` separated by
 * commas (see `parseNameValueOptions`), the verbosity level being null
 * \param[in] o: options
 */
inline mfem_mgis::Parameters parseBoomerAMGOptions(const std::string& o) {
  auto params = parseNameValueOptions(o, "BoomerAMG");
  if (!mfem_mgis::contains(params, "VerbosityLevel")) {
    params.replaceOrInsert("VerbosityLevel", 0);
  }
  return params;
}  // end of parseBoomerAMGOptions

/*!
 * \return the options of the solvers of AmgX: a configuration file, or
 * parameters of AmgX completing the default configuration of mfem-mgis
 * (`AmgXParameters`), given as `name=value` separated by commas, and the
 * precision
 * \param[in] file: configuration file, if not empty
 * \param[in] o: parameters of AmgX, if not empty
 * \param[in] single: computations of AmgX in single precision
 */
inline mfem_mgis::Parameters parseAmgXOptions(std::string_view file,
                                              const std::string& o,
                                              const bool single) {
  if ((!file.empty()) && (!o.empty())) {
    std::cerr << "a configuration file of AmgX and parameters of AmgX can't "
                 "be given together\n";
    std::exit(EXIT_FAILURE);
  }
  auto params = mfem_mgis::Parameters{};
  if (!file.empty()) {
    params.replaceOrInsert("ConfigurationFile", std::string{file});
  }
  if (!o.empty()) {
    params.replaceOrInsert("AmgXParameters", parseNameValueOptions(o, "AmgX"));
  }
  if (single) {
    params.replaceOrInsert("Precision", std::string{"Single"});
  }
  return params;
}  // end of parseAmgXOptions

/*!
 * \return the options of the conjugate gradient of AmgX (`AmgXPCG`)
 * \param[in] amgx: options of AmgX (see `parseAmgXOptions`)
 * \param[in] tolerance: relative tolerance
 * \param[in] n: maximum number of iterations
 */
inline mfem_mgis::Parameters amgxPCGOptions(const mfem_mgis::Parameters& amgx,
                                            const mfem_mgis::real tolerance,
                                            const int n) {
  auto params = amgx;
  params.replaceOrInsert("RelativeTolerance", tolerance);
  params.replaceOrInsert("MaximumNumberOfIterations", n);
  return params;
}  // end of amgxPCGOptions

#endif /* LIB_CAS_POUTRE_SOLVER_OPTIONS_HXX */
