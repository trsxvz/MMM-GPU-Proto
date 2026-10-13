/*!
 * \file   tests/MeshReadModeTest.cxx
 * \brief  This file tests the errors of the mesh read modes: they are
 * reported to the calling code, which may handle them, instead of aborting
 * the program
 * \date   13/10/2026
 */

#include <string>
#include <cstdlib>
#include <iostream>
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"

/*!
 * \return true if building a finite element discretization with the given
 * mesh read mode fails with an error message containing the given text
 * \param[in] mesh_file: mesh file
 * \param[in] mode: mesh read mode
 * \param[in] parallel: parallel computation
 * \param[in] expected: text expected in the error message
 */
static bool checkRefusal(const std::string& mesh_file,
                         const std::string& mode,
                         const bool parallel,
                         const std::string& expected) {
  auto ctx = mfem_mgis::Context{};
  const auto fed =
      mfem_mgis::make_shared<mfem_mgis::FiniteElementDiscretization>(
          ctx, mfem_mgis::Parameters{{"MeshFileName", mesh_file},
                                     {"MeshReadMode", mode},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", 1},
                                     {"UnknownsSize", 3},
                                     {"Parallel", parallel}});
  if (fed) {
    std::cerr << "mode '" << mode << "' accepted\n";
    return false;
  }
  const auto message = ctx.getErrorMessage();
  if (message.find(expected) == std::string::npos) {
    std::cerr << "mode '" << mode << "': unexpected message '" << message
              << "'\n";
    return false;
  }
  return true;
}  // end of checkRefusal

int main(int argc, char** argv) {
  mfem_mgis::initialize(argc, argv);
  if (argc != 2) {
    std::cerr << "usage: MeshReadModeTest mesh_file\n";
    return EXIT_FAILURE;
  }
  const auto mesh = std::string{argv[1]};
  auto success = true;
  // restart from checkpoint files which do not exist
  success = checkRefusal("does-not-exist", "Restart", true,
                         "checkpoint file 'does-not-exist") &&
            success;
  // restart in a sequential computation
  success = checkRefusal(mesh, "Restart", false, "'Restart' is not handled") &&
            success;
  // invalid mode
  success = checkRefusal(mesh, "FromNowhere", true, "MeshReadMode") && success;
  std::cout << (success ? "success" : "failure") << '\n';
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
