add_executable(assemblage assembly.cxx)
target_link_libraries(assemblage PRIVATE mfem-mgis::MFEMMGIS)
if(MFEM_MGIS_HAS_GLOBAL_STIFFNESS_MATRIX)
  target_compile_definitions(assemblage
    PRIVATE MFEM_MGIS_HAS_GLOBAL_STIFFNESS_MATRIX)
endif()
add_mfront_library(assemblage-behaviours
  "${MFEM_MGIS_SOURCE_DIR}/tests/IsotropicLinearHardeningPlasticity.mfront"
  "${MFEM_MGIS_SOURCE_DIR}/tests/SaintVenantKirchhoffElasticity.mfront"
  "${MFEM_MGIS_SOURCE_DIR}/tests/StationaryNonLinearHeatTransfer.mfront")
foreach(_mesh beam-tet.mesh cube.mesh)
  configure_file("${MFEM_MGIS_SOURCE_DIR}/tests/${_mesh}"
    "${CMAKE_CURRENT_BINARY_DIR}/${_mesh}" COPYONLY)
endforeach()

# short tests: the matrices assembled with and without the positions of
# the entries are the same, and no entry is searched with the positions, in
# mechanics (small and finite strain) and in heat transfer, for tetrahedra
# of order 1 and hexahedra of orders 1 and 2
add_test(NAME assemblage-mecanique
  COMMAND ./assemblage -r 1 -n 2 -t mechanics)
add_test(NAME assemblage-thermique
  COMMAND ./assemblage -r 1 -n 2 -t thermal)
add_test(NAME assemblage-grandes-transformations
  COMMAND ./assemblage -m cube.mesh -r 2 -n 2 -t finite-strain)
add_test(NAME assemblage-hexaedres-ordre-2
  COMMAND ./assemblage -m cube.mesh -r 1 -o 2 -n 2 -t mechanics)
# in parallel, with the construction of the global matrix
add_test(NAME assemblage-parallele
  COMMAND mpirun -n 2 ./assemblage -r 1 -n 2 -t mechanics -p)
