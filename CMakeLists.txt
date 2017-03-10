cmake_minimum_required (VERSION 2.8)

project (lady Fortran)

set (CMAKE_Fortran_FLAGS "-fdefault-real-8")

set (src
  "${PROJECT_SOURCE_DIR}/src/commons_mod.F90"
  "${PROJECT_SOURCE_DIR}/src/cartesian_domain_mod.F90"
  "${PROJECT_SOURCE_DIR}/src/cartesian_mesh_mod.F90"
  "${PROJECT_SOURCE_DIR}/src/shape_function_mod.F90"
  "${PROJECT_SOURCE_DIR}/src/quad_points_mod.F90"
  "${PROJECT_SOURCE_DIR}/src/parcels_mod.F90"
  "${PROJECT_SOURCE_DIR}/src/lady_mod.F90"
)

add_library (lady ${src})

add_executable (shallow_water_driver ${PROJECT_SOURCE_DIR}/src/shallow_water_driver.F90)
target_link_libraries (shallow_water_driver lady)
