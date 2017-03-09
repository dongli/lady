cmake_minimum_required (VERSION 2.8)

project (lady Fortran)

set (CMAKE_Fortran_FLAGS "-fdefault-real-8")

set (src
  "${PROJECT_SOURCE_DIR}/src/commons.F90"
  "${PROJECT_SOURCE_DIR}/src/cartesian_domain.F90"
  "${PROJECT_SOURCE_DIR}/src/cartesian_mesh.F90"
	# "${PROJECT_SOURCE_DIR}/src/kinetic_energy.F90"
  # "${PROJECT_SOURCE_DIR}/src/internal_energy.F90"
  # "${PROJECT_SOURCE_DIR}/src/potential_energy.F90"
  "${PROJECT_SOURCE_DIR}/src/shape_function.F90"
  "${PROJECT_SOURCE_DIR}/src/quad_points.F90"
  "${PROJECT_SOURCE_DIR}/src/parcels.F90"
  "${PROJECT_SOURCE_DIR}/src/lady.F90"
)

add_library (lady ${src})

add_executable (shallow_water_driver ${PROJECT_SOURCE_DIR}/src/shallow_water_driver.F90)
target_link_libraries (shallow_water_driver lady)
