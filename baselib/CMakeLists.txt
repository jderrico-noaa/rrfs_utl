set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include")

add_library(kinds STATIC kinds.f90)
add_library(${PROJECT_NAME}::kinds ALIAS kinds)
set_target_properties(kinds PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
target_include_directories(kinds PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                       $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
list(APPEND BASE_LIB_Targets kinds)

add_library(ncio STATIC module_ncio.f90)
add_library(${PROJECT_NAME}::ncio ALIAS ncio)
set_target_properties(ncio PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
target_link_libraries(ncio PUBLIC NetCDF::NetCDF_Fortran)
target_include_directories(ncio PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                       $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
list(APPEND BASE_LIB_Targets ncio)

add_library (map_util STATIC module_map_utils.f90)
add_library(${PROJECT_NAME}::map_util ALIAS map_util)
target_link_libraries(map_util PRIVATE ${PROJECT_NAME}::kinds)
set_target_properties(map_util PROPERTIES Fortran_MODULE_DIRECTORY "${module_dir}")
target_include_directories(map_util PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
                                           $<INSTALL_INTERFACE:include>)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include DESTINATION ${CMAKE_INSTALL_PREFIX})
list(APPEND BASE_LIB_Targets map_util)

install(
  TARGETS ${BASE_LIB_Targets}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
