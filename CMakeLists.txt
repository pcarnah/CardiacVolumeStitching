cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(CardiacVolumeStitching)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "http://slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/CardiacVolumeStitching")
set(EXTENSION_CATEGORY "Cardiac")
set(EXTENSION_CONTRIBUTORS "Patrick Carnahan (Robarts Research Institute)")
set(EXTENSION_DESCRIPTION "This is an example of a simple extension")
set(EXTENSION_ICONURL "http://www.example.com/Slicer/Extensions/CardiacVolumeStitching.png")
set(EXTENSION_SCREENSHOTURLS "http://www.example.com/Slicer/Extensions/CardiacVolumeStitching/Screenshots/1.png")
set(EXTENSION_DEPENDS "NA") # Specified as a space separated string, a list or 'NA' if any
set(EXTENSION_BUILD_SUBDIRECTORY inner-build)

set(SUPERBUILD_TOPLEVEL_PROJECT inner)

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})
mark_as_superbuild(Slicer_DIR)

find_package(Git REQUIRED)
mark_as_superbuild(GIT_EXECUTABLE)

#-----------------------------------------------------------------------------
# SuperBuild setup
option(${EXTENSION_NAME}_SUPERBUILD "Build ${EXTENSION_NAME} and the projects it depends on." ON)
mark_as_advanced(${EXTENSION_NAME}_SUPERBUILD)
if(${EXTENSION_NAME}_SUPERBUILD)
  include("${CMAKE_CURRENT_SOURCE_DIR}/SuperBuild.cmake")
  return()
endif()

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(CardiacVolumeStitching)
add_subdirectory(MonogenicFilter)
## NEXT_MODULE

#-----------------------------------------------------------------------------
# Install SimpleElastix
# 
find_package(SimpleElastix)
install(FILES 
	${SimpleElastix_DIR}/Wrapping/Python/SimpleElastix.py
	${SimpleElastix_DIR}/Wrapping/Python/_SimpleElastix.pyd
	${SimpleElastix_DIR}/Wrapping/Python/__init__.py
DESTINATION ${Slicer_THIRDPARTY_LIB_DIR}/qt-scripted-modules/SimpleElastix/
COMPONENT RuntimeLibraries
)

install(DIRECTORY
	${SimpleElastix_RUNTIME_DIRS}/$<CONFIG>/
DESTINATION ${Slicer_THIRDPARTY_LIB_DIR}/qt-scripted-modules/SimpleElastix/
COMPONENT RuntimeLibraries
)

#-----------------------------------------------------------------------------
# Install directory, install project name, install component, and install subdirectory.
set(EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS)
#list(APPEND EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS "${elastix_DIR};elastix;RuntimeLibraries;/")
set(${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS "${EXTENSION_CPACK_INSTALL_CMAKE_PROJECTS}" CACHE STRING "List of external projects to install" FORCE)

#-----------------------------------------------------------------------------
list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${CMAKE_BINARY_DIR};${EXTENSION_NAME};ALL;/")
list(APPEND CPACK_INSTALL_CMAKE_PROJECTS "${${EXTENSION_NAME}_CPACK_INSTALL_CMAKE_PROJECTS}")
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
