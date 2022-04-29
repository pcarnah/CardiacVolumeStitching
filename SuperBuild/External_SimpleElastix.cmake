set(proj SimpleElastix)
set(ep_common_cxx_flags "${CMAKE_CXX_FLAGS_INIT} ${ADDITIONAL_CXX_FLAGS}")
# Set dependency list
set(${proj}_DEPENDS elastix)

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj)

if(NOT DEFINED ${proj}_DIR AND NOT ${CMAKE_PROJECT_NAME}_USE_SYSTEM_${proj})


  set(git_protocol "https")

  set(${proj}_INSTALL_DIR ${CMAKE_BINARY_DIR}/${proj}-install)
  set(${proj}_DIR ${CMAKE_BINARY_DIR}/${proj}-build)
  set(${proj}_RUNTIME_DIR ${CMAKE_BINARY_DIR}/${Slicer_THIRDPARTY_BIN_DIR})
  
  set(${proj}_cxx_flags "${ep_common_cxx_flags}")
  
  if (APPLE)
    set(${proj}_cxx_flags "${ep_common_cxx_flags} -Wno-inconsistent-missing-override")
  endif()

  ExternalProject_Add(${proj}
    # Slicer
    ${${proj}_EP_ARGS}
    SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj}
    #SOURCE_SUBDIR src # requires CMake 3.7 or later
    BINARY_DIR ${proj}-build
    INSTALL_DIR ${${proj}_INSTALL_DIR}
    GIT_REPOSITORY "${git_protocol}://github.com/pcarnah/SimpleElastix.git"
    #GIT_TAG "4e3bf31b7997d26bbd904c18a9882188d4368790"
    #--Patch step-------------
    PATCH_COMMAND ""
    #--Configure step-------------
    CMAKE_CACHE_ARGS
      -DSubversion_SVN_EXECUTABLE:STRING=${Subversion_SVN_EXECUTABLE}
      -DGIT_EXECUTABLE:STRING=${GIT_EXECUTABLE}
	  -DSimpleITK_DIR:STRING=${SimpleITK_DIR}
      -DITK_DIR:STRING=${ITK_DIR}
	  -DElastix_DIR:STRING=${elastix_DIR}
      -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER}
      -DCMAKE_CXX_FLAGS:STRING=${${proj}_cxx_flags}
      -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER}
      -DCMAKE_C_FLAGS:STRING=${ep_common_c_flags}
      -DCMAKE_CXX_STANDARD:STRING=${CMAKE_CXX_STANDARD}
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      -DSWIG_EXECUTABLE:PATH=${SWIG_EXECUTABLE}
      -DPYTHON_EXECUTABLE:PATH=${PYTHON_EXECUTABLE}
      -DPYTHON_LIBRARY:FILEPATH=${PYTHON_LIBRARY}
      -DPYTHON_INCLUDE_DIR:PATH=${PYTHON_INCLUDE_DIR}
	  -DSimpleITK_PYTHON_USE_VIRTUALENV:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DBUILD_DOXYGEN:BOOL=OFF
	  -DBUILD_EXAMPLES:BOOL=OFF
      -DWRAP_CSHARP:BOOL=OFF 
	  -DWRAP_DEFAULT:BOOL=OFF 
	  -DWRAP_JAVA:BOOL=OFF 
	  -DSimpleITK_INT64_PIXELIDS:BOOL=OFF 
	  -DWRAP_PYTHON:BOOL=ON 
	  -DBUILD_SHARED_LIBS:BOOL=ON
      -DCMAKE_RUNTIME_OUTPUT_DIRECTORY:PATH=${CMAKE_BINARY_DIR}/${Slicer_THIRDPARTY_BIN_DIR}
      -DCMAKE_LIBRARY_OUTPUT_DIRECTORY:PATH=${CMAKE_BINARY_DIR}/${Slicer_THIRDPARTY_LIB_DIR}
      -DCMAKE_ARCHIVE_OUTPUT_DIRECTORY:PATH=${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}
    #--Build step-----------------
    #--Install step-----------------
    INSTALL_COMMAND ""
    DEPENDS
      ${${proj}_DEPENDS}
    )
  #set(${proj}_DIR ${${proj}_INSTALL_DIR})
  #if(UNIX)
  #  set(${proj}_DIR ${${proj}_INSTALL_DIR}/share/elastix)
  #endif()
  
else()
  ExternalProject_Add_Empty(${proj} DEPENDS ${${proj}_DEPENDS})
endif()

mark_as_superbuild(${proj}_DIR:PATH)
