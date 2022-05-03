#-----------------------------------------------------------------------------
set(MODULE_NAME CardiacVolumeStitching)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/Icons/${MODULE_NAME}.png
  Resources/RegistrationParameters/Par0001Rigid.txt
  Resources/RegistrationParameters/Par0001NonRigid.txt
  Resources/RegistrationParameters/Par0001RigidTest.txt
  Resources/RegistrationParameters/Par0003RigidGroup.txt
  Resources/RegistrationParameters/Par0003NonRigidGroup.txt
  )

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)

  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${MODULE_NAME}.py)

  # Additional build-time testing
  add_subdirectory(Testing)
endif()
