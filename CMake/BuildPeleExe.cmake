function(build_pele_exe pele_exe_name pele_physics_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)

  add_executable(${pele_exe_name} "")

  if(CLANG_TIDY_EXE)
    set_target_properties(${pele_exe_name} PROPERTIES CXX_CLANG_TIDY
                          "${CLANG_TIDY_EXE};--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy")
  endif()

  target_sources(${pele_exe_name}
     PRIVATE
       pelelmex_prob_parm.H
       pelelmex_prob.H
       pelelmex_prob.cpp
  )

  target_include_directories(${pele_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${pele_exe_name} PRIVATE ${SRC_DIR})
  target_include_directories(${pele_exe_name} PRIVATE ${CMAKE_BINARY_DIR})

  target_sources(${pele_exe_name}
     PRIVATE
       ${SRC_DIR}/PeleLMeX_DeriveUserDefined.cpp
       ${SRC_DIR}/PeleLMeX_DiffusionOp.H
       ${SRC_DIR}/PeleLMeX_DiffusionOp.cpp
       ${SRC_DIR}/PeleLMeX_EBUserDefined.H
       ${SRC_DIR}/PeleLMeX_FlowControllerData.H
       ${SRC_DIR}/PeleLMeX.H
       ${SRC_DIR}/PeleLMeX.cpp
       ${SRC_DIR}/PeleLMeX_Advance.cpp
       ${SRC_DIR}/PeleLMeX_Advection.cpp
       ${SRC_DIR}/PeleLMeX_BC.cpp
       ${SRC_DIR}/PeleLMeX_BCfill.H
       ${SRC_DIR}/PeleLMeX_Data.cpp
       ${SRC_DIR}/PeleLMeX_Derive.H
       ${SRC_DIR}/PeleLMeX_Derive.cpp
       ${SRC_DIR}/PeleLMeX_DeriveFunc.H
       ${SRC_DIR}/PeleLMeX_DeriveFunc.cpp
       ${SRC_DIR}/PeleLMeX_Diagnostics.cpp
       ${SRC_DIR}/PeleLMeX_Diffusion.cpp
       ${SRC_DIR}/PeleLMeX_EB.cpp
       ${SRC_DIR}/PeleLMeX_Eos.cpp
       ${SRC_DIR}/PeleLMeX_Evaluate.cpp
       ${SRC_DIR}/PeleLMeX_Evolve.cpp
       ${SRC_DIR}/PeleLMeX_FlowController.cpp
       ${SRC_DIR}/PeleLMeX_Forces.cpp
       ${SRC_DIR}/PeleLMeX_PatchFlowVariables.H
       ${SRC_DIR}/PeleLMeX_PatchFlowVariables.cpp
       ${SRC_DIR}/PeleLMeX_Init.cpp
       ${SRC_DIR}/PeleLMeX_Plot.cpp
       ${SRC_DIR}/PeleLMeX_ProblemSpecificFunctions.H
       ${SRC_DIR}/PeleLMeX_ProblemSpecificFunctions.cpp
       ${SRC_DIR}/PeleLMeX_Projection.cpp
       ${SRC_DIR}/PeleLMeX_Reactions.cpp
       ${SRC_DIR}/PeleLMeX_Regrid.cpp
       ${SRC_DIR}/PeleLMeX_Setup.cpp
       ${SRC_DIR}/PeleLMeX_Tagging.cpp
       ${SRC_DIR}/PeleLMeX_BPatch.H
       ${SRC_DIR}/PeleLMeX_BPatch.cpp
       ${SRC_DIR}/PeleLMeX_Temporals.cpp
       ${SRC_DIR}/PeleLMeX_Timestep.cpp
       ${SRC_DIR}/PeleLMeX_TransportProp.cpp
       ${SRC_DIR}/PeleLMeX_UMac.cpp
       ${SRC_DIR}/PeleLMeX_UserKeys.H
       ${SRC_DIR}/PeleLMeX_Utils.H
       ${SRC_DIR}/PeleLMeX_Utils.cpp
       ${SRC_DIR}/PeleLMeX_Index.H
       ${SRC_DIR}/PeleLMeX_K.H
       ${SRC_DIR}/main.cpp
  )

  if(PELE_PHYSICS_ENABLE_SOOT)
    target_sources(${pele_exe_name}
      PRIVATE
        ${SRC_DIR}/PeleLMeX_Soot.cpp)
  endif()

  if(PELE_PHYSICS_ENABLE_SPRAY)
    target_sources(${pele_exe_name}
      PRIVATE
        ${SRC_DIR}/PeleLMeX_SprayParticles.cpp
        SprayParticlesInitInsert.cpp)
  endif()

  if(PELE_PHYSICS_ENABLE_RADIATION)
    target_sources(${pele_exe_name}
      PRIVATE
        ${SRC_DIR}/PeleLMeX_Radiation.cpp)
  endif()

  if(NOT "${pele_exe_name}" STREQUAL "${PROJECT_NAME}-UnitTests")
    target_sources(${pele_exe_name}
       PRIVATE
         ${CMAKE_SOURCE_DIR}/Source/main.cpp
    )
  endif()

  if(PELE_ENABLE_CUDA)
    set(pctargets "${pele_exe_name};${pele_physics_lib_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(PELE_SOURCES ${tgt} SOURCES)
      list(FILTER PELE_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${PELE_SOURCES} PROPERTIES LANGUAGE CUDA)
    endforeach()
    set_target_properties(${pele_exe_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(${pele_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas --disable-optimizer-constants>)
  endif()

  target_link_libraries(${pele_exe_name} PRIVATE ${pele_physics_lib_name} AMReX-Hydro::amrex_hydro_api AMReX::amrex)

  install(TARGETS ${pele_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
