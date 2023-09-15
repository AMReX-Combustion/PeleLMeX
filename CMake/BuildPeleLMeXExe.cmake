function(build_pelelmex_exe pelelmex_exe_name pelelmex_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)

  add_executable(${pelelmex_exe_name} "")

  if(CLANG_TIDY_EXE)
    set_target_properties(${pelelmex_exe_name} PROPERTIES CXX_CLANG_TIDY 
                          "${CLANG_TIDY_EXE};--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy")
  endif()

  target_sources(${pelelmex_exe_name}
     PRIVATE
       pelelm_prob_parm.H
       pelelm_prob.H
       pelelm_prob.cpp
  )
  
  #PeleLMeX include directories
  target_include_directories(${pelelmex_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${pelelmex_exe_name} PRIVATE ${SRC_DIR})
  target_include_directories(${pelelmex_exe_name} PRIVATE ${CMAKE_BINARY_DIR})
  #target_include_directories(${pelelmex_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Source/Params/param_includes)

  # Set PeleMP flags
  set(PELEMP_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/PeleMP/Source)
  if(PELELMEX_ENABLE_AMREX_PARTICLES AND PELEMP_SPRAY_FUEL_NUM GREATER 0)
    target_compile_definitions(${pelelmex_exe_name} PRIVATE PELELMEX_USE_SPRAY)
    target_compile_definitions(${pelelmex_exe_name} PRIVATE SPRAY_FUEL_NUM=${PELEMP_SPRAY_FUEL_NUM})
    target_sources(${pelelmex_exe_name} PRIVATE
	           SprayParticlesInitInsert.cpp
                   ${SRC_DIR}/PeleLMSprayParticles.cpp
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayParticles.cpp
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayParticles.H
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayFuelData.H
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayInterpolation.H
                   ${PELEMP_SRC_DIR}/PP_Spray/Drag.H
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayInjection.H
                   ${PELEMP_SRC_DIR}/PP_Spray/SpraySetup.cpp
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayDerive.cpp
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayJet.H
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayJet.cpp
                   ${PELEMP_SRC_DIR}/PP_Spray/SprayIO.cpp
                   ${PELEMP_SRC_DIR}/PP_Spray/WallFunctions.H
                   ${PELEMP_SRC_DIR}/PP_Spray/Distribution/DistBase.H
                   ${PELEMP_SRC_DIR}/PP_Spray/Distribution/Distributions.H
                   ${PELEMP_SRC_DIR}/PP_Spray/Distribution/Distributions.cpp)
    target_include_directories(${pelelmex_exe_name} PRIVATE ${PELEMP_SRC_DIR}/PP_Spray)
    target_include_directories(${pelelmex_exe_name} PRIVATE ${PELEMP_SRC_DIR}/PP_Spray/Distribution)
  endif()
  if(PELELMEX_ENABLE_SOOT)
    target_compile_definitions(${pelelmex_exe_name} PRIVATE PELELMEX_USE_SOOT)
    target_compile_definitions(${pelelmex_exe_name} PRIVATE NUM_SOOT_MOMENTS=${PELEMP_NUM_SOOT_MOMENTS})
    set(SOOT_MOMENTS_VALUES 3 6)
    if(NOT PELEMP_NUM_SOOT_MOMENTS IN_LIST SOOT_MOMENTS_VALUES)
      message(FATAL_ERROR "NUM_SOOT_MOMENTS must be either 3 or 6")
    endif()
    target_sources(${pelelmex_exe_name} PRIVATE
                   ${SRC_DIR}/Soot.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel_react.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel_derive.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/Constants_Soot.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootData.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootReactions.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel_derive.H
                   ${SRC_DIR}/PeleLMSoot.cpp)
    target_include_directories(${pelelmex_exe_name} PRIVATE ${PELEMP_SRC_DIR}/Soot_Models)
  endif()

  target_sources(${pelelmex_exe_name}
     PRIVATE
       ${SRC_DIR}/DeriveUserDefined.cpp
       ${SRC_DIR}/DiffusionOp.H
       ${SRC_DIR}/DiffusionOp.cpp
       ${SRC_DIR}/EBUserDefined.H
       ${SRC_DIR}/PeleFlowControllerData.H
       ${SRC_DIR}/PeleLM.H
       ${SRC_DIR}/PeleLM.cpp
       ${SRC_DIR}/PeleLMAdvance.cpp
       ${SRC_DIR}/PeleLMAdvection.cpp
       ${SRC_DIR}/PeleLMBC.cpp
       ${SRC_DIR}/PeleLMBCfill.H
       ${SRC_DIR}/PeleLMData.cpp
       ${SRC_DIR}/PeleLMDerive.H
       ${SRC_DIR}/PeleLMDerive.cpp
       ${SRC_DIR}/PeleLMDeriveFunc.H
       ${SRC_DIR}/PeleLMDeriveFunc.cpp
       ${SRC_DIR}/PeleLMDiagnostics.cpp
       ${SRC_DIR}/PeleLMDiffusion.cpp
       ${SRC_DIR}/PeleLMEB.cpp
       ${SRC_DIR}/PeleLMEos.cpp
       ${SRC_DIR}/PeleLMEvaluate.cpp
       ${SRC_DIR}/PeleLMEvolve.cpp
       ${SRC_DIR}/PeleLMFlowController.cpp
       ${SRC_DIR}/PeleLMForces.cpp
       ${SRC_DIR}/PeleLMInit.cpp
       ${SRC_DIR}/PeleLMPlot.cpp
       ${SRC_DIR}/PeleLMProjection.cpp
       ${SRC_DIR}/PeleLMReactions.cpp
       ${SRC_DIR}/PeleLMRegrid.cpp
       ${SRC_DIR}/PeleLMSetup.cpp
       ${SRC_DIR}/PeleLMTagging.cpp
       ${SRC_DIR}/PeleLMTemporals.cpp
       ${SRC_DIR}/PeleLMTimestep.cpp
       ${SRC_DIR}/PeleLMTransportProp.cpp
       ${SRC_DIR}/PeleLMUMac.cpp
       ${SRC_DIR}/PeleLMUserKeys.H
       ${SRC_DIR}/PeleLMUtils.H
       ${SRC_DIR}/PeleLMUtils.cpp
       ${SRC_DIR}/PeleLM_Index.H
       ${SRC_DIR}/PeleLM_K.H
       ${SRC_DIR}/Utils.cpp
       ${SRC_DIR}/main.cpp
  )

  #if(PELELMEX_ENABLE_ASCENT)
  #  target_sources(${pelelmex_exe_name}
  #    PRIVATE
  #      ${SRC_DIR}/PeleAscent.H
  #      ${SRC_DIR}/PeleAscent.cpp
  #  )
  #endif()

  #if(PELELMEX_ENABLE_MASA)
  #  target_compile_definitions(${pelelmex_exe_name} PRIVATE PELELMEX_USE_MASA)
  #  target_sources(${pelelmex_exe_name} PRIVATE ${SRC_DIR}/MMS.cpp)
  #  target_link_libraries(${pelelmex_exe_name} PRIVATE MASA::MASA)
  #  if(PELELMEX_ENABLE_FPE_TRAP_FOR_TESTS)
  #    set_source_files_properties(${SRC_DIR}/PeleLMeX.cpp PROPERTIES COMPILE_DEFINITIONS PELELMEX_ENABLE_FPE_TRAP)
  #  endif()
  #endif()

  if(NOT "${pelelmex_exe_name}" STREQUAL "PeleLMeX-UnitTests")
    target_sources(${pelelmex_exe_name}
       PRIVATE
         ${CMAKE_SOURCE_DIR}/Source/main.cpp
    )
  endif()

  if(PELELMEX_ENABLE_CUDA)
    set(pctargets "${pelelmex_exe_name};${pelelmex_lib_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(PELELMEX_SOURCES ${tgt} SOURCES)
      list(FILTER PELELMEX_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${PELELMEX_SOURCES} PROPERTIES LANGUAGE CUDA)
    endforeach()
    set_target_properties(${pelelmex_exe_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(${pelelmex_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas --disable-optimizer-constants>)
  endif()
 
  target_link_libraries(${pelelmex_exe_name} PRIVATE ${pelelmex_lib_name} AMReX-Hydro::amrex_hydro_api AMReX::amrex)

  #Define what we want to be installed during a make install 
  install(TARGETS ${pelelmex_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
