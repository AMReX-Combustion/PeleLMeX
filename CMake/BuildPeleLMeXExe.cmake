function(build_pelelmex_exe pelelmex_exe_name pelelmex_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)

  add_executable(${pelelmex_exe_name} "")

  if(CLANG_TIDY_EXE)
    set_target_properties(${pelelmex_exe_name} PROPERTIES CXX_CLANG_TIDY 
                          "${CLANG_TIDY_EXE};--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy")
  endif()

  target_sources(${pelelmex_exe_name}
     PRIVATE
       pelelmex_prob_parm.H
       pelelmex_prob.H
       pelelmex_prob.cpp
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
                   ${SRC_DIR}/PeleLMeX_SprayParticles.cpp
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
                   ${SRC_DIR}/PeleLMeX_Soot.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel_react.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel_derive.cpp
                   ${PELEMP_SRC_DIR}/Soot_Models/Constants_Soot.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootData.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootReactions.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel.H
                   ${PELEMP_SRC_DIR}/Soot_Models/SootModel_derive.H
                   ${SRC_DIR}/PeleLMeXSoot.cpp)
    target_include_directories(${pelelmex_exe_name} PRIVATE ${PELEMP_SRC_DIR}/Soot_Models)
  endif()

  target_sources(${pelelmex_exe_name}
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
       ${SRC_DIR}/PeleLMeX_Init.cpp
       ${SRC_DIR}/PeleLMeX_Plot.cpp
       ${SRC_DIR}/PeleLMeX_Projection.cpp
       ${SRC_DIR}/PeleLMeX_Reactions.cpp
       ${SRC_DIR}/PeleLMeX_Regrid.cpp
       ${SRC_DIR}/PeleLMeX_Setup.cpp
       ${SRC_DIR}/PeleLMeX_Tagging.cpp
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
