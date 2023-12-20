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

  # Spray
  set(PELE_PHYSICS_SPRAY_DIR ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Source/Spray)

  if(PELELMEX_ENABLE_AMREX_PARTICLES AND PELE_SPRAY_FUEL_NUM GREATER 0)
     target_compile_definitions(${pelelmex_exe_name} PRIVATE PELELMEX_USE_SPRAY)
     target_compile_definitions(${pelelmex_exe_name} PRIVATE SPRAY_FUEL_NUM=${PELE_SPRAY_FUEL_NUM})
     target_sources(${pelelmex_exe_name} PRIVATE
                    SprayParticlesInitInsert.cpp
                    ${SRC_DIR}/Particle.cpp
                    ${PELE_PHYSICS_SPRAY_DIR}/Drag.H
                    ${PELE_PHYSICS_SPRAY_DIR}/SprayDerive.cpp
                    ${PELE_PHYSICS_SPRAY_DIR}/SprayFuelData.H
                    ${PELE_PHYSICS_SPRAY_DIR}/SprayIO.cpp
                    ${PELE_PHYSICS_SPRAY_DIR}/SprayInjection.H
                    ${PELE_PHYSICS_SPRAY_DIR}/SprayInterpolation.H
                    ${PELE_PHYSICS_SPRAY_DIR}/SprayJet.H
                    ${PELE_PHYSICS_SPRAY_DIR}/SprayJet.cpp
                    ${PELE_PHYSICS_SPRAY_DIR}/SprayParticles.H
                    ${PELE_PHYSICS_SPRAY_DIR}/SprayParticles.cpp
                    ${PELE_PHYSICS_SPRAY_DIR}/SpraySB.cpp
                    ${PELE_PHYSICS_SPRAY_DIR}/SpraySetup.cpp
                    ${PELE_PHYSICS_SPRAY_DIR}/WallFunctions.H
                    ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/AhamedSplash.H
                    ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/ReitzKHRT.H
                    ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/SBData.H
                    ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/TABBreakup.H
                    ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/WallFilm.H
                    ${PELE_PHYSICS_SPRAY_DIR}/Distribution/DistBase.H
                    ${PELE_PHYSICS_SPRAY_DIR}/Distribution/Distributions.H
                    ${PELE_PHYSICS_SPRAY_DIR}/Distribution/Distributions.cpp)
     target_include_directories(${pelelmex_exe_name} PUBLIC ${PELE_PHYSICS_SPRAY_DIR})
     target_include_directories(${pelelmex_exe_name} PUBLIC ${PELE_PHYSICS_SPRAY_DIR}/Distribution)
     target_include_directories(${pelelmex_exe_name} PUBLIC ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash)
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

  # Soot
  if(PELELMEX_ENABLE_SOOT)
    target_sources(${pelelmex_exe_name}
      PRIVATE
        ${SRC_DIR}/PeleLMeX_Soot.cpp)
    if(PELELMEX_ENABLE_RADIATION)
      target_sources(${pelelmex_exe_name}
        PRIVATE
          ${SRC_DIR}/PeleLMeX_Radiation.cpp)
    endif()
  endif()

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
