target_include_directories(pelelmex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(pelelmex
    PRIVATE
        main.cpp
        DeriveUserDefined.cpp
        DiffusionOp.H
        DiffusionOp.cpp
        PeleLM.cpp
        PeleLM.H
        PeleLMeX_Index.H
        PeleLMeX_K.H
        PeleLMAdvance.cpp
        PeleLMAdvection.cpp
        PeleLMBC.cpp
        PeleLMBCfill.H
        PeleLMData.cpp
        PeleLMDerive.cpp
        PeleLMDerive.H
        PeleLMDeriveFunc.cpp
        PeleLMDeriveFunc.H
        PeleLMDiagnostics.cpp
        PeleLMDiffusion.cpp
        PeleLMEB.cpp
        PeleLMEos.cpp
        PeleLMEvaluate.cpp
        PeleLMEvolve.cpp
        PeleLMForces.cpp
        PeleLMInit.cpp
        PeleLMPlot.cpp
        PeleLMProjection.cpp
        PeleLMReactions.cpp
        PeleLMRegrid.cpp
        PeleLMSetup.cpp
        PeleLMTagging.cpp
        PeleLMTemporals.cpp
        PeleLMTimestep.cpp
        PeleLMTransportProp.cpp
        PeleLMUMac.cpp
        PeleLMUtils.H
        PeleLMUtils.cpp
        PeleFlowControllerData.H
        PeleLMFlowController.cpp
        Utils.cpp
    )

if (PELELMEX_USE_EFIELD)
    add_subdirectory(Efield)
endif ()

if (PELELMEX_USE_PARTICLES)
    add_subdirectory(Spray)
endif ()

if (PELELMEX_USE_SOOT)
    add_subdirectory(Soot)
endif ()
