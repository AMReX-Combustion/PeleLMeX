target_include_directories(pelelmex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(pelelmex
    PRIVATE
        AMReX_MLABecCecLap_${AMReX_SPACEDIM}D_K.H
        AMReX_MLABecCecLap_K.H
        AMReX_MLABecCecLaplacian.cpp
        AMReX_MLABecCecLaplacian.H
        AMReX_MLCellABecCecLap.cpp
        AMReX_MLCellABecCecLap.H
        AMReX_MLCellABecCecLap_${AMReX_SPACEDIM}D_K.H
        AMReX_MLCellABecCecLap_K.H
        PrecondOp.H
        PrecondOp.cpp
    )
