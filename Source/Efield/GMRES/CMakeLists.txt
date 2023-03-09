target_include_directories(pelelmex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(pelelmex
    PRIVATE
        MLGMRES.H
        MLGMRES.cpp
    )
