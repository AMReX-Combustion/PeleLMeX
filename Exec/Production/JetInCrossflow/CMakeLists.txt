target_include_directories(pelelmex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

set(PELELMEX_CHEMISTRY_MODEL LiDryer PARENT_SCOPE)
set(PELELMEX_EOS_MODEL       Fuego PARENT_SCOPE)
set(PELELMEX_TRANSPORT_MODEL Simple PARENT_SCOPE)
set(PELELMEX_DIM "3" PARENT_SCOPE)
set(PELELMEX_MPI ON PARENT_SCOPE)

target_sources(pelelmex
    PRIVATE
        pelelmex_prob.cpp
        pelelmex_prob.H
        pelelmex_prob_parm.H
    )
