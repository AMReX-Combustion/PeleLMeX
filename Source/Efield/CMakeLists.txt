target_include_directories(pelelmex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(pelelmex
    PRIVATE
        EOS_Extension.H
        PeleLMEF.H
        PeleLMEF_Constants.H
        PeleLMEF_K.H
        PeleLMEFDeriveFunc.cpp
        PeleLMEFDeriveFunc.H
        PeleLMEFIonDrift.cpp
        PeleLMEFNLSolve.cpp
        PeleLMEFPoisson.cpp
        PeleLMEFReactions.cpp
        PeleLMEFTimeStep.cpp
        PeleLMEFTransport.cpp
        PeleLMEFUtils.cpp
    )

add_subdirectory(GMRES)
add_subdirectory(LinOps)
