target_include_directories(pelelmex PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(pelelmex
    PRIVATE
        DiagBase.cpp
        DiagBase.H
        DiagConditional.cpp
        DiagConditional.H
        DiagFilter.cpp
        DiagFilter.H
        DiagFramePlane.cpp
        DiagFramePlane.H
        DiagPDF.cpp
        DiagPDF.H
    )
