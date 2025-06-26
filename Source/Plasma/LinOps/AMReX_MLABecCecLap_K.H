#ifndef AMREX_MLABECCECLAP_K_H
#define AMREX_MLABECCECLAP_K_H

#include <AMReX_FArrayBox.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLABecCecLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLABecCecLap_2D_K.H>
#else
#include <AMReX_MLABecCecLap_3D_K.H>
#endif

#endif
