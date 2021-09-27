#include <PeleLM.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

using namespace amrex;

void PeleLM::makeEBGeometry()
{
    int max_coarsening_level = 100;
    int req_coarsening_level = 2;
    EB2::Build(geom.back(),req_coarsening_level,max_coarsening_level);
}

#endif
