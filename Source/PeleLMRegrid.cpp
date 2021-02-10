#include <PeleLM.H>

using namespace amrex;

void PeleLM::MakeNewLevelFromCoarse( int lev,
                                           amrex::Real time,
                                     const amrex::BoxArray& ba,
                                     const amrex::DistributionMapping& dm) {
   BL_PROFILE_VAR("PeleLM::MakeNewLevelFromCoarse()", MakeNewLevelFromCoarse);
}

void PeleLM::RemakeLevel( int lev,
                                amrex::Real time,
                          const amrex::BoxArray& ba,
                          const amrex::DistributionMapping& dm) {
   BL_PROFILE_VAR("PeleLM::RemakeLevel()", RemakeLevel);
}

void PeleLM::ClearLevel(int lev) {
   BL_PROFILE_VAR("PeleLM::ClearLevel()", ClearLevel);
}
