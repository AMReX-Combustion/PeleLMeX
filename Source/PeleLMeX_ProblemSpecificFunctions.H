#ifndef PROBLEMSPECIFICFUNCTIONS_H
#define PROBLEMSPECIFICFUNCTIONS_H

#include <PeleLMeX_Index.H>
#include <PeleLMeX.H>

using namespace amrex;

/* 
Problem specific functions:
- This file must be copied locally to the case directory
- Add the following to GNUmakefile: CEXE_headers += ProblemSpecificFunctions.H
- Only copy the functions you want to use
*/

static void set_ode_names(Vector<std::string>& a_ode_names)
{
#if NUM_ODE > 0
    a_ode_names.resize(NUM_ODE);
    for (int n = 0; n < NUM_ODE; n++) {
      a_ode_names[n] = "ODE_" + std::to_string(n);
    }
#endif
}

static void problem_modify_ext_sources(
    Real /*time*/,
    Real /*dt*/,
    int /*lev*/,
    MultiArray4<const Real> const& /*state_old*/,
    MultiArray4<const Real> const& /*state_new*/,
    Vector<std::unique_ptr<MultiFab>>& /*a_extSource*/,
    const GeometryData& /*geomdata*/,
    ProbParm const& /*prob_parm*/)
{
  /*
  Notes: 
    1) a_extSource contains sources from velocity forcing coming in.
       This function should add to rather than overwrite a_extSource.
    2) Requires compilation with USE_MODIFIED_SOURCES = TRUE
     
  TODO: 
    1) Test geomdata, prob_parm
    2) Add and test capabilities with species, velocity and enthalpy.

  // Example: adding a exponential decay source term to an ode_qty 
  auto ext_source_arr = a_extSource[lev]->arrays();
  amrex::ParallelFor(
    *a_extSource[lev], 
    [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      for (int n = 0; n < NUM_ODE; n++){
        amrex::Real B_n = state_old[box_no](i, j, k, FIRSTODE + n);
        ext_source_arr[box_no](i, j, k, FIRSTODE + n) += -1000.*B_n;
      }
    });
  amrex::Gpu::streamSynchronize();
  */
}

#endif