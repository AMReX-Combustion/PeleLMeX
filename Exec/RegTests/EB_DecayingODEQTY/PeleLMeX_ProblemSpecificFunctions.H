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
      a_ode_names[n] = "MY_ODE_" + std::to_string(n);
    }
#endif
}

static void problem_modify_ext_sources(
    Real /*time*/,
    Real /*dt*/,
    int lev,
    MultiArray4<const Real> const& state_old_arr,
    MultiArray4<const Real> const& /*state_new_arr*/,
    Vector<std::unique_ptr<MultiFab>>& a_extSource,
    const GeometryData& /*geomdata*/,
    ProbParm const& prob_parm)
{
  /* 
  Notes: 
    1) a_extSource contains sources from velocity forcing coming in.
       This function should add to rather than overwrite a_extSource.
    2) Requires peleLM.user_defined_ext_sources = true in input file
     
  TODO: 
    1) Add and test capabilities with species, velocity and enthalpy.
  */

  auto ext_source_arr = a_extSource[lev]->arrays();
  amrex::ParallelFor(
    *a_extSource[lev], 
    [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      for (int n = 0; n < NUM_ODE; n++){
        amrex::Real B_n = state_old_arr[box_no](i, j, k, FIRSTODE + n);
        ext_source_arr[box_no](i, j, k, FIRSTODE + n) += prob_parm.ode_srcstrength*B_n;
      }
    });
  amrex::Gpu::streamSynchronize();
}

#endif