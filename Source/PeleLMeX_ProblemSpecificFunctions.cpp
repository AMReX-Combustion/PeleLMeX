#include <PeleLMeX.H>
#include <PeleLMeX_K.H>
#include <PeleLMeX_ProblemSpecificFunctions.H>

using namespace amrex;

/*
Problem specific functions:
- This file must be copied locally to the case directory
- Add the following to GNUmakefile:
          CEXE_sources += PeleLMeX_ProblemSpecificFunctions.cpp
- Modify as needed
*/

#if NUM_ODE > 0
void
set_ode_names(Vector<std::string>& a_ode_names)
{
  a_ode_names.resize(NUM_ODE);
  for (int n = 0; n < NUM_ODE; n++) {
    a_ode_names[n] = "ODE_" + std::to_string(n);
  }
}
#endif

void
problem_modify_ext_sources(
  Real /*time*/,
  Real /*dt*/,
  int /*lev*/,
  MultiArray4<const Real> const& /*state_old_arr*/,
  MultiArray4<const Real> const& /*state_new_arr*/,
  Vector<std::unique_ptr<MultiFab>>& /*a_extSource*/,
  const GeometryData& /*geomdata*/,
  ProbParm const& /*prob_parm*/)
{
  /*
  Notes:
    1) a_extSource contains sources from velocity forcing coming in.
       This function should add to rather than overwrite a_extSource.
    2) Requires peleLM.user_defined_ext_sources = true in input file

  // Example: Exponential decay ode quantity
  auto ext_source_arr = a_extSource[lev]->arrays();
  ParallelFor(
    *a_extSource[lev],
    [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      for (int n = 0; n < NUM_ODE; n++){
        Real B_n = state_old_arr[box_no](i, j, k, FIRSTODE + n);
        Real src_strength = -1.0 * pow(10.0,n+1);
        ext_source_arr[box_no](i, j, k, FIRSTODE + n) += src_strength * B_n;
      }
    });
  Gpu::streamSynchronize();
  */
}