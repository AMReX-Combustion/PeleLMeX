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
    a_ode_names[n] = "MY_ODE_" + std::to_string(n);
  }
}
#endif

void
problem_modify_ext_sources(
  Real /*time*/,
  Real /*dt*/,
  const MultiFab& state_old,
  const MultiFab& /*state_new*/,
  std::unique_ptr<MultiFab>& ext_src,
  const GeometryData& /*geomdata*/,
  const ProbParm& prob_parm)
{
  /*
  Notes:
    1) ext_src contains sources from velocity forcing coming in.
       This function should add to rather than overwrite ext_src.
    2) Requires "peleLM.user_defined_ext_sources = true" in input file
  */

  auto ext_src_arr = ext_src->arrays();
  auto const& state_old_arr = state_old.const_arrays();

  ParallelFor(
    *ext_src, [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      for (int n = 0; n < NUM_ODE; n++) {
        Real B_n = state_old_arr[box_no](i, j, k, FIRSTODE + n);
        Real src = prob_parm.ode_srcstrength * pow(10.0, n + 1) * B_n;
        ext_src_arr[box_no](i, j, k, FIRSTODE + n) += src;
      }
    });
  Gpu::streamSynchronize();
}