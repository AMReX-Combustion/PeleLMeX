#include <PeleLMeX.H>
#include <PeleLMeX_K.H>
#include <PeleLMeX_PatchFlowVariables.H>

using namespace amrex;

void
patchFlowVariables(
  const amrex::Geometry& /*geom*/,
  ProbParm const& /*prob_parm*/,
  amrex::MultiFab& /*a_mf*/)
{
  Abort("Using patchFlowVariables requires providing a definition in local "
        "PatchFlowVariables.cpp");
}
