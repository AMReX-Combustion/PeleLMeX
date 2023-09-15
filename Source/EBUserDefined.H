#ifndef EBUSERDEFINED_H
#define EBUSERDEFINED_H

using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
void
EBUserDefined(
  const Geometry& /*geom*/,
  const int /*required_coarsening_level*/,
  const int /*max_coarsening_level*/)
{
  // ParmParse your geometry parameters

  // Build geometry pieces using EB2::* methods

  // Build your geometry shop using EB2::makeShop

  // Build geom using EB2::Build

  // We shoulnd't be here, copy this file in your run folder
  // and implement your geometry
  Abort("In default EBUserDefined function! Shouldn't be here. Copy and edit "
        "this file for your needs");
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBState(
  const amrex::Real* /*xEBface[AMREX_SPACEDIM]*/,
  amrex::Real* /*s_ext[NVAR]*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/)
{
  // Provide the state on the EB centroid provided
  // by xEBface
  // Not the entire has to be specified. For instance, if using isothermal
  // boundary, only s_ext[TEMP] will be used.
  // Aborting here (message will not show on GPUs)
  Abort("Default setEBState called. Copy and edit this file for your needs");
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBType(
  const amrex::Real* /*xEBface[AMREX_SPACEDIM]*/,
  amrex::Real& /*EBflagType*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/)
{
  // Provide a factor applied on the EB face diffusion
  // coefficient: 0.0 will turn off the diffusion (homogeneous Neumann BC)
  // while 1.0 will let it unchanged. Intermediate values can be specified
  // for partially Dirichlet BCs
  Abort("Default setEBType called. Copy and edit this file for your needs");
}
#endif
#endif
