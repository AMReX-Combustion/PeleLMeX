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
  const amrex::Real xEBface[AMREX_SPACEDIM],
  amrex::Real s_ext[NVAR],
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/)
{
  if (xEBface[1] > 0.0) {
    s_ext[TEMP] = 800.0;
  } else {
    s_ext[TEMP] = 300.0;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBType(
  const amrex::Real* /*xEBface[AMREX_SPACEDIM]*/,
  amrex::Real& EBflagType,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/)
{
  EBflagType = 1.0;
}
#endif
#endif
