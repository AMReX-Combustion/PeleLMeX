#ifndef PROBLEMSPECIFICFUNCTIONS_H
#define PROBLEMSPECIFICFUNCTIONS_H

using namespace amrex;

#include <PeleLMeX_Index.H>
#include <PeleLMeX.H>

void set_ode_names(amrex::Vector<std::string>& a_ode_names)
{
    a_ode_names.resize(NUM_ODE);
#if NUM_ODE > 0
    for (int n = 0; n < NUM_ODE; n++) {
      a_ode_names[n] = "ODE_" + std::to_string(n);
    }
#endif
}

void problem_modify_ext_sources(
    amrex::Real /*time*/,
    amrex::Real /*dt*/,
    const amrex::MultiFab& /*state_old*/,
    const amrex::MultiFab& /*state_new*/,
    amrex::MultiFab& /*ext_src*/,
    int /*ng*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParm const& /*prob_parm*/)
{
  /* Notes: ext_src contains sources from velocity forcing coming in
      This function should add to rather than overwrite ext_src.
  */
}

#endif