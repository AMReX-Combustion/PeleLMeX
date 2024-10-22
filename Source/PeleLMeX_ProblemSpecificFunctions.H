#ifndef PROBLEMSPECIFICFUNCTIONS_H
#define PROBLEMSPECIFICFUNCTIONS_H

using namespace amrex;

class PeleLM;

void set_ode_names(Vector<std::string>& a_ode_names);

void problem_modify_ext_sources(
  Real time,
  Real dt,
  const MultiFab& state_old,
  const MultiFab& state_new,
  std::unique_ptr<MultiFab>& ext_src,
  const GeometryData& geomdata,
  const ProbParm& prob_parm);
#endif