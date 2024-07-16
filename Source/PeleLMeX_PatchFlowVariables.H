#ifndef PELELM_PATCHFLOWVARIABLES_H
#define PELELM_PATCHFLOWVARIABLES_H

class PeleLM;

void patchFlowVariables(
  const amrex::Geometry& geom,
  ProbParm const& prob_parm,
  amrex::MultiFab& a_mf);
#endif
