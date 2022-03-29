#ifndef _DIAGBASE_H_
#define _DIAGBASE_H_

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include "Factory.H"

class DiagBase : public pele::physics::Factory<DiagBase>
{
public:
  static std::string base_identifier() { return "DiagBase"; }

  virtual void init(const std::string &a_prefix) = 0;

  virtual void close() = 0; 

  virtual bool doDiag(const amrex::Real &a_time, int a_nstep);

  virtual void prepare(int a_nlevels,
                       const amrex::Vector<amrex::Geometry> &a_geoms,
                       const amrex::Vector<amrex::BoxArray> &a_grids) = 0;

  virtual void processDiag(const amrex::Real &a_time,
                           const amrex::Vector<const amrex::MultiFab*> &a_state) = 0;

protected:
  int m_verbose = 0;
  amrex::Real m_per = -1.0;
  int m_freq = -1;
};

#endif
