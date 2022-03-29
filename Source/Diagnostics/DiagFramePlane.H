#ifndef _DIAGFRAMEPLANE_H_
#define _DIAGFRAMEPLANE_H_

#include "DiagBase.H"

class DiagFramePlane : public DiagBase::Register<DiagFramePlane>
{
public:
  static std::string identifier() { return "DiagFramePlane"; }

  void init(const std::string &a_prefix) override;

  void prepare(int a_nlevels,
               const amrex::Vector<amrex::Geometry> &a_geoms,
               const amrex::Vector<amrex::BoxArray> &a_grids) override;

  void processDiag(const amrex::Real &a_time,
                   const amrex::Vector<const amrex::MultiFab*> &a_state) override;

  void close() override {}

private:
  int m_normal;
  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> m_center;
};

#endif
