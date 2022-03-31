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
               const amrex::Vector<amrex::BoxArray> &a_grids,
               const amrex::Vector<amrex::DistributionMapping> &a_dmap) override;

  void processDiag(const amrex::Real &a_time,
                   const amrex::Vector<const amrex::MultiFab*> &a_state) override;

  void close() override {}

private:
  int m_normal;
  amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> m_center;

  // Interpolation data
  amrex::Vector<amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> > m_intwgt;
  amrex::Vector<int> m_k0;

  // 2D-plane boxArray vector
  amrex::Geometry m_geomLev0;
  amrex::Vector<amrex::BoxArray> m_sliceBA;
  amrex::Vector<amrex::Vector<int> > m_dmConvert;
  amrex::Vector<amrex::DistributionMapping> m_sliceDM;
};

#endif
