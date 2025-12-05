#ifndef PELELM_UTILS_H
#define PELELM_UTILS_H

#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

template <class T>
amrex::Gpu::DeviceVector<T>
convertToDeviceVector(const amrex::Vector<T>& v)
{
  static_assert(
    std::is_trivially_copyable_v<T>, "T must be trivially copyable");
  int ncomp = v.size();
  amrex::Gpu::DeviceVector<T> v_d(ncomp);
#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(v_d.data(), v.data(), sizeof(T) * ncomp);
#else
  std::memcpy(v_d.data(), v.data(), sizeof(T) * ncomp);
#endif
  return v_d;
}

void
pltMLMGResidual(
  amrex::MLMG& a_mlmg,
  const amrex::Vector<amrex::MultiFab*>& a_sol,
  const amrex::Vector<const amrex::MultiFab*>& a_rhs,
  const std::string& a_solver_name,
  const std::string& a_plot_file_root,
  int a_step,
  const amrex::Vector<amrex::Geometry>& a_geom,
  class PeleLM* a_pelelm);

#endif
