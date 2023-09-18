#ifndef PRECOND_OP_H
#define PRECOND_OP_H

#include <AMReX_MLMG.H>
#include <AMReX_BC_TYPES.H>

#ifdef AMREX_USE_EB
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLEBTensorOp.H>
#endif
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLABecCecLaplacian.H>
#include <AMReX_MLTensorOp.H>

class PeleLM;

class PrecondOp
{
public:
  PrecondOp(PeleLM* a_pelelm);

  // Data and parameters
  PeleLM* m_pelelm;

  void setDiffOpScalars(
    const amrex::Real& aScal,
    const amrex::Real& bScal,
    const amrex::Real& cScal);
  void setDriftOpScalars(const amrex::Real& aScal, const amrex::Real& bScal);
  void setStildaOpScalars(const amrex::Real& aScal, const amrex::Real& bScal);

  void setDiffOpBCs(const amrex::BCRec& a_bcrec);
  void setDriftOpBCs(const amrex::BCRec& a_bcrec);
  void setStildaOpBCs(const amrex::BCRec& a_bcrec);

  void setDiffOpLevelBC(int lev, amrex::MultiFab* a_MF);
  void setDriftOpLevelBC(int lev, amrex::MultiFab* a_MF);
  void setStildaOpLevelBC(int lev, amrex::MultiFab* a_MF);

  void setDiffOpRelaxation(const amrex::Real& a_omega);
  void setDiffOpACoeff(int lev, const amrex::Real& acoeff);
  void setDiffOpBCoeff(
    int lev,
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_bCoeff);
  void setDiffOpCCoeff(
    int lev,
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_cCoeff);

  void getDiffOpDiagonal(int lev, amrex::MultiFab& a_diffOpDiag);

  void setDriftOpBCoeff(
    int lev,
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_bCoeff);
  void setStildaOpBCoeff(
    int lev,
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_bCoeff);

  void diffOpSolve(
    const amrex::Vector<amrex::MultiFab*>& a_sol,
    const amrex::Vector<const amrex::MultiFab*>& a_rhs,
    const amrex::Real& rtol,
    const amrex::Real& atol);

  void StildaOpSolve(
    const amrex::Vector<amrex::MultiFab*>& a_sol,
    const amrex::Vector<const amrex::MultiFab*>& a_rhs,
    const amrex::Real& rtol,
    const amrex::Real& atol);

  void driftOpApply(
    const amrex::Vector<amrex::MultiFab*>& a_Ax,
    const amrex::Vector<amrex::MultiFab*>& a_x);

  void readParameters();

  amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>
  getOpBC(amrex::Orientation::Side a_side, const amrex::BCRec& a_bc);

  // PC linear ops
  std::unique_ptr<amrex::MLABecCecLaplacian> m_diff;
  std::unique_ptr<amrex::MLABecLaplacian> m_drift;
  std::unique_ptr<amrex::MLABecLaplacian> m_Stilda;

  // Associated MLMGs
  std::unique_ptr<amrex::MLMG> m_mlmg_diff;
  std::unique_ptr<amrex::MLMG> m_mlmg_drift;
  std::unique_ptr<amrex::MLMG> m_mlmg_Stilda;

  // Verbose for each Ops
  int m_diff_verbose = 0;
  int m_drift_verbose = 0;
  int m_Stilda_verbose = 0;

  // Options to control MLMG behavior
  int m_fixed_mg_it = -1;
  int m_mg_verbose = 0;
  int m_mg_bottom_verbose = 0;
  int m_mg_max_iter = 100;
  int m_mg_bottom_maxiter = 100;
  int m_mg_max_fmg_iter = 0;
  int m_mg_max_coarsening_level_diff = 100;
  int m_mg_max_coarsening_level_Stilda = 100;
  int m_mg_maxorder = 2;
  int m_num_pre_smooth = 2;
  int m_num_post_smooth = 2;

  amrex::Real m_mg_rtol = 1.0e-11;
  amrex::Real m_mg_atol = 1.0e-14;

  std::string m_mg_bottom_solver = "bicgstab";
};
#endif
