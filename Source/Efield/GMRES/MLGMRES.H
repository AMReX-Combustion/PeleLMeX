#ifndef MLGMRES_H
#define MLGMRES_H
#include <AMReX_REAL.H>
#include <PeleLMeX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>

#define MEMBER_FUNC_PTR(object, ptrToMember) ((object).*(ptrToMember))

typedef void (PeleLM::*MLJtimesVFunc)(
  const amrex::Vector<amrex::MultiFab*>& a_x,
  const amrex::Vector<amrex::MultiFab*>& a_Ax);
typedef void (PeleLM::*MLPrecondFunc)(
  const amrex::Vector<amrex::MultiFab*>& v,
  const amrex::Vector<amrex::MultiFab*>& Pv);
typedef void (PeleLM::*MLNormFunc)(
  const amrex::Vector<amrex::MultiFab*>& MFin, amrex::Real& norm);

class MLGMRESSolver
{
public:
  MLGMRESSolver();
  virtual ~MLGMRESSolver();

  MLGMRESSolver(const MLGMRESSolver&) = delete;
  MLGMRESSolver& operator=(const MLGMRESSolver&) = delete;
  MLGMRESSolver& operator=(MLGMRESSolver&&) = delete;

  void readParameters();
  void define(PeleLM* a_level, const int a_nComp, const int a_nGrow);

  void setJtimesV(MLJtimesVFunc a_jtv);
  void setPrecond(MLPrecondFunc a_prec);
  void setNorm(MLNormFunc a_norm);

  void setVerbose(int a_v) { m_verbose = a_v; };
  void setMaxRestart(int a_maxRestart) { m_restart = a_maxRestart; };

  MLJtimesVFunc jtimesv() const noexcept;
  MLPrecondFunc precond() const noexcept;
  MLNormFunc norm() const noexcept;

  int solve(
    const amrex::Vector<amrex::MultiFab*>& sol,
    const amrex::Vector<amrex::MultiFab*>& rhs,
    amrex::Real a_abs_tol,
    amrex::Real a_rel_tol);

private:
  void one_restart(
    const amrex::Vector<amrex::MultiFab*>& a_x,
    const amrex::Vector<amrex::MultiFab*>& a_rhs);

  void one_iter(const int iter, amrex::Real& resNorm);

  void appendBasisVector(
    const int iter, amrex::Vector<amrex::Vector<amrex::MultiFab>>& Base);

  void gramSchmidtOrtho(
    const int iter, amrex::Vector<amrex::Vector<amrex::MultiFab>>& Base);

  amrex::Real givensRotation(const int iter);

  void prepareForSolve();

  void computeMLResidual(
    const amrex::Vector<amrex::MultiFab*>& a_x,
    const amrex::Vector<amrex::MultiFab*>& a_rhs,
    const amrex::Vector<amrex::MultiFab*>& a_res);

  amrex::Real computeMLResidualNorm(
    const amrex::Vector<amrex::MultiFab*>& a_x,
    const amrex::Vector<amrex::MultiFab*>& a_rhs);

  amrex::Real computeMLNorm(const amrex::Vector<amrex::MultiFab*>& a_vec);

  void
  updateSolution(int final_iter, const amrex::Vector<amrex::MultiFab*>& a_x);

  amrex::Real MFVecDot(
    const amrex::Vector<const amrex::MultiFab*>& a_mf1,
    int mf1comp,
    const amrex::Vector<const amrex::MultiFab*>& a_mf2,
    int mf2comp,
    int nComp,
    int nGrow);

  void MFVecSaxpy(
    const amrex::Vector<amrex::MultiFab*>& a_mfdest,
    amrex::Real a_a,
    const amrex::Vector<const amrex::MultiFab*>& a_mfsrc,
    int destComp,
    int srcComp,
    int nComp,
    int nGrow);

  int m_verbose;

  //   GMRES attributes
  int check_GramSchmidtOrtho = 1;
  bool m_converged;
  int m_krylovSize = 10;
  int m_restart = 2;
  int m_nComp = -1;
  amrex::Real initResNorm;
  amrex::Real target_relResNorm;
  amrex::Real target_absResNorm;
  int iter_count;
  int restart_count;

  //   GMRES data
  //   MultiFab data
  amrex::Vector<amrex::Vector<amrex::MultiFab>> KspBase; // Krylov basis
  amrex::Vector<amrex::MultiFab> Ax;                     // A container for A*x
  amrex::Vector<amrex::MultiFab> res; // A container for residual

  //   Real data
  amrex::Vector<amrex::Vector<amrex::Real>> H;      // Hessenberg matrix
  amrex::Vector<amrex::Real> y;                     // Solution vector
  amrex::Vector<amrex::Real> g;                     // Residual
  amrex::Vector<amrex::Vector<amrex::Real>> givens; // Givens rotation

  //   AMReX attributes
  PeleLM* m_pelelm;
  int m_nGrow;
  amrex::Vector<amrex::Geometry> m_geom;
  amrex::Vector<amrex::BoxArray> m_grids;
  amrex::Vector<amrex::DistributionMapping> m_dmap;

  //   GMRES user-provided functions
  MLJtimesVFunc m_jtv = nullptr;  // JtimesV function. Always needed
  MLPrecondFunc m_prec = nullptr; // Preconditioner
  MLNormFunc m_norm = nullptr; // User-prefered MF norm. L2 used if not provided
};

#endif
