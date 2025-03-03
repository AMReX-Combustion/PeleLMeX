#ifndef PELELM_EF_H
#define PELELM_EF_H
// included in the PeleLM.H

// E-field with EB capability is not fully implemented
#ifdef AMREX_USE_EB
static_assert(
  !(PELE_USE_PLASMA),
  "E-field with EB capability is not fully implemented; cannot compile with "
  "both PLASMA and EB");
#endif

// Algorithm int/real data holder
int ef_verbose = 1;
int ef_substep = 1;       // number of sub time steps
amrex::Real dtsub;        // sub time step size
amrex::Real curtime;      // time of current substep
amrex::Real nl_stateNorm; // norm of the non-linear state
amrex::Real nl_residNorm; // norm of the non-linear residual

// Restart options
int m_restart_nonEF = 0;
int m_restart_electroneutral = 1;
int m_restart_resetTime = 1;

// Scheme
int m_nEAdvOrder = 2;

// Electron transport
int m_electronKappaTab = 0;
amrex::Real m_fixedKappaE = 0.4;

// Scaling of non-linear system
amrex::Real nE_scale;          // nE state scale
amrex::Real phiV_scale;        // phiV state scale
amrex::Real FnE_scale = 1.0;   // nE residual scale
amrex::Real FphiV_scale = 1.0; // phiV residual scale

// non-linear / linear solve
int m_ef_maxNewtonIter = 5;
int m_ef_use_PETSC_direct = 0;
int m_ef_PC_approx = 1;
int m_ef_diffT_jfnk = 1;
amrex::Real m_ABecCecOmega = 0.9;
amrex::Real m_ef_lambda_jfnk = 1.0e-7;
amrex::Real m_ef_newtonTol = std::pow(1.0e-13, 2.0 / 3.0);
amrex::Real m_ef_GMRES_reltol = 1.0e-6;
amrex::Real m_ef_GMRES_abstol = 1.0e-12;
amrex::Real m_ef_PC_MG_Tol = 1.0e-6;

// Species charge per unit mass
amrex::GpuArray<amrex::Real, NUM_SPECIES> zk;

// NLSolve level data container
struct LevelDataNLSolve
{
  LevelDataNLSolve() = default;
  LevelDataNLSolve(
    amrex::BoxArray const& ba,
    amrex::DistributionMapping const& dm,
    amrex::FabFactory<amrex::FArrayBox> const& factory,
    int a_nGrow);
  amrex::MultiFab nlState;          // Non-linear state
  amrex::MultiFab nlResid;          // Non-linear residual
  amrex::MultiFab backgroundCharge; // Background charge from the ions
  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>
    gPhiVOld; // Gradient of old electro-static pot.
  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>
    uEffnE; // Effective velocity of the electrons
  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>
    umac; // Need a duplicate of umac ... TODO find better way
};

// Preconditioner
std::unique_ptr<PrecondOp> m_precond_op;
PrecondOp* getPrecondOp();

amrex::Vector<std::unique_ptr<LevelDataNLSolve>> m_leveldatanlsolve;
LevelDataNLSolve*
getLevelDataNLSolvePtr(int lev)
{
  return m_leveldatanlsolve[lev].get();
};

amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>
getNLgradPhiVVect();
amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>> getUeffVect();
amrex::Vector<amrex::MultiFab*> getNLresidVect();
amrex::Vector<amrex::MultiFab*> getNLstateVect();
amrex::Vector<amrex::MultiFab*> getNLBGChargeVect();

// Temporals
int m_do_ionsBalance = 0;
std::ofstream tmpIonsFile;

// Extra diagnostics
int m_do_extraEFdiags = 0;
amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_ionsFluxes;
#endif
