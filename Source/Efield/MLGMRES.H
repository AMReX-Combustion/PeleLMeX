#ifndef MLGMRES_H_
#define MLGMRES_H_
#include <AMReX_REAL.H>
#include <PeleLM.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>

#define MEMBER_FUNC_PTR(object,ptrToMember)  ((object).*(ptrToMember))

//typedef void (PeleLM::*JtimesVFunc) (const amrex::MultiFab& a_x, amrex::MultiFab& a_Ax);
//typedef void (PeleLM::*PrecondFunc) (const amrex::MultiFab& v, amrex::MultiFab& Pv);
//typedef void (PeleLM::*NormFunc)    (const amrex::MultiFab& MFin, amrex::Real& norm);
typedef void (PeleLM::*MLJtimesVFunc) (const amrex::Vector<amrex::MultiFab*> &a_x, const amrex::Vector<amrex::MultiFab*> &a_Ax);
typedef void (PeleLM::*MLPrecondFunc) (const amrex::Vector<amrex::MultiFab*> &v, const amrex::Vector<amrex::MultiFab*> &Pv);
typedef void (PeleLM::*MLNormFunc)    (const amrex::Vector<amrex::MultiFab*> &MFin, amrex::Real& norm);

class MLGMRESSolver
{
   public:
     MLGMRESSolver();
     virtual ~MLGMRESSolver();

     MLGMRESSolver (const MLGMRESSolver&) = delete;
     MLGMRESSolver& operator= (const MLGMRESSolver&) = delete;
     MLGMRESSolver& operator= (MLGMRESSolver&&) = delete;

     void define(PeleLM* a_level,
                 const int a_KrylovSize,
                 const int a_nComp,
                 const int a_nGrow);

     void setJtimesV(MLJtimesVFunc a_jtv);
     void setPrecond(MLPrecondFunc a_prec);
     void setNorm(MLNormFunc a_norm);

     void setVerbose(int a_v){ m_verbose = a_v;};
     void setMaxRestart(int a_maxRestart){ m_restart = a_maxRestart;};

     MLJtimesVFunc jtimesv() const noexcept;
     MLPrecondFunc precond() const noexcept;
     MLNormFunc    norm() const noexcept;

     int solve(const amrex::Vector<amrex::MultiFab*> &sol,
               const amrex::Vector<amrex::MultiFab*> &rhs,
               amrex::Real a_abs_tol,
               amrex::Real a_rel_tol);

   private:

     void one_restart(const amrex::Vector<amrex::MultiFab*> &a_x,
                      const amrex::Vector<amrex::MultiFab*> &a_rhs);

     void one_iter(const int iter, amrex::Real &resNorm);

     void appendBasisVector(const int iter, amrex::Vector<amrex::MultiFab>& Base);

     void gramSchmidtOrtho(const int iter, amrex::Vector<amrex::MultiFab>& Base);

     amrex::Real givensRotation(const int iter);

     void prepareForSolve();

     void computeResidual(const amrex::MultiFab& a_x,
                          const amrex::MultiFab& a_rhs,
                          amrex::MultiFab& a_res);

     void computeMLResidual(const amrex::Vector<amrex::MultiFab*> &a_x,
                            const amrex::Vector<amrex::MultiFab*> &a_rhs,
                            const amrex::Vector<amrex::MultiFab*> &a_res);

     amrex::Real computeMLResidualNorm(const amrex::Vector<amrex::MultiFab*> &a_x,
                                       const amrex::Vector<amrex::MultiFab*> &a_rhs);

     amrex::Real computeResidualNorm(const amrex::MultiFab& a_x,
                                     const amrex::MultiFab& a_rhs);

     amrex::Real computeMLNorm(const amrex::Vector<amrex::MultiFab*> &a_vec);

     amrex::Real computeNorm(const amrex::MultiFab& a_vec);

     int m_verbose;

//   GMRES attributes
     int check_GramSchmidtOrtho = 1;
     bool m_converged;
     int m_krylovSize;
     int m_restart;
     int m_nComp;
     amrex::Real initResNorm;
     amrex::Real targetResNorm;
     int iter_count;
     int restart_count;

//   GMRES data
//   MultiFab data
     amrex::Vector<amrex::Vector<amrex::MultiFab>> KspBase;  // Krylov basis
     amrex::Vector<amrex::MultiFab> Ax;                      // A container for A*x
     amrex::Vector<amrex::MultiFab> res;                     // A container for residual

//   Real data
     amrex::Vector<amrex::Vector<amrex::Real>> H;        // Hessenberg matrix
     amrex::Vector<amrex::Real>                y;        // Solution vector
     amrex::Vector<amrex::Real>                g;        // Residual
     amrex::Vector<amrex::Vector<amrex::Real>> givens;   // Givens rotation

//   AMReX attributes
     PeleLM* m_pelelm;
     int m_nGrow;
     amrex::Vector<amrex::Geometry> m_geom;
     amrex::Vector<amrex::BoxArray> m_grids;
     amrex::Vector<amrex::DistributionMapping> m_dmap;

//   GMRES user-provided functions
     MLJtimesVFunc m_jtv = nullptr;     // JtimesV function. Always needed
     MLPrecondFunc m_prec = nullptr;    // Preconditioner
     MLNormFunc m_norm = nullptr;       // User-prefered MF norm. L2 used if not provided
};

#endif
