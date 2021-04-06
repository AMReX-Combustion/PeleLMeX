#include <MLGMRES.H>

using namespace amrex;

MLGMRESSolver::MLGMRESSolver() {
   m_verbose = 0;
   m_restart = 1;
}

MLGMRESSolver::~MLGMRESSolver() {}

void
MLGMRESSolver::define(PeleLM* a_pelelm,
                    const int a_KrylovSize, 
                    const int a_nComp,
                    const int a_nGrow)
{
   BL_PROFILE("MLGMRESSolver::define()");

   m_pelelm = a_pelelm;

   int finest_level = m_pelelm->finestLevel();

   // Resize level vector
   m_grids.resize(finest_level+1);
   m_dmap.resize(finest_level+1);
   Ax.resize(finest_level+1);
   res.resize(finest_level+1);

   // Grab stuff from AmrCore
   m_geom = m_pelelm->Geom(0,finest_level);
   m_grids = m_pelelm->boxArray(0,finest_level);
   m_dmap = m_pelelm->DistributionMap(0,finest_level);

   // Internal GMRES size
   m_nComp = a_nComp;
   m_nGrow = a_nGrow;
   m_krylovSize = a_KrylovSize;

   // Build krylov base memory and work MFs
   KspBase.resize(m_krylovSize+1);
   for (int n = 0; n <= m_krylovSize ; ++n) {
      KspBase[n].resize(finest_level+1);
      for (int lev = 0; lev <= finest_level; lev++) {
         KspBase[n][lev].define(m_grids[lev],m_dmap[lev],m_nComp,m_nGrow);
      }
   }
   for (int lev = 0; lev <= finest_level; lev++) {
      Ax[lev].define(m_grids[lev],m_dmap[lev],m_nComp,m_nGrow);
      res[lev].define(m_grids[lev],m_dmap[lev],m_nComp,m_nGrow);
   }

   // Work Reals
   H.resize(m_krylovSize+1);
   givens.resize(m_krylovSize+1);
   for (int n = 0; n <= m_krylovSize ; ++n) {
      H[n].resize(m_krylovSize,0.0);
      givens[n].resize(2,0.0);
   }
   y.resize(m_krylovSize+1,0.0);
   g.resize(m_krylovSize+1,0.0);
}

MLJtimesVFunc
MLGMRESSolver::jtimesv() const noexcept
{
   return m_jtv;
}

MLPrecondFunc
MLGMRESSolver::precond() const noexcept
{
   return m_prec;
}

MLNormFunc
MLGMRESSolver::norm() const noexcept
{
   return m_norm;
}

void
MLGMRESSolver::setJtimesV(MLJtimesVFunc a_jtv)
{
   m_jtv = a_jtv;
}

void
MLGMRESSolver::setPrecond(MLPrecondFunc a_prec)
{
   m_prec = a_prec;
}

void
MLGMRESSolver::setNorm(MLNormFunc a_norm)
{
   m_norm = a_norm;
}

int
MLGMRESSolver::solve(const Vector<MultiFab*> &a_sol,
                     const Vector<MultiFab*> &a_rhs,
                     amrex::Real a_abs_tol,
                     amrex::Real a_rel_tol)
{
   BL_PROFILE("MLGMRESSolver::solve()");

// Checks
   AMREX_ALWAYS_ASSERT(m_jtv != nullptr);
   AMREX_ALWAYS_ASSERT(a_sol[0]->nComp() == m_nComp && a_rhs[0]->nComp() == m_nComp);
   if ( m_prec == nullptr ) Print() << "Using unpreconditioned GMRES ! Might take a while to converge ...\n";


   initResNorm = computeMLResidualNorm(a_sol,a_rhs);                      // Initial resisual norm
   Real rhsNorm  = computeMLNorm(a_rhs);                                  // RHS norm
   targetResNorm = initResNorm * a_rel_tol;                             // Target relative tolerance

   if ( m_verbose > 0 ) {
      amrex::Print() << "GMRES: Initial rhs = " << rhsNorm << "\n";
      amrex::Print() << "GMRES: Initial residual = " << initResNorm << "\n";
   }
   if ( initResNorm < a_abs_tol ) {
      amrex::Print() << "GMRES: no need for iterations \n";
      return 0;
   }

   iter_count = 0;
   /*
   restart_count = 0;
   m_converged = false;
   do {
//    Prepare for solve
      prepareForSolve();
      one_restart(a_sol,a_rhs);
      restart_count++;
   } while( !m_converged && restart_count < m_restart );

   Real finalResNorm = computeResidualNorm(a_sol,a_rhs);                // Final resisual norm
   if ( m_verbose > 0 ) amrex::Print() << "GMRES: Final residual, resid/resid0 = " << finalResNorm << ", "
                                       << finalResNorm/initResNorm << "\n";
   */
   return iter_count;
}

void
MLGMRESSolver::prepareForSolve()
{
   // Initilize or cleanup GMRES data
   for ( int k = 0 ; k <= m_krylovSize ; ++k ) {
      g[k] = 0.0;
      y[k] = 0.0;
      for ( int n = 0 ; n < m_krylovSize ; ++n ) {
         H[k][n] = 0.0;
      }
      givens[k][0] = 0.0;
      givens[k][1] = 0.0;
   }
}

void
MLGMRESSolver::one_restart(const Vector<MultiFab*> &a_x, const Vector<MultiFab*> &a_rhs)
{
   computeMLResidual(a_x,a_rhs,GetVecOfPtrs(res));
   Real resNorm_0 = computeMLNorm(GetVecOfPtrs(res));
   if ( m_verbose > 1 ) amrex::Print() << "     [Restart:"<< restart_count << "] initial relative res: " << resNorm_0/initResNorm << "\n";

   // Initialize KspBase with normalized residual
   g[0] = resNorm_0;
   // TODO: funct scale Vector<MultiFab> res.mult(1.0/resNorm_0);
   // TODO: funct Copy Vector<MultiFab> MultiFab::Copy(KspBase[0],res, 0 ,0, m_nComp, 0);

   int k_end;
   Real resNorm = resNorm_0;

   for ( int k = 0 ; k < m_krylovSize ; ++k ) {
      // Do one GMRES iteration, update the residual norm
      one_iter(k,resNorm);
      iter_count++;

      // Test exit condition
      if ( resNorm < targetResNorm ) {
         m_converged = true;
         k_end = k;
         break;
      }

      // Last iteration and not converged yet, just set k_end
      if ( k == m_krylovSize - 1 ) k_end = m_krylovSize - 1;
   }

   // Solve the minimization problem H.y = g
   y[k_end] = g[k_end]/H[k_end][k_end];
   for ( int k = k_end-1; k >= 0; --k ) {
      Real sum_tmp = 0.0;
      for ( int j = k+1; j <= k_end; ++j ) {
         sum_tmp += H[k][j] * y[j];
      }
      y[k] = ( g[k] - sum_tmp ) / H[k][k];
   }

   // Compute solution update
   // TODO: update function
   /*
   MultiFab update(m_grids,m_dmap,m_nComp,0);
   update.setVal(0.0);
   for ( int i = k_end; i >= 0; --i ) {
      MultiFab::Saxpy(update,y[i],KspBase[i], 0, 0, m_nComp, 0);
   }

   // Update the solution MF
   MultiFab::Add(a_x,update,0,0,m_nComp,0);
   */
}

void
MLGMRESSolver::one_iter(const int iter, Real &resNorm)
{
   BL_PROFILE("MLGMRESSolver::one_iter()");
   if ( m_verbose > 1 ) amrex::Print() << "     [Iter:" << iter_count << "] residual norm: " << resNorm / initResNorm << "\n";
   /*
   appendBasisVector(iter,KspBase);
   gramSchmidtOrtho(iter,KspBase);
   */
   resNorm = givensRotation(iter);
}

void
MLGMRESSolver::appendBasisVector(const int iter, Vector<MultiFab>& Base)
{
   /*
   if ( m_prec == nullptr ) {
      MEMBER_FUNC_PTR(*m_pelelm,m_jtv)(Base[iter],Base[iter+1]);
   } else {
      MEMBER_FUNC_PTR(*m_pelelm,m_jtv)(Base[iter],Ax);
      MEMBER_FUNC_PTR(*m_pelelm,m_prec)(Ax,Base[iter+1]);
   }
   */
}

void
MLGMRESSolver::gramSchmidtOrtho(const int iter, Vector<MultiFab>& Base)
{
   /*
   for ( int row = 0; row <= iter; ++row ) {
      H[row][iter] = MultiFab::Dot(Base[iter+1],0,Base[row],0,m_nComp,0);
      Real GS_corr = - H[row][iter];
      MultiFab::Saxpy(Base[iter+1],GS_corr,Base[row],0,0,m_nComp,0);
      if ( check_GramSchmidtOrtho ) {
         Real Hcorr = MultiFab::Dot(Base[iter+1],0,Base[row],0,m_nComp,0);
         if ( std::fabs(Hcorr) > 1.0e-15 ) {
            H[row][iter] += Hcorr;
            GS_corr = - Hcorr;
            MultiFab::Saxpy(Base[iter+1],GS_corr,Base[row],0,0,m_nComp,0);
         }
      }
   }
   Real normNewVec = computeNorm(Base[iter+1]);
   H[iter+1][iter] = normNewVec;
   if ( normNewVec > 0 ) Base[iter+1].mult(1.0/normNewVec);
   */
}

Real
MLGMRESSolver::givensRotation(const int iter)
{
   for ( int row = 0; row < iter; ++row ) {
      Real v1 =   givens[row][0] * H[row][iter] - givens[row][1] * H[row+1][iter];
      Real v2 =   givens[row][1] * H[row][iter] + givens[row][0] * H[row+1][iter];
      H[row][iter] = v1;
      H[row+1][iter] = v2;
   }
   Real norm_Hlast = std::sqrt(H[iter][iter]*H[iter][iter] + H[iter+1][iter]*H[iter+1][iter]);
   if ( norm_Hlast > 0.0 ) {
      givens[iter][0] =  H[iter][iter] / norm_Hlast;
      givens[iter][1] = -H[iter+1][iter] / norm_Hlast;
      H[iter][iter] = givens[iter][0] * H[iter][iter] - givens[iter][1] * H[iter+1][iter];
      H[iter+1][iter] = 0.0;
      Real v1 = givens[iter][0] * g[iter] - givens[iter][1] * g[iter+1];
      Real v2 = givens[iter][1] * g[iter] + givens[iter][0] * g[iter+1];
      g[iter] = v1;
      g[iter+1] = v2;
   }
   return std::fabs(g[iter+1]);
}

Real
MLGMRESSolver::computeMLResidualNorm(const Vector<MultiFab*> &a_x,
                                     const Vector<MultiFab*> &a_rhs)
{
   computeMLResidual(a_x,a_rhs,GetVecOfPtrs(res));
   return computeMLNorm(GetVecOfPtrs(res));
}

Real
MLGMRESSolver::computeResidualNorm(const MultiFab& a_x, const MultiFab& a_rhs)
{
   //computeResidual(a_x,a_rhs,res);
   //return computeNorm(res);
   return 0.0;
}

Real
MLGMRESSolver::computeNorm(const MultiFab& a_vec)
{
   Real norm = 0.0;
   /*
   if ( m_norm != nullptr ) {
      MEMBER_FUNC_PTR(*m_pelelm,m_norm)(a_vec,norm);
   } else {
      for ( int comp = 0; comp < m_nComp; comp++ ) {
         norm += MultiFab::Dot(a_vec,comp,a_vec,comp,1,0);
      }
      norm = std::sqrt(norm);
   }
   */
   return norm;
}

Real
MLGMRESSolver::computeMLNorm(const Vector<MultiFab*> &a_vec)
{
   Real r = 0.0;
   if ( m_norm != nullptr ) {
      MEMBER_FUNC_PTR(*m_pelelm,m_norm)(a_vec,r);
   } else {
      for ( int comp = 0; comp < m_nComp; comp++ ) {
         Real norm = 0.0;
         for (int lev = 0; lev <= a_vec.size(); ++lev) {
            norm += MultiFab::Dot(*a_vec[lev],comp,*a_vec[lev],comp,1,0);
         }
         r += norm;
      }
      r = std::sqrt(r);
   }
   return r;
}

void
MLGMRESSolver::computeResidual(const MultiFab& a_x, const MultiFab& a_rhs, MultiFab& a_res)
{
   /*
   BL_PROFILE("MLGMRESSolver::computeResidual()");
   MEMBER_FUNC_PTR(*m_pelelm,m_jtv)(a_x,Ax);
   if ( m_prec != nullptr ) {
      MultiFab::Xpay(Ax, -1.0, a_rhs, 0, 0, m_nComp, 0);
      MEMBER_FUNC_PTR(*m_pelelm,m_prec)(Ax,a_res);
   } else {
      MultiFab::LinComb(a_res,1.0,Ax,0,-1.0,a_rhs,0,0,m_nComp,0);
   }
   */
}

void
MLGMRESSolver::computeMLResidual(const Vector<MultiFab*> &a_x, const Vector<MultiFab*> &a_rhs, const Vector<MultiFab*> &a_res)
{
   BL_PROFILE("MLGMRESSolver::computeMLResidual()");
   MEMBER_FUNC_PTR(*m_pelelm,m_jtv)(a_x,GetVecOfPtrs(Ax));
   if ( m_prec != nullptr ) {
      //MultiFab::Xpay(Ax, -1.0, a_rhs, 0, 0, m_nComp, 0);
      MEMBER_FUNC_PTR(*m_pelelm,m_prec)(GetVecOfPtrs(Ax),a_res);
   } else {
      //MultiFab::LinComb(a_res,1.0,Ax,0,-1.0,a_rhs,0,0,m_nComp,0);
   }
}
