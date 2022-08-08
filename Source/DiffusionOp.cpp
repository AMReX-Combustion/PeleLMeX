#include <PeleLM.H>
#include <DiffusionOp.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB_utils.H>
#include <AMReX_EBFArrayBox.H>
#endif

using namespace amrex;

//---------------------------------------------------------------------------------------
// Diffusion Operator

DiffusionOp::DiffusionOp (PeleLM* a_pelelm, int ncomp)
                        : m_pelelm(a_pelelm), m_ncomp(ncomp)
{
   BL_PROFILE("DiffusionOp::DiffusionOp()");
   readParameters();

   // Solve LPInfo
   LPInfo info_solve;
   info_solve.setAgglomeration(1);
   info_solve.setConsolidation(1);
   info_solve.setMetricTerm(false);
   info_solve.setMaxCoarseningLevel(m_mg_max_coarsening_level);

   // Apply LPInfo (no coarsening)
   LPInfo info_apply;
   info_apply.setMaxCoarseningLevel(0);

#ifdef AMREX_USE_EB
   // Get vector of EB Factory
   Vector<EBFArrayBoxFactory const*> ebfactVec;
   for (int lev = 0; lev <= m_pelelm->finestLevel(); ++lev) {
       ebfactVec.push_back(&(m_pelelm->EBFactory(lev)));
   }
#endif

   // Scalar apply op.
#ifdef AMREX_USE_EB
   m_scal_apply_op.reset(new MLEBABecLap(m_pelelm->Geom(0,m_pelelm->finestLevel()),
                                         m_pelelm->boxArray(0,m_pelelm->finestLevel()),
                                         m_pelelm->DistributionMap(0,m_pelelm->finestLevel()),
                                         info_apply, ebfactVec, m_ncomp));
#else
   m_scal_apply_op.reset(new MLABecLaplacian(m_pelelm->Geom(0,m_pelelm->finestLevel()),
                                             m_pelelm->boxArray(0,m_pelelm->finestLevel()),
                                             m_pelelm->DistributionMap(0,m_pelelm->finestLevel()),
                                             info_apply, {}, m_ncomp));
#endif
   m_scal_apply_op->setMaxOrder(m_mg_maxorder);

   // Scalar solve op.
#ifdef AMREX_USE_EB
   m_scal_solve_op.reset(new MLEBABecLap(m_pelelm->Geom(0,m_pelelm->finestLevel()),
                                         m_pelelm->boxArray(0,m_pelelm->finestLevel()),
                                         m_pelelm->DistributionMap(0,m_pelelm->finestLevel()),
                                         info_solve, ebfactVec, m_ncomp));
#else
   m_scal_solve_op.reset(new MLABecLaplacian(m_pelelm->Geom(0,m_pelelm->finestLevel()),
                                             m_pelelm->boxArray(0,m_pelelm->finestLevel()),
                                             m_pelelm->DistributionMap(0,m_pelelm->finestLevel()),
                                             info_solve, {}, m_ncomp));
#endif
   m_scal_solve_op->setMaxOrder(m_mg_maxorder);

   // Gradient op. : scalar/coefficient already preset
#ifdef AMREX_USE_EB
   m_gradient_op.reset(new MLEBABecLap(m_pelelm->Geom(0,m_pelelm->finestLevel()),
                                       m_pelelm->boxArray(0,m_pelelm->finestLevel()),
                                       m_pelelm->DistributionMap(0,m_pelelm->finestLevel()),
                                       info_apply, ebfactVec, m_ncomp));
#else
   m_gradient_op.reset(new MLABecLaplacian(m_pelelm->Geom(0,m_pelelm->finestLevel()),
                                           m_pelelm->boxArray(0,m_pelelm->finestLevel()),
                                           m_pelelm->DistributionMap(0,m_pelelm->finestLevel()),
                                           info_apply, {}, m_ncomp));
#endif
   m_gradient_op->setMaxOrder(m_mg_maxorder);
   m_gradient_op->setScalars(0.0,1.0);
   for (int lev = 0; lev <= m_pelelm->finestLevel(); ++lev) {
      m_gradient_op->setBCoeffs(lev,-1.0);
   }

}

void DiffusionOp::diffuse_scalar(Vector<MultiFab*> const& a_phi, int phi_comp,
                                 Vector<MultiFab const*> const& a_rhs, int rhs_comp,
                                 Vector<Array<MultiFab*,AMREX_SPACEDIM>> const& a_flux, int flux_comp,
                                 Vector<MultiFab const*> const& a_acoeff,
                                 Vector<MultiFab const*> const& a_density,
                                 Vector<MultiFab const*> const& a_bcoeff, int bcoeff_comp,
                                 Vector<BCRec> a_bcrec,
                                 int ncomp,
                                 int isPoissonSolve,
                                 Real a_dt)
{
   BL_PROFILE("DiffusionOp::diffuse_scalar()");

   //----------------------------------------------------------------
   // What are we dealing with ?
   int have_density = (a_density.empty()) ? 0 : 1;
   int have_fluxes  = (a_flux.empty()) ? 0 : 1;
   int have_acoeff  = (a_acoeff.empty()) ? 0 : 1;
   int have_bcoeff  = (a_bcoeff.empty()) ? 0 : 1;

   //----------------------------------------------------------------
   // Checks
   AMREX_ASSERT(m_ncomp == 1 || m_ncomp == ncomp);
   AMREX_ASSERT(a_phi[0]->nComp() >= phi_comp+ncomp);
   AMREX_ASSERT(a_rhs[0]->nComp() >= rhs_comp+ncomp);
   if ( have_fluxes ) {
      AMREX_ASSERT(a_flux[0][0]->nComp() >= flux_comp+ncomp);
   }
   if ( have_bcoeff ) {
      AMREX_ASSERT(a_bcoeff[0]->nComp() >= bcoeff_comp+ncomp);
      AMREX_ASSERT(a_bcrec.size() >= ncomp);
   }

   int finest_level = m_pelelm->finestLevel();

   //----------------------------------------------------------------
   // Duplicate phi_old to include rho scaling
   // include 1 ghost cell to provide levelBC
   // NOTE: this is a bit weird for species: since we already updated the density after adv.,
   // when we divide by \rho, it is inconsistent. But it only matters if it screws
   // up the ghost cell values 'cause interiors are just an initial solution for the solve.
   Vector<MultiFab> phi(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      phi[lev].define(a_phi[lev]->boxArray(),a_phi[lev]->DistributionMap(),
                      ncomp, 1, MFInfo(), a_phi[lev]->Factory());
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         const Box& gbx     = mfi.growntilebox();
         auto const& a_phi_arr = a_phi[lev]->const_array(mfi,phi_comp);
         auto const& a_rho_arr = (have_density) ? a_density[lev]->const_array(mfi)
                                                : a_phi[lev]->const_array(mfi);     // Get dummy Array4 if no density
         auto const& phi_arr   = phi[lev].array(mfi);
         amrex::ParallelFor(gbx, ncomp, [a_phi_arr,a_rho_arr,phi_arr,have_density]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            if ( have_density ) {
               phi_arr(i,j,k,n) = a_phi_arr(i,j,k,n) / a_rho_arr(i,j,k);
            } else {
               phi_arr(i,j,k,n) = a_phi_arr(i,j,k,n);
            }
         });
      }
   }

   //----------------------------------------------------------------
   // Setup solve LinearOp coefficients
   // LinOp is \alpha A \phi - \beta \nabla \cdot B \nabla \phi = rhs
   // => \alpha = 1.0, A is a_acoeff if provided, 1.0 otherwise
   // => \beta = a_dt, B face centered diffusivity bcoeff^{np1,k}

   Real alpha = (isPoissonSolve) ? 0.0 : 1.0;
   Real beta  = a_dt;
   m_scal_solve_op->setScalars(alpha,beta);
   for (int lev = 0; lev <= finest_level; ++lev) {
      if (have_acoeff) {
         m_scal_solve_op->setACoeffs(lev, *a_acoeff[lev]);
      } else {
         m_scal_solve_op->setACoeffs(lev, 1.0);
      }
   }

   //----------------------------------------------------------------
   // Solve and get fluxes on a m_ncomp component basis
   for (int comp = 0; comp < ncomp; comp+=m_ncomp) {

      // Aliases
      Vector<Array<MultiFab*,AMREX_SPACEDIM>> fluxes(finest_level+1);
      Vector<MultiFab> component;
      Vector<MultiFab> rhs;

      // Allow for component specific LinOp BC
      m_scal_solve_op->setDomainBC(m_pelelm->getDiffusionLinOpBC(Orientation::low,a_bcrec[comp]),
                                   m_pelelm->getDiffusionLinOpBC(Orientation::high,a_bcrec[comp]));

      // Set aliases and bcoeff comp
      for (int lev = 0; lev <= finest_level; ++lev) {
         if (have_fluxes) {
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++ ) {
               fluxes[lev][idim] = new MultiFab(*a_flux[lev][idim],amrex::make_alias,flux_comp+comp,m_ncomp);
            }
         }

         if (have_bcoeff) {
            int doZeroVisc = 1;
            Vector<BCRec> subBCRec = {a_bcrec.begin()+comp,a_bcrec.begin()+comp+m_ncomp};
            Array<MultiFab,AMREX_SPACEDIM> bcoeff_ec = m_pelelm->getDiffusivity(lev, bcoeff_comp+comp, m_ncomp,
                                                                                doZeroVisc, subBCRec, *a_bcoeff[lev]);
#ifdef AMREX_USE_EB
            m_scal_solve_op->setBCoeffs(lev, GetArrOfConstPtrs(bcoeff_ec), MLMG::Location::FaceCentroid);
#else
            m_scal_solve_op->setBCoeffs(lev, GetArrOfConstPtrs(bcoeff_ec));
#endif
         } else {
            m_scal_solve_op->setBCoeffs(lev, 1.0);
         }

         component.emplace_back(phi[lev],amrex::make_alias,comp,m_ncomp);
         rhs.emplace_back(*a_rhs[lev],amrex::make_alias,rhs_comp+comp,m_ncomp);
         m_scal_solve_op->setLevelBC(lev, &component[lev]);
      }

      // Setup linear solver
      MLMG mlmg(*m_scal_solve_op);

      // Maximum iterations
      mlmg.setMaxIter(m_mg_max_iter);
      mlmg.setMaxFmgIter(m_mg_max_fmg_iter);
      mlmg.setBottomMaxIter(m_mg_bottom_maxiter);

      // Verbosity
      mlmg.setVerbose(m_mg_verbose);
      mlmg.setBottomVerbose(m_mg_bottom_verbose);

      mlmg.setPreSmooth(m_num_pre_smooth);
      mlmg.setPostSmooth(m_num_post_smooth);

      // Solve
      mlmg.solve(GetVecOfPtrs(component), GetVecOfConstPtrs(rhs), m_mg_rtol, m_mg_atol);

      // Need to get the fluxes
      if ( have_fluxes ) {
#ifdef AMREX_USE_EB
         mlmg.getFluxes(fluxes, MLMG::Location::FaceCentroid);
#else
         mlmg.getFluxes(fluxes, MLMG::Location::FaceCenter);
#endif

         for (int lev = 0; lev <= finest_level; ++lev) {
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++ ) {
               delete fluxes[lev][idim];
            }
         }
      }
   }

   //----------------------------------------------------------------
   // Copy the results of the solve back into a_phi
   // Times rho{np1,kp1} if needed
   // Don't touch the ghost cells
   for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         const Box& bx     = mfi.tilebox();
         auto const& a_phi_arr = a_phi[lev]->array(mfi,phi_comp);
         auto const& a_rho_arr = (have_density) ? a_density[lev]->const_array(mfi)
                                                : a_phi[lev]->const_array(mfi);     // Get dummy Array4 if no density
         auto const& phi_arr   = phi[lev].const_array(mfi);
         amrex::ParallelFor(bx, ncomp, [a_phi_arr,a_rho_arr,phi_arr,have_density]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            if ( have_density ) {
               a_phi_arr(i,j,k,n) = phi_arr(i,j,k,n) * a_rho_arr(i,j,k);
            } else {
               a_phi_arr(i,j,k,n) = phi_arr(i,j,k,n);
            }
         });
      }
   }
}

void DiffusionOp::computeDiffLap(Vector<MultiFab*> const& a_laps, int lap_comp,
                                 Vector<MultiFab const*> const& a_phi, int phi_comp,
                                 Vector<MultiFab const*> const& a_bcoeff, int bcoeff_comp,
                                 Vector<BCRec> a_bcrec,
                                 int ncomp)
{
   BL_PROFILE("DiffusionOp::computeDiffLap()");

   //----------------------------------------------------------------
   // Checks
   AMREX_ASSERT(m_ncomp == 1 || m_ncomp == ncomp);
   AMREX_ASSERT(a_laps[0]->nComp() >= lap_comp+ncomp);
   AMREX_ASSERT(a_phi[0]->nComp() >= phi_comp+ncomp);
   AMREX_ASSERT(a_bcoeff[0]->nComp() >= bcoeff_comp+ncomp);
   AMREX_ASSERT(a_bcrec.size() >= ncomp);

   int finest_level = m_pelelm->finestLevel();

   // Copy phi with 1 ghost cell
   Vector<MultiFab> phi(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      phi[lev].define(a_phi[lev]->boxArray(),a_phi[lev]->DistributionMap(),
                      ncomp, 1, MFInfo(), a_phi[lev]->Factory());
      MultiFab::Copy(phi[lev],*a_phi[lev],phi_comp,0,ncomp,1);
   }

   // LinOp is \alpha A \phi - \beta \nabla \cdot B \nabla \phi
   // => \alpha = 0, A doesn't matter
   // => \beta = -1.0, B face centered diffusivity a_bcoeff

   // Set scalars \alpha & \beta
   Real alpha = 0.0;
   Real beta  = -1.0;
   m_scal_apply_op->setScalars(alpha, beta);

   for (int comp = 0; comp < ncomp; comp+=m_ncomp) {

      // Component based vector of data
      Vector<MultiFab> laps;
      Vector<MultiFab> component;

      for (int lev = 0; lev <= finest_level; ++lev) {
          laps.emplace_back(*a_laps[lev],amrex::make_alias,lap_comp+comp,m_ncomp);
          component.emplace_back(phi[lev],amrex::make_alias,comp,m_ncomp);
          int doZeroVisc = 0;
          Vector<BCRec> subBCRec = {a_bcrec.begin()+comp,a_bcrec.begin()+comp+m_ncomp};
          Array<MultiFab,AMREX_SPACEDIM> bcoeff_ec = m_pelelm->getDiffusivity(lev, bcoeff_comp+comp, m_ncomp,
                                                                              doZeroVisc, subBCRec, *a_bcoeff[lev]);

#ifdef AMREX_USE_EB
          m_scal_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(bcoeff_ec), MLMG::Location::FaceCentroid);
#else
          m_scal_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(bcoeff_ec));
#endif
          m_scal_apply_op->setLevelBC(lev, &component[lev]);
      }

      MLMG mlmg(*m_scal_apply_op);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(component));
   }
}

void DiffusionOp::computeDiffFluxes(Vector<Array<MultiFab*,AMREX_SPACEDIM>> const& a_flux, int flux_comp,
                                    Vector<MultiFab const*> const& a_phi, int phi_comp,
                                    Vector<MultiFab const*> const& a_density,
                                    Vector<MultiFab const*> const& a_bcoeff, int bcoeff_comp,
                                    Vector<BCRec> a_bcrec,
                                    int ncomp,
                                    Real scale,
                                    int do_avgDown)
{
   BL_PROFILE("DiffusionOp::computeDiffFluxes()");

   // TODO: how come this is not used ?
   amrex::ignore_unused(scale);

   //----------------------------------------------------------------
   // Checks
   AMREX_ASSERT(m_ncomp == 1 || m_ncomp == ncomp);
   AMREX_ASSERT(a_flux[0][0]->nComp() >= flux_comp+ncomp);
   AMREX_ASSERT(a_phi[0]->nComp() >= phi_comp+ncomp);
   AMREX_ASSERT(a_bcoeff[0]->nComp() >= bcoeff_comp+ncomp);
   AMREX_ASSERT(a_bcrec.size() >= ncomp);

   int finest_level = m_pelelm->finestLevel();

   int have_density = (a_density.empty()) ? 0 : 1;

   // Duplicate phi since it is modified by the LinOp
   // and if have_density -> divide by density
   Vector<MultiFab> phi(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      phi[lev].define(a_phi[lev]->boxArray(),a_phi[lev]->DistributionMap(),
                      ncomp, 1, MFInfo(), a_phi[lev]->Factory());
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         const Box& gbx     = mfi.growntilebox();
         auto const& a_phi_arr = a_phi[lev]->const_array(mfi,phi_comp);
         auto const& a_rho_arr = (have_density) ? a_density[lev]->const_array(mfi)
                                                : a_phi[lev]->const_array(mfi);     // Get dummy Array4 if no density
         auto const& phi_arr   = phi[lev].array(mfi);
         amrex::ParallelFor(gbx, ncomp, [a_phi_arr,a_rho_arr,phi_arr,have_density]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            if ( have_density ) {
               phi_arr(i,j,k,n) = a_phi_arr(i,j,k,n) / a_rho_arr(i,j,k);
            } else {
               phi_arr(i,j,k,n) = a_phi_arr(i,j,k,n);
            }
         });
      }
   }

   // LinOp is \alpha A \phi - \beta \nabla \cdot B \nabla \phi
   // => \alpha = 0, A doesn't matter
   // => \beta = -1.0, B face centered diffusivity a_bcoeff

   // Set scalars \alpha & \beta
   Real alpha = 0.0;
   Real beta  = -1.0;
   m_scal_apply_op->setScalars(alpha, beta);

   // Get fluxes on a m_ncomp component(s) basis
   for (int comp = 0; comp < ncomp; comp+=m_ncomp) {

      // Component based vector of data
      Vector<Array<MultiFab*,AMREX_SPACEDIM>> fluxes(finest_level+1);
      Vector<MultiFab> component;
      Vector<MultiFab> laps;

      // Allow for component specific LinOp BC
      m_scal_apply_op->setDomainBC(m_pelelm->getDiffusionLinOpBC(Orientation::low,a_bcrec[comp]),
                                   m_pelelm->getDiffusionLinOpBC(Orientation::high,a_bcrec[comp]));

      for (int lev = 0; lev <= finest_level; ++lev) {
         for (int idim = 0; idim < AMREX_SPACEDIM; idim++ ) {
            fluxes[lev][idim] = new MultiFab(*a_flux[lev][idim],amrex::make_alias,flux_comp+comp,m_ncomp);
         }
         component.emplace_back(phi[lev],amrex::make_alias,comp,m_ncomp);
         int doZeroVisc = 1;
         Vector<BCRec> subBCRec = {a_bcrec.begin()+comp,a_bcrec.begin()+comp+m_ncomp};
         Array<MultiFab,AMREX_SPACEDIM> bcoeff_ec = m_pelelm->getDiffusivity(lev, bcoeff_comp+comp, m_ncomp,
                                                                             doZeroVisc, subBCRec, *a_bcoeff[lev]);
         laps.emplace_back(a_phi[lev]->boxArray(), a_phi[lev]->DistributionMap(),
                           m_ncomp, 1, MFInfo(), a_phi[lev]->Factory());
#ifdef AMREX_USE_EB
         m_scal_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(bcoeff_ec),  MLMG::Location::FaceCentroid);
#else
         m_scal_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(bcoeff_ec));
#endif
         m_scal_apply_op->setLevelBC(lev, &component[lev]);
      }

      MLMG mlmg(*m_scal_apply_op);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(component));
#ifdef AMREX_USE_EB
      mlmg.getFluxes(fluxes, GetVecOfPtrs(component),MLMG::Location::FaceCentroid);
#else
      mlmg.getFluxes(fluxes, GetVecOfPtrs(component),MLMG::Location::FaceCenter);
#endif
      for (int lev = 0; lev <= finest_level; ++lev) {
         for (int idim = 0; idim < AMREX_SPACEDIM; idim++ ) {
            delete fluxes[lev][idim];
         }
      }
   }

   // Average down if requested
   if (do_avgDown) avgDownFluxes(a_flux,flux_comp,ncomp);
}

void
DiffusionOp::computeGradient(const Vector<Array<MultiFab*,AMREX_SPACEDIM>> &a_grad,
                             const Vector<MultiFab*> &a_laps,
                             const Vector<MultiFab const*> &a_phi,
                             const BCRec &a_bcrec,
                             int do_avgDown)
{
   BL_PROFILE("DiffusionOp::computeGradient()");

   // Do I need the Laplacian out ?
   int need_laplacian = (a_laps.empty()) ? 0 : 1;

   // Checks: one components only and 1 ghost cell at least
   AMREX_ASSERT(a_phi[0]->nComp() == 1);
   AMREX_ASSERT(a_phi[0]->nGrow() >= 1);

   int finest_level = m_pelelm->finestLevel();

   // Set domainBCs
   m_gradient_op->setDomainBC(m_pelelm->getDiffusionLinOpBC(Orientation::low,a_bcrec),
                              m_pelelm->getDiffusionLinOpBC(Orientation::high,a_bcrec));

   // Duplicate phi since it is modified by the LinOp
   // and setup level BCs
   Vector<MultiFab> phi(finest_level+1);
   Vector<MultiFab> laps;
   for (int lev = 0; lev <= finest_level; ++lev) {
      phi[lev].define(a_phi[lev]->boxArray(),a_phi[lev]->DistributionMap(),
                      1, 1, MFInfo(), a_phi[lev]->Factory());
      MultiFab::Copy(phi[lev], *a_phi[lev], 0, 0, 1, 1);
      m_gradient_op->setLevelBC(lev, &phi[lev]);
      if (need_laplacian) {
         laps.emplace_back(*a_laps[lev],amrex::make_alias,0,1);
      } else {
         laps.emplace_back(a_phi[lev]->boxArray(), a_phi[lev]->DistributionMap(),
                           1, 1, MFInfo(), a_phi[lev]->Factory());
      }
   }

   MLMG mlmg(*m_gradient_op);
   mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
#ifdef AMREX_USE_EB
   mlmg.getFluxes(a_grad, GetVecOfPtrs(phi),MLMG::Location::FaceCentroid);
#else
   mlmg.getFluxes(a_grad, GetVecOfPtrs(phi),MLMG::Location::FaceCenter);
#endif
   if (do_avgDown) avgDownFluxes(a_grad, 0, 1);
}

void
DiffusionOp::avgDownFluxes(const Vector<Array<MultiFab*,AMREX_SPACEDIM>> &a_fluxes,
                           int flux_comp,
                           int ncomp)
{

   int finest_level = m_pelelm->finestLevel();

   for (int lev = finest_level; lev > 0; --lev) {
      // Get the requested components only
      Array<MultiFab*,AMREX_SPACEDIM> flux_fine;
      Array<MultiFab*,AMREX_SPACEDIM> flux_crse;
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++ ) {
         flux_fine[idim] = new MultiFab(*a_fluxes[lev][idim],amrex::make_alias,flux_comp,ncomp);
         flux_crse[idim] = new MultiFab(*a_fluxes[lev-1][idim],amrex::make_alias,flux_comp,ncomp);
      }
#ifdef AMREX_USE_EB
      EB_average_down_faces(GetArrOfConstPtrs(flux_fine),
                            flux_crse,
                            m_pelelm->refRatio(lev-1),flux_crse[0]->nGrow());
#else
      average_down_faces(GetArrOfConstPtrs(flux_fine),
                         flux_crse,
                         m_pelelm->refRatio(lev-1),flux_crse[0]->nGrow());
#endif
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++ ) {
         delete flux_fine[idim];
         delete flux_crse[idim];
      }
   }
}

void
DiffusionOp::readParameters ()
{
   ParmParse pp("diffusion");

   pp.query("verbose", m_mg_verbose);
   pp.query("atol", m_mg_atol);
   pp.query("rtol", m_mg_rtol);
   pp.query("max_iter", m_mg_max_iter);
   pp.query("bottom_solver", m_mg_bottom_solver);
}

//---------------------------------------------------------------------------------------
// Tensor Operator

DiffusionTensorOp::DiffusionTensorOp (PeleLM* a_pelelm)
   : m_pelelm(a_pelelm)
{

   readParameters();

   int finest_level = m_pelelm->finestLevel();

   auto bcRecVel = m_pelelm->fetchBCRecArray(VELX,AMREX_SPACEDIM);

   // Solve LPInfo
   LPInfo info_solve;
   info_solve.setMaxCoarseningLevel(m_mg_max_coarsening_level);

#ifdef AMREX_USE_EB
   // Get vector of EB Factory
   Vector<EBFArrayBoxFactory const*> ebfactVec;
   for (int lev = 0; lev <= finest_level; ++lev) {
       ebfactVec.push_back(&(m_pelelm->EBFactory(lev)));
   }
#endif

#ifdef AMREX_USE_EB
   m_solve_op.reset(new MLEBTensorOp(m_pelelm->Geom(0,finest_level),
                                     m_pelelm->boxArray(0,finest_level),
                                     m_pelelm->DistributionMap(0,finest_level),
                                     info_solve,ebfactVec));
#else
   m_solve_op.reset(new MLTensorOp(m_pelelm->Geom(0,finest_level),
                                   m_pelelm->boxArray(0,finest_level),
                                   m_pelelm->DistributionMap(0,finest_level),
                                   info_solve));
#endif
   m_solve_op->setMaxOrder(m_mg_maxorder);
   m_solve_op->setDomainBC(m_pelelm->getDiffusionTensorOpBC(Orientation::low,bcRecVel),
                           m_pelelm->getDiffusionTensorOpBC(Orientation::high,bcRecVel));

   // Apply LPInfo (no coarsening)
   LPInfo info_apply;
   info_apply.setMaxCoarseningLevel(0);

#ifdef AMREX_USE_EB
   m_apply_op.reset(new MLEBTensorOp(m_pelelm->Geom(0,finest_level),
                                     m_pelelm->boxArray(0,finest_level),
                                     m_pelelm->DistributionMap(0,finest_level),
                                     info_apply,ebfactVec));
#else
   m_apply_op.reset(new MLTensorOp(m_pelelm->Geom(0,finest_level),
                                   m_pelelm->boxArray(0,finest_level),
                                   m_pelelm->DistributionMap(0,finest_level),
                                   info_apply));
#endif
   m_apply_op->setMaxOrder(m_mg_maxorder);
   m_apply_op->setDomainBC(m_pelelm->getDiffusionTensorOpBC(Orientation::low,bcRecVel),
                           m_pelelm->getDiffusionTensorOpBC(Orientation::high,bcRecVel));

}

void DiffusionTensorOp::compute_divtau (Vector<MultiFab*> const& a_divtau,
                                        Vector<MultiFab const*> const& a_vel,
                                        Vector<MultiFab const*> const& a_density,
                                        Vector<MultiFab const*> const& a_beta,
                                        const BCRec &a_bcrec,
                                        Real scale)
{
   int finest_level = m_pelelm->finestLevel();

   int have_density = (a_density.empty()) ? 0 : 1;

   // Duplicate vel since it is modified by the TensorOp
   Vector<MultiFab> vel(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      vel[lev].define(a_vel[lev]->boxArray(), a_vel[lev]->DistributionMap(),
                     AMREX_SPACEDIM, 1, MFInfo(), a_vel[lev]->Factory());
      MultiFab::Copy(vel[lev], *a_vel[lev], 0, 0, AMREX_SPACEDIM, 1);
   }

#ifdef AMREX_USE_EB
   // Need a tempory divTau to apply redistribution
   Vector<MultiFab> divtau_tmp(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      divtau_tmp[lev].define(a_divtau[lev]->boxArray(),
                             a_divtau[lev]->DistributionMap(),
                             AMREX_SPACEDIM, 2, MFInfo(),
                             a_divtau[lev]->Factory());
      divtau_tmp[lev].setVal(0.0);
   }

   m_apply_op->setScalars(0.0, -scale);
   for (int lev = 0; lev <= finest_level; ++lev) {
       if (have_density) { // alpha being zero, not sure that this does anything.
          m_apply_op->setACoeffs(lev, *a_density[lev]);
       }
       int doZeroVisc = 0;
       Array<MultiFab,AMREX_SPACEDIM> beta_ec = m_pelelm->getDiffusivity(lev, 0, 1,
                                                                         doZeroVisc, {a_bcrec}, *a_beta[lev]);
       m_apply_op->setShearViscosity(lev, GetArrOfConstPtrs(beta_ec), MLMG::Location::FaceCentroid);
       m_apply_op->setEBShearViscosity(lev, *a_beta[lev]);
       m_apply_op->setLevelBC(lev, &vel[lev]);
   }

   MLMG mlmg(*m_apply_op);
   mlmg.apply(GetVecOfPtrs(divtau_tmp), GetVecOfPtrs(vel));

   // Flux redistribute explicit diffusion fluxes into outgoing a_divtau
   for(int lev = 0; lev <= finest_level; ++lev)
   {
      amrex::single_level_redistribute( divtau_tmp[lev], *a_divtau[lev], 0, AMREX_SPACEDIM, m_pelelm->Geom(lev));
   }

#else
   m_apply_op->setScalars(0.0, -scale);
   for (int lev = 0; lev <= finest_level; ++lev) {
       if (have_density) { // alpha being zero, not sure that this does anything.
          m_apply_op->setACoeffs(lev, *a_density[lev]);
       }
       int doZeroVisc = 0;
       Array<MultiFab,AMREX_SPACEDIM> beta_ec = m_pelelm->getDiffusivity(lev, 0, 1,
                                                                         doZeroVisc, {a_bcrec}, *a_beta[lev]);
       m_apply_op->setShearViscosity(lev, GetArrOfConstPtrs(beta_ec));
       m_apply_op->setLevelBC(lev, &vel[lev]);
   }

   MLMG mlmg(*m_apply_op);
   mlmg.apply(a_divtau, GetVecOfPtrs(vel));
#endif

   if (have_density) {
      for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(*a_divtau[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            auto const& divtau_arr = a_divtau[lev]->array(mfi);
            auto const& rho_arr = a_density[lev]->const_array(mfi);
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real rhoinv = 1.0/rho_arr(i,j,k);
                AMREX_D_TERM(divtau_arr(i,j,k,0) *= rhoinv;,
                             divtau_arr(i,j,k,1) *= rhoinv;,
                             divtau_arr(i,j,k,2) *= rhoinv;);
            });
         }
      }
   }
}

void DiffusionTensorOp::diffuse_velocity (Vector<MultiFab*> const& a_vel,
                                          Vector<MultiFab const*> const& a_density,
                                          Vector<MultiFab const*> const& a_beta,
                                          const BCRec &a_bcrec,
                                          Real a_dt)
{

   const int finest_level = m_pelelm->finestLevel();

   int have_density = (a_density.empty()) ? 0 : 1;

   AMREX_ASSERT( (!m_pelelm->m_incompressible && have_density) ||
                 (m_pelelm->m_incompressible && !have_density) );

   m_solve_op->setScalars(1.0, a_dt);
   for (int lev = 0; lev <= finest_level; ++lev) {
       if ( have_density ) {
          m_solve_op->setACoeffs(lev, *a_density[lev]);
       } else {
          m_solve_op->setACoeffs(lev, m_pelelm->m_rho);
       }
       int doZeroVisc = 0;
       Array<MultiFab,AMREX_SPACEDIM> beta_ec = m_pelelm->getDiffusivity(lev, 0, 1,
                                                                         doZeroVisc, {a_bcrec}, *a_beta[lev]);
#ifdef AMREX_USE_EB
       m_solve_op->setShearViscosity(lev, GetArrOfConstPtrs(beta_ec), MLMG::Location::FaceCentroid);
       m_solve_op->setEBShearViscosity(lev, *a_beta[lev]);
#else
       m_solve_op->setShearViscosity(lev, GetArrOfConstPtrs(beta_ec));
#endif
       m_solve_op->setLevelBC(lev, a_vel[lev]);
   }

   // TODO check Why * rho ??
   Vector<MultiFab> rhs(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
       rhs[lev].define(a_vel[lev]->boxArray(),
                       a_vel[lev]->DistributionMap(), AMREX_SPACEDIM, 0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(rhs[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
           Box const& bx = mfi.tilebox();
           auto const& rhs_a = rhs[lev].array(mfi);
           auto const& vel_a = a_vel[lev]->const_array(mfi);
           auto const& rho_a = (have_density) ? a_density[lev]->const_array(mfi)
                                              : a_vel[lev]->const_array(mfi);          // Dummy unused Array4
           amrex::ParallelFor(bx, AMREX_SPACEDIM, [=,rho_incomp=m_pelelm->m_rho,
                                                     is_incomp=m_pelelm->m_incompressible]
           AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
           {
               if (is_incomp) {
                  rhs_a(i,j,k,n) = rho_incomp * vel_a(i,j,k,n);
               } else {
                  rhs_a(i,j,k,n) = rho_a(i,j,k) * vel_a(i,j,k,n);
               }
           });
       }
    }

    MLMG mlmg(*m_solve_op);

    // Maximum iterations for MultiGrid / ConjugateGradients
    mlmg.setMaxIter(m_mg_max_iter);
    mlmg.setMaxFmgIter(m_mg_max_fmg_iter);
    mlmg.setBottomMaxIter(m_mg_bottom_maxiter);

    // Verbosity for MultiGrid / ConjugateGradients
    mlmg.setVerbose(m_mg_verbose);
    mlmg.setBottomVerbose(m_mg_bottom_verbose);

    mlmg.setPreSmooth(m_num_pre_smooth);
    mlmg.setPostSmooth(m_num_post_smooth);

    mlmg.solve(a_vel, GetVecOfConstPtrs(rhs), m_mg_rtol, m_mg_atol);
}

void
DiffusionTensorOp::readParameters ()
{
   ParmParse pp("tensor_diffusion");

   pp.query("verbose", m_mg_verbose);
   pp.query("atol", m_mg_atol);
   pp.query("rtol", m_mg_rtol);
   pp.query("max_iter", m_mg_max_iter);
   pp.query("bottom_solver", m_mg_bottom_solver);
}

