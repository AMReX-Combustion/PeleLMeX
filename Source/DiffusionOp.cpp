#include <PeleLM.H>
#include <DiffusionOp.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB_utils.H>
#endif

using namespace amrex;

//---------------------------------------------------------------------------------------
// Diffusion Operator

DiffusionOp::DiffusionOp (PeleLM* a_pelelm)
                        : m_pelelm(a_pelelm)
{
   readParameters();

   // Solve LPInfo
   LPInfo info_solve;
   info_solve.setMaxCoarseningLevel(m_mg_max_coarsening_level);

   // Apply LPInfo (no coarsening)
   LPInfo info_apply;
   info_apply.setMaxCoarseningLevel(0);

   // Scalar apply op.
   m_scal_apply_op.reset(new MLABecLaplacian(m_pelelm->Geom(0,m_pelelm->finestLevel()),
                                             m_pelelm->boxArray(0,m_pelelm->finestLevel()),
                                             m_pelelm->DistributionMap(0,m_pelelm->finestLevel()),
                                             info_apply));
   m_scal_apply_op->setMaxOrder(m_mg_maxorder);

   // Gradient op. : scalar/coefficient already preset
   m_gradient_op.reset(new MLABecLaplacian(m_pelelm->Geom(0,m_pelelm->finestLevel()),
                                           m_pelelm->boxArray(0,m_pelelm->finestLevel()), 
                                           m_pelelm->DistributionMap(0,m_pelelm->finestLevel()),
                                           info_apply));
   m_gradient_op->setMaxOrder(m_mg_maxorder);
   m_gradient_op->setScalars(0.0,-1.0);
   for (int lev = 0; lev <= m_pelelm->finestLevel(); ++lev) {
      m_gradient_op->setBCoeffs(lev,1.0);
   }

}

void DiffusionOp::computeDiffLap(Vector<MultiFab*> const& a_laps,
                                 Vector<MultiFab const*> const& a_phi,
                                 Vector<MultiFab const*> const& a_density,
                                 Vector<MultiFab const*> const& a_eta)
{

   // TODO: EB version

}

void DiffusionOp::computeDiffFluxes(Vector<Array<MultiFab*,AMREX_SPACEDIM>> const& a_flux,
                                    int flux_comp,
                                    Vector<MultiFab const*> const& a_phi,
                                    int phi_comp,
                                    Vector<MultiFab const*> const& a_density,
                                    Vector<MultiFab const*> const& a_beta,
                                    int beta_comp,
                                    Vector<BCRec> a_bcrec,
                                    int ncomp,
                                    Real scale,
                                    int do_avgDown)
{
   BL_PROFILE_VAR("DiffusionOp::computeDiffFluxes()", computeDiffFluxes);

   AMREX_ASSERT(a_flux[0][0]->nComp() >= flux_comp+ncomp);
   AMREX_ASSERT(a_phi[0]->nComp() >= phi_comp+ncomp);
   AMREX_ASSERT(a_beta[0]->nComp() >= beta_comp+ncomp);

   int finest_level = m_pelelm->finestLevel();

   int have_density = (a_density.empty()) ? 0 : 1;

   // Duplicate phi since it is modified by the LinOp
   // and if have_density -> divide by density
   Vector<MultiFab> phi(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      phi[lev].define(a_phi[lev]->boxArray(),a_phi[lev]->DistributionMap(),
                      ncomp, 1, MFInfo(), a_phi[lev]->Factory());
#ifdef _OPENMP
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
   // => \beta = -1.0, B face centered diffusivity a_beta

   // Set scalars \alpha & \beta
   m_scal_apply_op->setScalars(0.0, -1.0);

   // Get fluxes on a per component basis
   for (int comp = 0; comp < ncomp; ++comp) {

      // Component based vector of data
      Vector<Array<MultiFab*,AMREX_SPACEDIM>> fluxes(finest_level+1);
      Vector<MultiFab> component;
      Vector<MultiFab> laps;

      // Allow for component specific LinOp BC
      m_scal_apply_op->setDomainBC(m_pelelm->getDiffusionLinOpBC(Orientation::low,a_bcrec[comp]),
                                   m_pelelm->getDiffusionLinOpBC(Orientation::high,a_bcrec[comp]));

      for (int lev = 0; lev <= finest_level; ++lev) {
         for (int idim = 0; idim < AMREX_SPACEDIM; idim++ ) {
            fluxes[lev][idim] = new MultiFab(*a_flux[lev][idim],amrex::make_alias,flux_comp+comp,1);
         }
         component.emplace_back(phi[lev],amrex::make_alias,comp,1);
         Array<MultiFab,AMREX_SPACEDIM> beta_ec = m_pelelm->getDiffusivity(lev, beta_comp, 1, {a_bcrec[comp]}, *a_beta[lev]);
         laps.emplace_back(a_phi[lev]->boxArray(), a_phi[lev]->DistributionMap(),
                           1, 1, MFInfo(), a_phi[lev]->Factory());

         m_scal_apply_op->setBCoeffs(lev, GetArrOfConstPtrs(beta_ec));
         m_scal_apply_op->setLevelBC(lev, &component[lev]);
      }

      MLMG mlmg(*m_scal_apply_op);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(component));
#ifdef AMREX_USE_EB
      mlmg.getFluxes(fluxes, GetVecOfPtrs(component),MLMG::Location::FaceCentroid);
#else
      mlmg.getFluxes(fluxes, GetVecOfPtrs(component),MLMG::Location::FaceCenter);
#endif
   }

   // Scale the
   scaleExtensiveFluxes(a_flux,flux_comp,ncomp,-1.0*scale);
   if (do_avgDown) avgDownFluxes(a_flux,flux_comp,ncomp);
}

void
DiffusionOp::computeGradient(const Vector<Array<MultiFab*,AMREX_SPACEDIM>> &a_grad,
                             const Vector<MultiFab const*> &a_phi,
                             const BCRec &a_bcrec,
                             int do_avgDown)
{
   BL_PROFILE_VAR("DiffusionOp::computeGradient()", computeGradient);

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
      laps.emplace_back(a_phi[lev]->boxArray(), a_phi[lev]->DistributionMap(),
                        1, 1, MFInfo(), a_phi[lev]->Factory());
   }

   MLMG mlmg(*m_gradient_op);
   mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
#ifdef AMREX_USE_EB
   mlmg.getFluxes(a_grad, GetVecOfPtrs(phi),MLMG::Location::FaceCentroid);
#else
   mlmg.getFluxes(a_grad, GetVecOfPtrs(phi),MLMG::Location::FaceCenter);
#endif
   scaleExtensiveFluxes(a_grad, 0, 1, -1.0);
   if (do_avgDown) avgDownFluxes(a_grad, 0, 1);
}

void
DiffusionOp::scaleExtensiveFluxes(const Vector<Array<MultiFab*,AMREX_SPACEDIM>> &a_fluxes,
                                  int flux_comp,
                                  int ncomp,
                                  Real fac)
{

   int finest_level = m_pelelm->finestLevel();

   for (int lev = 0; lev <= finest_level; ++lev) {

      MultiFab dummy(m_pelelm->boxArray(lev),m_pelelm->DistributionMap(lev),1,0);
      const Real* dx = m_pelelm->Geom(lev).CellSize();
#if ( AMREX_SPACEDIM == 2 )
      Real areax = dx[1];
      Real areay = dx[0];
#elif ( AMREX_SPACEDIM == 3 )
      Real areax = dx[1]*dx[2];
      Real areay = dx[0]*dx[2];
      Real areaz = dx[0]*dx[1];
#endif

#ifdef AMREX_USE_EB
      auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&> a_fluxes[lev][0]->Factory();
      Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
      areafrac  = ebfactory.getAreaFrac();
#endif

      MFItInfo mfi_info;
      if (Gpu::notInLaunchRegion()) mfi_info.EnableTiling().SetDynamic(true);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(dummy, mfi_info); mfi.isValid(); ++mfi)
      {
         AMREX_D_TERM(const auto& fx = a_fluxes[lev][0]->array(mfi,flux_comp);,
                      const auto& fy = a_fluxes[lev][1]->array(mfi,flux_comp);,
                      const auto& fz = a_fluxes[lev][2]->array(mfi,flux_comp););

         AMREX_D_TERM(const Box ubx = mfi.nodaltilebox(0);,
                      const Box vbx = mfi.nodaltilebox(1);,
                      const Box wbx = mfi.nodaltilebox(2););

#ifdef AMREX_USE_EB
         Box bx = mfi.tilebox();

         const EBFArrayBox&      cc_fab = static_cast<EBFArrayBox const&>(Soln[mfi]);
         const EBCellFlagFab&    flags = cc_fab.getEBCellFlagFab();

         if ( flags.getType(amrex::grow(bx,0)) == FabType::covered ) {
            //
            // For now, set to very large num so we know if you accidentally use it
            // MLMG will set covered fluxes to zero
            //
            AMREX_D_TERM(AMREX_PARALLEL_FOR_4D(ubx, i, j, k, n, {fx(i,j,k,n) = COVERED_VAL;});,
                         AMREX_PARALLEL_FOR_4D(vbx, i, j, k, n, {fy(i,j,k,n) = COVERED_VAL;});,
                         AMREX_PARALLEL_FOR_4D(wbx, i, j, k, n, {fz(i,j,k,n) = COVERED_VAL;}););
         } else if ( flags.getType(amrex::grow(bx,0)) != FabType::regular ) {
            AMREX_D_TERM( const auto& afrac_x = areafrac[0]->array(mfi);,
                          const auto& afrac_y = areafrac[1]->array(mfi);,
                          const auto& afrac_z = areafrac[2]->array(mfi););

            AMREX_D_TERM(AMREX_PARALLEL_FOR_4D(ubx, ncomp, i, j, k, n, {fx(i,j,k,n) *= fac*areax*afrac_x(i,j,k);});,
                         AMREX_PARALLEL_FOR_4D(vbx, ncomp, i, j, k, n, {fy(i,j,k,n) *= fac*areay*afrac_y(i,j,k);});,
                         AMREX_PARALLEL_FOR_4D(wbx, ncomp, i, j, k, n, {fz(i,j,k,n) *= fac*areaz*afrac_z(i,j,k);}););
         } else
#endif
         {
            AMREX_D_TERM(AMREX_PARALLEL_FOR_4D(ubx, ncomp, i, j, k, n, {fx(i,j,k,n) *= fac*areax;});,
                         AMREX_PARALLEL_FOR_4D(vbx, ncomp, i, j, k, n, {fy(i,j,k,n) *= fac*areay;});,
                         AMREX_PARALLEL_FOR_4D(wbx, ncomp, i, j, k, n, {fz(i,j,k,n) *= fac*areaz;}););
         }
      }
   }
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
   }
}

void
DiffusionOp::readParameters ()
{
   ParmParse pp("diffusion");

   pp.query("verbose", m_verbose);
   // TODO: add all the user-defined options
}

//---------------------------------------------------------------------------------------
// Tensor Operator

DiffusionTensorOp::DiffusionTensorOp (PeleLM* a_pelelm)
   : m_pelelm(a_pelelm)
{

   readParameters();
   
   int finest_level = m_pelelm->finestLevel();

   // Solve LPInfo
   LPInfo info_solve;
   info_solve.setMaxCoarseningLevel(m_mg_max_coarsening_level);

   // Apply LPInfo (no coarsening)
   LPInfo info_apply;
   info_apply.setMaxCoarseningLevel(0);

   m_apply_op.reset(new MLTensorOp(m_pelelm->Geom(0,finest_level),
                                   m_pelelm->boxArray(0,finest_level),
                                   m_pelelm->DistributionMap(0,finest_level),
                                   info_apply));
   m_apply_op->setMaxOrder(m_mg_maxorder);
   auto bcRecVel = m_pelelm->fetchBCRecArray(VELX,AMREX_SPACEDIM);  
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

   //TODO EB

   m_apply_op->setScalars(0.0, -scale);
   for (int lev = 0; lev <= finest_level; ++lev) {
       if (have_density) { // alpha being zero, not sure that this does anything.
          m_apply_op->setACoeffs(lev, *a_density[lev]);  
       }
       Array<MultiFab,AMREX_SPACEDIM> beta_ec = m_pelelm->getDiffusivity(lev, 0, 1, {a_bcrec}, *a_beta[lev]);
       m_apply_op->setShearViscosity(lev, GetArrOfConstPtrs(beta_ec));
       m_apply_op->setLevelBC(lev, &vel[lev]);
   }

   MLMG mlmg(*m_apply_op);
   mlmg.apply(a_divtau, GetVecOfPtrs(vel));

   if (have_density) {
      for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef _OPENMP
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



void
DiffusionTensorOp::readParameters ()
{
   ParmParse pp("tensor_diffusion");

   pp.query("verbose", m_verbose);
   // TODO: add all the user-defined options
}

