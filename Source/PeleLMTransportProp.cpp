#include <PeleLM.H>
#include <PeleLM_K.H>
#include <pelelm_prob.H>
#include <DiffusionOp.H>
#ifdef PELE_USE_EFIELD
#include <PeleLMEF_K.H>
#endif

using namespace amrex;

void PeleLM::calcTurbViscosity(const TimeStamp &a_time) {
   BL_PROFILE("PeleLM::calcTurbViscosity()");

   // We shouldn't be here unless we're doing LES
   AMREX_ALWAYS_ASSERT(m_do_les);

   if (m_les_verbose > 0) {
     amrex::Print() << "   Computing Turbulent Viscosity with LES model: " << m_les_model
                    << " for time " << getTime(0, a_time) << std::endl;
   }

   // Create temporary multifab to store velocity gradient tensor
   amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> GradVel(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
     auto ldata_p = getLevelDataPtr(lev,a_time);
     const auto& ba = ldata_p->state.boxArray();
     const auto& dm = ldata_p->state.DistributionMap();
     const auto& factory = ldata_p->state.Factory();
     constexpr int ncomp = AMREX_SPACEDIM*AMREX_SPACEDIM;
     for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
       GradVel[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm,
                                 ncomp, 0, MFInfo(), factory);
     }
   }

   // compute the velocity gradient
   getDiffusionTensorOp()->computeGradientTensor(GetVecOfArrOfPtrs(GradVel),
                                                 GetVecOfConstPtrs(getVelocityVect(a_time)));

   for (int lev = 0; lev <= finest_level; ++lev) {
     auto ldata_p = getLevelDataPtr(lev,a_time);

     // Get density and cp (if needed) at faces for computing turbulent transport properties
     amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> dens_fc;
     amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> cp_fc;
     amrex::MultiFab cp_cc;
     const auto& ba = ldata_p->state.boxArray();
     const auto& dm = ldata_p->state.DistributionMap();
     const auto& factory = ldata_p->state.Factory();
     if (m_incompressible) {
       // just set density as a constant; don't need to worry about cp
       for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
         dens_fc[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm,
                              1, 0, MFInfo(), factory);
         dens_fc[idim].setVal(m_rho);
       }
     } else {
       // get cp_cc (valid in 1 grow cell for interpolation to FCs)
       int ngrow = 1;
       cp_cc.define(ba, dm, 1, ngrow, MFInfo(), factory);
       auto const& state_arr      = ldata_p->state.const_arrays();
       auto const& cp_arr       = cp_cc.arrays();
       amrex::ParallelFor(cp_cc, cp_cc.nGrowVect(), [=]
                          AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                          {
                            getCpmixGivenRYT( i, j, k,
                                              Array4<Real const>(state_arr[box_no], DENSITY),
                                              Array4<Real const>(state_arr[box_no], FIRSTSPEC),
                                              Array4<Real const>(state_arr[box_no], TEMP),
                                              Array4<Real      >(cp_arr[box_no]) );
                          });
       Gpu::streamSynchronize();

       // this function really just interpolates CCs to FCs in this case
       int doZeroVisc = 0;
       auto bcRec = fetchBCRecArray(DENSITY,1);
       dens_fc = getDiffusivity(lev,DENSITY,1,doZeroVisc,{bcRec},ldata_p->state);
       cp_fc   = getDiffusivity(lev,0      ,1,doZeroVisc,{bcRec},cp_cc);
     }

     // Now we compute the turbulent viscosity at faces
     for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
       auto const& velgrad_arr = GradVel[lev][idim].const_arrays();
       auto const& dens_arr = dens_fc[idim].const_arrays();
       auto const& mut_arr = ldata_p->visc_turb_fc[idim].arrays();
       const amrex::Real vol = AMREX_D_TERM(  geom[lev].CellSize(0),
                                            * geom[lev].CellSize(1),
                                            * geom[lev].CellSize(2));
       const amrex::Real l_scale = (AMREX_SPACEDIM == 2) ? std::sqrt(vol)
                                                         : std::cbrt(vol);

#ifdef AMREX_USE_EB
       auto const& ebfact = EBFactory(lev);
       auto const vfrac = ebfact.getVolFrac().const_arrays();
#endif
       if (m_les_model == "Smagorinsky") {
         const amrex::Real prefact = m_les_cs_smag * m_les_cs_smag * l_scale * l_scale;
         amrex::ParallelFor(ldata_p->visc_turb_fc[idim], ldata_p->visc_turb_fc[idim].nGrowVect(), [=]
                            AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                            {
                              getTurbViscSmagorinsky( i, j, k, prefact,
                                                      Array4<Real const>(velgrad_arr[box_no]),
                                                      Array4<Real const>(dens_arr[box_no]),
                                                      Array4<Real      >(mut_arr[box_no]) );
#ifdef AMREX_USE_EB
                              if (idim==0) {
                                 const amrex::Real vfr_m = (AMREX_SPACEDIM==2) ? vfrac[box_no](i-1,j,k)
                                                                               : std::cbrt(vfrac[box_no](i-1,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i-1,j,k));
                                 const amrex::Real vfr_p = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k)
                                                                               : std::cbrt(vfrac[box_no](i,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k));
                                 mut_arr[box_no](i,j,k) *= amrex::min(vfr_m,vfr_p);
                              } else if (idim==1) {
                                 const amrex::Real vfr_m = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j-1,k)
                                                                               : std::cbrt(vfrac[box_no](i,j-1,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j-1,k));
                                 const amrex::Real vfr_p = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k)
                                                                               : std::cbrt(vfrac[box_no](i,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k));
                                 mut_arr[box_no](i,j,k) *= amrex::min(vfr_m,vfr_p);
                              } else if (idim==2) {
                                 const amrex::Real vfr_m = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k-1)
                                                                               : std::cbrt(vfrac[box_no](i,j,k-1)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k-1));
                                 const amrex::Real vfr_p = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k)
                                                                               : std::cbrt(vfrac[box_no](i,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k));
                                 mut_arr[box_no](i,j,k) *= amrex::min(vfr_m,vfr_p);
                              }
#endif
                            });
       } else if (m_les_model == "WALE") {
         const amrex::Real prefact = m_les_cm_wale * m_les_cm_wale * l_scale * l_scale;
         amrex::ParallelFor(ldata_p->visc_turb_fc[idim], ldata_p->visc_turb_fc[idim].nGrowVect(), [=]
                            AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                            {
                              getTurbViscWALE( i, j, k, prefact,
                                                Array4<Real const>(velgrad_arr[box_no]),
                                                Array4<Real const>(dens_arr[box_no]),
                                                Array4<Real      >(mut_arr[box_no]) );
#ifdef AMREX_USE_EB
                              if (idim==0) {
                                 const amrex::Real vfr_m = (AMREX_SPACEDIM==2) ? vfrac[box_no](i-1,j,k)
                                                                               : std::cbrt(vfrac[box_no](i-1,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i-1,j,k));
                                 const amrex::Real vfr_p = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k)
                                                                               : std::cbrt(vfrac[box_no](i,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k));
                                 mut_arr[box_no](i,j,k) *= amrex::min(vfr_m,vfr_p);
                              } else if (idim==1) {
                                 const amrex::Real vfr_m = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j-1,k)
                                                                               : std::cbrt(vfrac[box_no](i,j-1,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j-1,k));
                                 const amrex::Real vfr_p = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k)
                                                                               : std::cbrt(vfrac[box_no](i,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k));
                                 mut_arr[box_no](i,j,k) *= amrex::min(vfr_m,vfr_p);
                              } else if (idim==2) {
                                 const amrex::Real vfr_m = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k-1)
                                                                               : std::cbrt(vfrac[box_no](i,j,k-1)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k-1));
                                 const amrex::Real vfr_p = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k)
                                                                               : std::cbrt(vfrac[box_no](i,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k));
                                 mut_arr[box_no](i,j,k) *= amrex::min(vfr_m,vfr_p);
                              }
#endif
                            });
       } else if (m_les_model == "Sigma") {
         const amrex::Real prefact = m_les_cs_sigma * m_les_cs_sigma * l_scale * l_scale;
         amrex::ParallelFor(ldata_p->visc_turb_fc[idim], ldata_p->visc_turb_fc[idim].nGrowVect(), [=]
                            AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                            {
                              getTurbViscSigma( i, j, k, prefact,
                                                Array4<Real const>(velgrad_arr[box_no]),
                                                Array4<Real const>(dens_arr[box_no]),
                                                Array4<Real      >(mut_arr[box_no]) );
#ifdef AMREX_USE_EB
                              if (idim==0) {
                                 const amrex::Real vfr_m = (AMREX_SPACEDIM==2) ? vfrac[box_no](i-1,j,k)
                                                                               : std::cbrt(vfrac[box_no](i-1,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i-1,j,k));
                                 const amrex::Real vfr_p = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k)
                                                                               : std::cbrt(vfrac[box_no](i,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k));
                                 mut_arr[box_no](i,j,k) *= amrex::min(vfr_m,vfr_p);
                              } else if (idim==1) {
                                 const amrex::Real vfr_m = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j-1,k)
                                                                               : std::cbrt(vfrac[box_no](i,j-1,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j-1,k));
                                 const amrex::Real vfr_p = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k)
                                                                               : std::cbrt(vfrac[box_no](i,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k));
                                 mut_arr[box_no](i,j,k) *= amrex::min(vfr_m,vfr_p);
                              } else if (idim==2) {
                                 const amrex::Real vfr_m = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k-1)
                                                                               : std::cbrt(vfrac[box_no](i,j,k-1)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k-1));
                                 const amrex::Real vfr_p = (AMREX_SPACEDIM==2) ? vfrac[box_no](i,j,k)
                                                                               : std::cbrt(vfrac[box_no](i,j,k)) *
                                                                                 std::cbrt(vfrac[box_no](i,j,k));
                                 mut_arr[box_no](i,j,k) *= amrex::min(vfr_m,vfr_p);
                              }
#endif
                            });
       }
       Gpu::streamSynchronize();

       // Compute lambda_turb = alpha_t * cp = mu_t / Pr_t * cp
       if (!m_incompressible) {
         amrex::MultiFab::Copy(ldata_p->lambda_turb_fc[idim], ldata_p->visc_turb_fc[idim], 0, 0, 1, 0);
         amrex::MultiFab::Multiply(ldata_p->lambda_turb_fc[idim], cp_fc[idim], 0, 0, 1, 0);
         ldata_p->lambda_turb_fc[idim].mult(m_Prandtl_inv);
       }
     }
   }
}

void PeleLM::calcViscosity(const TimeStamp &a_time) {
   BL_PROFILE("PeleLM::calcViscosity()");

   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldata_p = getLevelDataPtr(lev,a_time);

      if (m_incompressible) {
         ldata_p->visc_cc.setVal(m_mu);
      } else {

         // Transport data pointer
         auto const* ltransparm = trans_parms.device_trans_parm();

         // MultiArrays
         auto const& sma = ldata_p->state.const_arrays();
         auto const& vma = ldata_p->visc_cc.arrays();

         amrex::ParallelFor(ldata_p->visc_cc, ldata_p->visc_cc.nGrowVect(), [=]
         AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
         {
            getVelViscosity( i, j, k,
                             Array4<Real const>(sma[box_no],FIRSTSPEC),
                             Array4<Real      >(sma[box_no],TEMP),
                             Array4<Real      >(vma[box_no],0),
                             ltransparm);
         });
      }
   }
   Gpu::streamSynchronize();
}

void PeleLM::calcDiffusivity(const TimeStamp &a_time) {
   BL_PROFILE("PeleLM::calcDiffusivity()");

   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldata_p = getLevelDataPtr(lev,a_time);

      // Transport data pointer
      auto const* ltransparm = trans_parms.device_trans_parm();

      // MultiArrays
      auto const& sma = ldata_p->state.const_arrays();
      auto const& dma = ldata_p->diff_cc.arrays();
#ifdef PELE_USE_EFIELD
      auto const& kma = ldata_p->mob_cc.arrays();
      GpuArray<Real,NUM_SPECIES> mwt{0.0};
      {
        auto eos = pele::physics::PhysicsType::eos();
        eos.molecular_weight(mwt.arr);
      }
#endif

      const amrex::Real Sc_inv = m_Schmidt_inv;
      const amrex::Real Pr_inv = m_Prandtl_inv;
      const int do_unity_le = m_unity_Le;
      const int do_soret = m_use_soret;
      amrex::ParallelFor(ldata_p->diff_cc, ldata_p->diff_cc.nGrowVect(), [=]
      AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
      {
         if (do_soret) {
           getTransportCoeffSoret( i, j, k,
                                  Array4<Real const>(sma[box_no],FIRSTSPEC),
                                  Array4<Real const>(sma[box_no],TEMP),
                                  Array4<Real      >(dma[box_no],0),
                                  Array4<Real      >(dma[box_no],NUM_SPECIES+2),
                                  Array4<Real      >(dma[box_no],NUM_SPECIES),
                                  Array4<Real      >(dma[box_no],NUM_SPECIES+1),
                                  ltransparm);

         } else {
           if (do_unity_le) {
             getTransportCoeffUnityLe( i, j, k, Sc_inv, Pr_inv,
                                       Array4<Real const>(sma[box_no],FIRSTSPEC),
                                       Array4<Real const>(sma[box_no],TEMP),
                                       Array4<Real      >(dma[box_no],0),
                                       Array4<Real      >(dma[box_no],NUM_SPECIES),
                                       Array4<Real      >(dma[box_no],NUM_SPECIES+1),
                                       ltransparm);
           } else {
             getTransportCoeff( i, j, k,
                                Array4<Real const>(sma[box_no],FIRSTSPEC),
                                Array4<Real const>(sma[box_no],TEMP),
                                Array4<Real      >(dma[box_no],0),
                                Array4<Real      >(dma[box_no],NUM_SPECIES),
                                Array4<Real      >(dma[box_no],NUM_SPECIES+1),
                                ltransparm);
           }
        }
#ifdef PELE_USE_EFIELD
        getKappaSp( i, j, k, mwt.arr, zk,
                    Array4<Real const>(sma[box_no],FIRSTSPEC),
                    Array4<Real      >(dma[box_no],0),
                    Array4<Real const>(sma[box_no],TEMP),
                    Array4<Real      >(kma[box_no],0));
#endif
      });
   }
   Gpu::streamSynchronize();
}

Array<MultiFab,AMREX_SPACEDIM>
PeleLM::getDiffusivity(int lev, int beta_comp, int ncomp, int doZeroVisc,
                       Vector<BCRec> bcrec,
                       MultiFab const& beta_cc,
                       int addTurbContrib)
{
   BL_PROFILE("PeleLM::getDiffusivity()");

   AMREX_ASSERT(bcrec.size() >= ncomp);
   AMREX_ASSERT(beta_cc.nComp() >= beta_comp+ncomp);

   const auto& ba = beta_cc.boxArray();
   const auto& dm = beta_cc.DistributionMap();
   const auto& factory = beta_cc.Factory();
   Array<MultiFab,AMREX_SPACEDIM> beta_ec{AMREX_D_DECL(MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                                                dm, ncomp, 0, MFInfo(), factory),
                                                       MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                                                dm, ncomp, 0, MFInfo(), factory),
                                                       MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                                                dm, ncomp, 0, MFInfo(), factory))};

   const Box& domain = geom[lev].Domain();

#ifdef AMREX_USE_EB
   // EB : use EB CCentroid -> FCentroid
   EB_interp_CellCentroid_to_FaceCentroid(beta_cc, GetArrOfPtrs(beta_ec), beta_comp, 0, ncomp, geom[lev], bcrec);
   EB_set_covered_faces(GetArrOfPtrs(beta_ec),1.234e40);
#else
   // NON-EB : use cen2edg_cpp
   bool use_harmonic_avg = m_harm_avg_cen2edge ? true : false;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(beta_cc,TilingIfNotGPU()); mfi.isValid();++mfi)
   {
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
      {
         const Box ebx = mfi.nodaltilebox(idim);
         const Box& edomain = amrex::surroundingNodes(domain,idim);
         const auto& diff_c  = beta_cc.const_array(mfi,beta_comp);
         const auto& diff_ec = beta_ec[idim].array(mfi);
         const auto bc_lo = bcrec[0].lo(idim);
         const auto bc_hi = bcrec[0].hi(idim);
         amrex::ParallelFor(ebx, [idim, ncomp, bc_lo, bc_hi, use_harmonic_avg, diff_c, diff_ec, edomain]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            int idx[3] = {i,j,k};
            bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) && ( idx[idim] <= edomain.smallEnd(idim) ) );
            bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) && ( idx[idim] >= edomain.bigEnd(idim) ) );
            cen2edg_cpp( i, j, k, idim, ncomp, use_harmonic_avg, on_lo, on_hi, diff_c, diff_ec);
         });
      }
   }
#endif

   // Add the turbulent component if desired
   // Use the following relationships to translate from the ncomp and beta_comp to the component
   // ncomp = NUM_SPECIES, beta_comp = 0           --> SPECIES DIFFUSIVITY
   // ncomp = 1          , beta_comp = NUM_SPECIES --> THERMAL CONDUCTIVITY
   // ncomp = 1          , beta_comp = 0           --> VISCOSITY
   // If PELE_USE_EFIELD is active, these relationships will not hold and LES is not supported
   if (addTurbContrib and m_do_les) {

     // If initializing the simulation, always recompute turbulent viscosity
     // otherwise, only recompute once per level per timestep (at old time)
     // calcTurbViscosity computes for all levels, so only call from the base level
     TimeStamp tstamp;
     if (getTime(lev, AmrNewTime) == 0.0 ) {
       tstamp = AmrNewTime;
       if (lev == 0) {
         calcTurbViscosity(tstamp);
       }
     } else if (lev == 0 and getTime(lev, AmrOldTime) > m_turb_visc_time[lev]) {
       tstamp = AmrOldTime;
       calcTurbViscosity(tstamp);
       m_turb_visc_time[lev] = getTime(lev, AmrOldTime);
     } else {
       tstamp = AmrOldTime;
     }
     auto ldata_p = getLevelDataPtr(lev,tstamp);

     // Identify and add the correct turbulent contribution
     for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
       if ((ncomp == 1) and (beta_comp == 0)) { // Viscosity
         amrex::MultiFab::Add(beta_ec[idim], ldata_p->visc_turb_fc[idim], 0, 0, 1, 0);
       } else if ((ncomp == NUM_SPECIES) and (beta_comp == 0)) { // Species diffusivity
         for (int ispec = 0; ispec < NUM_SPECIES; ispec++)
           amrex::MultiFab::Saxpy(beta_ec[idim], m_Schmidt_inv, ldata_p->visc_turb_fc[idim], 0, ispec, 1, 0);
       } else if ((ncomp == 1) and (beta_comp == NUM_SPECIES)) { // Thermal conductivity
         amrex::MultiFab::Add(beta_ec[idim], ldata_p->lambda_turb_fc[idim], 0, 0, 1, 0);
       } else  { // Invalid
         amrex::Abort("getDiffusivity(): LES model is on but cannot provide a turbulent transport coefficient");
       }
     }
   }

   // Enable zeroing diffusivity on faces to produce walls
   if (doZeroVisc) {
      const auto geomdata = geom[lev].data();
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         const Box& edomain = amrex::surroundingNodes(domain,idim);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(beta_ec[idim],TilingIfNotGPU()); mfi.isValid();++mfi) {
            const Box ebx = mfi.tilebox();
            const auto& diff_ec = beta_ec[idim].array(mfi);
            amrex::ParallelFor(ebx, [=]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                zero_visc(i, j, k, diff_ec, geomdata, edomain, idim, beta_comp, ncomp);
            });
         }
      }
   }

   return beta_ec;
}
