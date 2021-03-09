#include <PeleLM.H>
#include <PeleLMUtils.H>
#include <Godunov.H>

using namespace amrex;

void PeleLM::predictVelocity(std::unique_ptr<AdvanceAdvData>  &advData,
                        std::unique_ptr<AdvanceDiffData> &diffData)
{
   BL_PROFILE_VAR("PeleLM::predictVelocity()", predictVelocity);

   // set umac boundaries to zero 
   if ( advData->umac[0][0].nGrow() > 0 ) {
      for (int lev=0; lev <= finest_level; ++lev)
      {   
          AMREX_D_TERM(advData->umac[lev][0].setBndry(0.0);,
                       advData->umac[lev][1].setBndry(0.0);,
                       advData->umac[lev][2].setBndry(0.0););
      }   
   }

   //----------------------------------------------------------------
   // Get viscous forces
   int nGrow_force = 1;
   Vector<MultiFab> divtau(finest_level+1);
   Vector<MultiFab> velForces(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      divtau[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
      velForces[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, nGrow_force);
   }
   int use_density = 0;
   computeDivTau(AmrOldTime,GetVecOfPtrs(divtau),use_density);

   //----------------------------------------------------------------
   // Gather all the velocity forces
   // F = [ (gravity+...) - gradP + divTau ] / rho
   int add_gradP = 1;
   getVelForces(AmrOldTime,GetVecOfPtrs(divtau),GetVecOfPtrs(velForces),nGrow_force,add_gradP);

   //----------------------------------------------------------------
   // Predict face velocities with Godunov
   auto bcRecVel = fetchBCRecArray(VELX,AMREX_SPACEDIM); 
   auto bcRecVel_d = convertToDeviceVector(bcRecVel);
   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);
      Real time = getTime(lev,AmrOldTime);

#ifdef AMREX_USE_EB
#else
      godunov::predict_godunov(time,
                               AMREX_D_DECL(advData->umac[lev][0],
                                            advData->umac[lev][1],
                                            advData->umac[lev][2]),
                               ldata_p->velocity,
                               velForces[lev],
                               bcRecVel, bcRecVel_d.dataPtr(),
                               geom[lev], m_dt,
                               m_Godunov_ppm, m_Godunov_ForceInTrans);
#endif
   }
}

void PeleLM::createMACRHS()
{
   BL_PROFILE_VAR("PeleLM::createMACRHS()", createMACRHS);
}

void PeleLM::macProject(const TimeStamp &a_time,
                        std::unique_ptr<AdvanceAdvData>  &advData,
                        const Vector<MultiFab*> &a_divu)
{
   BL_PROFILE_VAR("PeleLM::macProject()", macProject);

   int has_divu = (!a_divu.empty());
   
   // Get face rho inv
   auto bcRec = fetchBCRecArray(DENSITY,1);
   Vector<Array<MultiFab,AMREX_SPACEDIM>> rho_inv(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev)
   {   
      if (m_incompressible) {
         Real rhoInv = 1.0/m_rho;
         for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            rho_inv[lev][idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim)),
                                      dmap[lev], 1, 0, MFInfo(), Factory(lev));
            rho_inv[lev][idim].setVal(rhoInv);
         }
      } else {
         auto ldata_p = getLevelDataPtr(lev,a_time);
         rho_inv[lev] = getDiffusivity(lev,0,1,{bcRec},ldata_p->density);
         for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            rho_inv[lev][idim].invert(1.0,0);
         }
      }
   }

   if (macproj->needInitialization()) {
      LPInfo lpInfo;
      lpInfo.setMaxCoarseningLevel(m_mac_mg_max_coarsening_level);
      macproj->initProjector(lpInfo, GetVecOfArrOfConstPtrs(rho_inv));  
      macproj->setDomainBC(getMACProjectionBC(Orientation::low),
                           getMACProjectionBC(Orientation::high));
   } else {
      macproj->updateBeta(GetVecOfArrOfConstPtrs(rho_inv));
   }

   // set MAC velocity and projection RHS
   macproj->setUMAC(GetVecOfArrOfPtrs(advData->umac));
   if (has_divu) macproj->setDivU(GetVecOfConstPtrs(a_divu));

   // Project
   macproj->project(m_mac_mg_rtol,m_mac_mg_atol);

   // FillBoundary umac
   for (int lev = 0; lev <= finest_level; ++lev) {
      AMREX_D_TERM(advData->umac[lev][0].FillBoundary(geom[lev].periodicity());,
                   advData->umac[lev][1].FillBoundary(geom[lev].periodicity());,
                   advData->umac[lev][2].FillBoundary(geom[lev].periodicity()));
   }
}

Array<LinOpBCType,AMREX_SPACEDIM>
PeleLM::getMACProjectionBC(Orientation::Side a_side) {

   Array<LinOpBCType,AMREX_SPACEDIM> r;
   for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      if (Geom(0).isPeriodic(idim)) {
         r[idim] = LinOpBCType::Periodic;
      } else {
         auto physbc = (a_side==Orientation::low) ? m_phys_bc.lo(idim) : m_phys_bc.hi(idim);
         if (physbc == Outflow) {
            r[idim] = LinOpBCType::Dirichlet;
         } else {
            r[idim] = LinOpBCType::Neumann;
         }
      }
   }
   return r;
}
