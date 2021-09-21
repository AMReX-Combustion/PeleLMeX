#include <PeleLM.H>

using namespace amrex;

PeleLM::PeleLM() = default;

PeleLM::~PeleLM()
{
   for (int lev = 0; lev <= finest_level; ++lev) {
      ClearLevel(lev);
   }
   prob_parm.reset();
   trans_parms.deallocate();
   m_reactor->close();

   closeTempFile();
}

PeleLM::LevelData*
PeleLM::getLevelDataPtr(int lev, const PeleLM::TimeStamp &a_time, int useUMac)
{
   AMREX_ASSERT(a_time==AmrOldTime || a_time==AmrNewTime || a_time==AmrHalfTime);
   if ( a_time == AmrOldTime ) { 
      return m_leveldata_old[lev].get();
   } else if ( a_time == AmrNewTime ) {
      return m_leveldata_new[lev].get();
   } else {
      m_leveldata_floating.reset( new LevelData(grids[lev], dmap[lev], *m_factory[lev],
                                  m_incompressible, m_has_divu,
                                  m_nAux, m_nGrowState, m_nGrowMAC));
      Real time = getTime(lev,a_time);
      if (useUMac) {
         // TODO: find a way to get U^{n+1/2} from Umac
         // For now get old time
         Real oldtime = getTime(lev,AmrOldTime);
         fillpatch_velocity(lev, oldtime, m_leveldata_floating->velocity, m_nGrowState);
      } else {
         fillpatch_velocity(lev, time, m_leveldata_floating->velocity, m_nGrowState);
      }
      if (!m_incompressible) {
         fillpatch_density(lev, time, m_leveldata_floating->density, m_nGrowState);
         fillpatch_species(lev, time, m_leveldata_floating->species, m_nGrowState);
         fillpatch_energy(lev, time, m_leveldata_floating->rhoh, m_leveldata_floating->temp, m_nGrowState);
#ifdef PLM_USE_EFIELD
         fillpatch_phiV(lev, time, m_leveldata_floating->phiV, m_nGrowState);
         fillpatch_nE(lev, time, m_leveldata_floating->nE, m_nGrowState);
#endif
      }
      return m_leveldata_floating.get();
   }
}

PeleLM::LevelDataReact*
PeleLM::getLevelDataReactPtr(int lev)
{
   if (m_do_react) {
      return m_leveldatareact[lev].get();
   } else {
      return nullptr;
   }
}

Vector<MultiFab *>
PeleLM::getVelocityVect(const TimeStamp &a_time) {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->velocity));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->velocity));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getSpeciesVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->species));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->species));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getDensityVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->density));
      }
   } else if ( a_time == AmrNewTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->density));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         Real time = getTime(lev,a_time);
         m_halfTimeDensity[lev].reset( new MultiFab(grids[lev], dmap[lev], 1, m_nGrowState) );
         fillpatch_density(lev,time,*(m_halfTimeDensity[lev].get()),m_nGrowState);
         r.push_back(m_halfTimeDensity[lev].get());
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getTempVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->temp));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->temp));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getDiffusivityVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->diff_cc));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->diff_cc));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getViscosityVect(const TimeStamp &a_time) {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->visc_cc));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->visc_cc));
      }
   }
   return r;
}

void
PeleLM::averageDownState(const PeleLM::TimeStamp &a_time)
{
   for (int lev = finest_level; lev > 0; --lev) {
      std::unique_ptr<MultiFab> stateFine = fillPatchState(lev, getTime(lev,a_time), 0);
      std::unique_ptr<MultiFab> stateCrse = fillPatchState(lev-1, getTime(lev,a_time), 0);
#ifdef AMREX_USE_EB
      EB_average_down(*stateFine,
                      *stateCrse,
                      0,NVAR,refRatio(lev-1));
#else
      average_down(*stateFine,
                   *stateCrse,
                   0,NVAR,refRatio(lev-1));
#endif
   }
}

void
PeleLM::averageDownDensity(const PeleLM::TimeStamp &a_time)
{
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->density,
                      ldataCrse_p->density,
                      0,1,refRatio(lev-1));
#else
      average_down(ldataFine_p->density,
                   ldataCrse_p->density,
                   0,1,refRatio(lev-1));
#endif
   }
}

void
PeleLM::averageDownSpecies(const PeleLM::TimeStamp &a_time)
{
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->species,
                      ldataCrse_p->species,
                      0,NUM_SPECIES,refRatio(lev-1));
#else
      average_down(ldataFine_p->species,
                   ldataCrse_p->species,
                   0,NUM_SPECIES,refRatio(lev-1));
#endif
   }
}

void
PeleLM::averageDownVelocity(const PeleLM::TimeStamp &a_time)
{
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->velocity,
                      ldataCrse_p->velocity,
                      0,AMREX_SPACEDIM,refRatio(lev-1));
#else
      average_down(ldataFine_p->velocity,
                   ldataCrse_p->velocity,
                   0,AMREX_SPACEDIM,refRatio(lev-1));
#endif
   }
}

void
PeleLM::averageDownTemp(const PeleLM::TimeStamp &a_time)
{
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->temp,
                      ldataCrse_p->temp,
                      0,1,refRatio(lev-1));
#else
      average_down(ldataFine_p->temp,
                   ldataCrse_p->temp,
                   0,1,refRatio(lev-1));
#endif
   }
}

void
PeleLM::averageDownEnthalpy(const PeleLM::TimeStamp &a_time)
{
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->rhoh,
                      ldataCrse_p->rhoh,
                      0,1,refRatio(lev-1));
#else
      average_down(ldataFine_p->rhoh,
                   ldataCrse_p->rhoh,
                   0,1,refRatio(lev-1));
#endif
   }
}

void
PeleLM::averageDownRhoRT(const PeleLM::TimeStamp &a_time)
{
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->rhoRT,
                      ldataCrse_p->rhoRT,
                      0,1,refRatio(lev-1));
#else
      average_down(ldataFine_p->rhoRT,
                   ldataCrse_p->rhoRT,
                   0,1,refRatio(lev-1));
#endif
   }
}

/*
Array<Real,3>
PeleLM::MFStat (const Vector<const MultiFab*> &a_mf, int comp)
{
   // Get the min/max/mean of a given component, not including the fine-covered cells
}
*/

Real 
PeleLM::MFSum (const Vector<const MultiFab*> &a_mf, int comp)
{
    // Get the integral of the MF, not including the fine-covered cells

    Real  volwgtsum = 0.0;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const Real* dx = geom[lev].CellSize();

       // Use amrex::ReduceSum
       Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
#ifdef AMREX_USE_EB
       const EBFArrayBoxFactory* ebfact = &EBFactory(lev);
       auto const& vfrac = ebfact->getVolFrac();
   
       Real sm = amrex::ReduceSum(*a_mf[lev], vfrac, 0, [vol, comp]
       AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr, Array4<Real const> const& vf_arr) -> Real
       {
           Real sum = 0.0;
           AMREX_LOOP_3D(bx, i, j, k,
           {
               sum += mf_arr(i,j,k,comp) * vf_arr(i,j,k) * vol;
           });
           return sum;
       });
#else
       Real sm = 0.0;
       if ( lev != finest_level ) {
          sm = amrex::ReduceSum(*a_mf[lev], *m_coveredMask[lev], 0, [vol, comp]
          AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr, Array4<int const> const& covered_arr) -> Real
          {
              Real sum = 0.0;
              AMREX_LOOP_3D(bx, i, j, k,
              {
                  sum += mf_arr(i,j,k,comp) * vol * covered_arr(i,j,k);
              });
              return sum;
          });
       } else {
          sm = amrex::ReduceSum(*a_mf[lev], 0, [vol, comp]
          AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr) -> Real
          {
              Real sum = 0.0;
              AMREX_LOOP_3D(bx, i, j, k,
              {
                  sum += mf_arr(i,j,k,comp) * vol;
              });
              return sum;
          });
       }
#endif

        volwgtsum += sm;
    } // lev

    ParallelDescriptor::ReduceRealSum(volwgtsum);

    return volwgtsum;
}


#ifdef PLM_USE_EFIELD
void
PeleLM::averageDownnE(const PeleLM::TimeStamp &a_time)
{
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->nE,
                      ldataCrse_p->nE,
                      0,1,refRatio(lev-1));
#else
      average_down(ldataFine_p->nE,
                   ldataCrse_p->nE,
                   0,1,refRatio(lev-1));
#endif
   }
}

void
PeleLM::averageDownPhiV(const PeleLM::TimeStamp &a_time)
{
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->phiV,
                      ldataCrse_p->phiV,
                      0,1,refRatio(lev-1));
#else
      average_down(ldataFine_p->phiV,
                   ldataCrse_p->phiV,
                   0,1,refRatio(lev-1));
#endif
   }
}

Vector<MultiFab *>
PeleLM::getPhiVVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->phiV));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->phiV));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getnEVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->nE));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->nE));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getnEDiffusivityVect(const TimeStamp &a_time) {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->diffE_cc));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->diffE_cc));
      }
   }
   return r;
}
#endif
