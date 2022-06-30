#include <PeleLM.H>

using namespace amrex;

pele::physics::transport::TransportParams<
  pele::physics::PhysicsType::transport_type>
  PeleLM::trans_parms;

PeleLM::PeleLM() = default;

PeleLM::~PeleLM()
{
   for (int lev = 0; lev <= finest_level; ++lev) {
      ClearLevel(lev);
   }

   if (!m_incompressible) {
      trans_parms.deallocate();
      m_reactor->close();
   }

   closeTempFile();
   typical_values.clear();

   delete prob_parm;
   The_Arena()->free(prob_parm_d);
   m_initial_ba.clear();
   m_regrid_ba.clear();
#ifdef PELELM_USE_SOOT
   cleanupSootModel();
#endif
}

PeleLM::LevelData*
PeleLM::getLevelDataPtr(int lev, const PeleLM::TimeStamp &a_time, int /*useUMac*/)
{
   AMREX_ASSERT(a_time==AmrOldTime || a_time==AmrNewTime || a_time==AmrHalfTime);
   if ( a_time == AmrOldTime ) {
      return m_leveldata_old[lev].get();
   } else if ( a_time == AmrNewTime ) {
      return m_leveldata_new[lev].get();
   } else {
      m_leveldata_floating.reset( new LevelData(grids[lev], dmap[lev], *m_factory[lev],
                                  m_incompressible, m_has_divu,
                                  m_nAux, m_nGrowState));
      Real time = getTime(lev,a_time);
      fillpatch_state(lev, time, m_leveldata_floating->state, m_nGrowState);
      //if (useUMac) {
      //   // TODO: find a way to get U^{n+1/2} from Umac
      //   // For now get old time
      //   Real oldtime = getTime(lev,AmrOldTime);
      //   fillpatch_velocity(lev, oldtime, m_leveldata_floating->state, VELX, m_nGrowState);
      //}
#ifdef PELE_USE_EFIELD
      if (!m_incompressible) {
         fillpatch_phiV(lev, time, m_leveldata_floating->phiV, m_nGrowState);
         fillpatch_nE(lev, time, m_leveldata_floating->nE, m_nGrowState);
      }
#endif
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

Vector<std::unique_ptr<MultiFab> >
PeleLM::getStateVect(const TimeStamp &a_time) {
   Vector<std::unique_ptr<MultiFab> > r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      if (m_incompressible) {
         for (int lev = 0; lev <= finest_level; ++lev) {
            r.push_back(std::make_unique<MultiFab> (m_leveldata_old[lev]->state,amrex::make_alias,0,AMREX_SPACEDIM));
         }
      } else {
         for (int lev = 0; lev <= finest_level; ++lev) {
            r.push_back(std::make_unique<MultiFab> (m_leveldata_old[lev]->state,amrex::make_alias,0,NVAR));
         }
      }
   } else {
      if (m_incompressible) {
         for (int lev = 0; lev <= finest_level; ++lev) {
            r.push_back(std::make_unique<MultiFab> (m_leveldata_new[lev]->state,amrex::make_alias,0,AMREX_SPACEDIM));
         }
      } else {
         for (int lev = 0; lev <= finest_level; ++lev) {
            r.push_back(std::make_unique<MultiFab> (m_leveldata_new[lev]->state,amrex::make_alias,0,NVAR));
         }
      }
   }
   return r;
}

Vector<std::unique_ptr<MultiFab> >
PeleLM::getVelocityVect(const TimeStamp &a_time) {
   Vector<std::unique_ptr<MultiFab> > r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(std::make_unique<MultiFab> (m_leveldata_old[lev]->state,amrex::make_alias,VELX,AMREX_SPACEDIM));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(std::make_unique<MultiFab> (m_leveldata_new[lev]->state,amrex::make_alias,VELX,AMREX_SPACEDIM));
      }
   }
   return r;
}

Vector<std::unique_ptr<MultiFab> >
PeleLM::getSpeciesVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<std::unique_ptr<MultiFab> > r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(std::make_unique<MultiFab> (m_leveldata_old[lev]->state,amrex::make_alias,FIRSTSPEC,NUM_SPECIES));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(std::make_unique<MultiFab> (m_leveldata_new[lev]->state,amrex::make_alias,FIRSTSPEC,NUM_SPECIES));
      }
   }
   return r;
}

Vector<std::unique_ptr<MultiFab> >
PeleLM::getDensityVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<std::unique_ptr<MultiFab> > r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(std::make_unique<MultiFab> (m_leveldata_old[lev]->state,amrex::make_alias,DENSITY,1));
      }
   } else if ( a_time == AmrNewTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(std::make_unique<MultiFab> (m_leveldata_new[lev]->state,amrex::make_alias,DENSITY,1));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         Real time = getTime(lev,a_time);
         r.push_back(std::make_unique<MultiFab> (grids[lev], dmap[lev], 1, m_nGrowState) );
         fillpatch_density(lev,time,*(r[lev].get()),0,m_nGrowState);
      }
   }
   return r;
}

Vector<std::unique_ptr<MultiFab> >
PeleLM::getTempVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<std::unique_ptr<MultiFab> > r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(std::make_unique<MultiFab> (m_leveldata_old[lev]->state,amrex::make_alias,TEMP,1));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(std::make_unique<MultiFab> (m_leveldata_new[lev]->state,amrex::make_alias,TEMP,1));
      }
   }
   return r;
}

Vector<std::unique_ptr<MultiFab> >
PeleLM::getRhoHVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<std::unique_ptr<MultiFab> > r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(std::make_unique<MultiFab> (m_leveldata_old[lev]->state,amrex::make_alias,RHOH,1));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(std::make_unique<MultiFab> (m_leveldata_new[lev]->state,amrex::make_alias,RHOH,1));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getDivUVect(const TimeStamp &a_time) {
   AMREX_ASSERT(!m_incompressible);
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->divu));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->divu));
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

Vector<MultiFab *>
PeleLM::getIRVect() {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(&(m_leveldatareact[lev]->I_R));
   }
   return r;
}

void
PeleLM::averageDownState(const PeleLM::TimeStamp &a_time)
{
   int nCompState = ( m_incompressible ) ? AMREX_SPACEDIM : NVAR;
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->state,
                      ldataCrse_p->state,
                      0,nCompState,refRatio(lev-1));
#else
      average_down(ldataFine_p->state,
                   ldataCrse_p->state,
                   0,nCompState,refRatio(lev-1));
#endif
   }
}

void
PeleLM::averageDownScalars(const PeleLM::TimeStamp &a_time)
{
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->state,
                      ldataCrse_p->state,
                      DENSITY,NUM_SPECIES+3,refRatio(lev-1));
#else
      average_down(ldataFine_p->state,
                   ldataCrse_p->state,
                   DENSITY,NUM_SPECIES+3,refRatio(lev-1));
#endif
   }
}

void
PeleLM::averageDown(const PeleLM::TimeStamp &a_time,
                    const int state_comp,
                    const int ncomp)
{
   for (int lev = finest_level; lev > 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev,a_time);
      auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
      EB_average_down(ldataFine_p->state,
                      ldataCrse_p->state,
                      state_comp,ncomp,refRatio(lev-1));
#else
      average_down(ldataFine_p->state,
                   ldataCrse_p->state,
                   state_comp,ncomp,refRatio(lev-1));
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
      EB_average_down(ldataFine_p->state,
                      ldataCrse_p->state,
                      DENSITY,1,refRatio(lev-1));
#else
      average_down(ldataFine_p->state,
                   ldataCrse_p->state,
                   DENSITY,1,refRatio(lev-1));
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
      EB_average_down(ldataFine_p->state,
                      ldataCrse_p->state,
                      VELX,AMREX_SPACEDIM,refRatio(lev-1));
#else
      average_down(ldataFine_p->state,
                   ldataCrse_p->state,
                   VELX,AMREX_SPACEDIM,refRatio(lev-1));
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
      EB_average_down(ldataFine_p->state,
                      ldataCrse_p->state,
                      RHORT,1,refRatio(lev-1));
#else
      average_down(ldataFine_p->state,
                   ldataCrse_p->state,
                   RHORT,1,refRatio(lev-1));
#endif
   }
}

#ifdef PELE_USE_EFIELD
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
