#include <PeleLMeX.H>
#include <memory>

#ifdef PELE_USE_SPRAY
#include "SprayParticles.H"
#endif

pele::physics::PeleParams<pele::physics::transport::TransParm<
  pele::physics::PhysicsType::eos_type,
  pele::physics::PhysicsType::transport_type>>
  PeleLM::trans_parms;

pele::physics::PeleParams<
  pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>>
  PeleLM::eos_parms;

PeleLM::PeleLM() = default;

PeleLM::~PeleLM()
{
  if (m_incompressible == 0) {
    trans_parms.deallocate();
    eos_parms.deallocate();
    m_reactor->close();
  }

  closeTempFile();
  typical_values.clear();
  freeProbParm();
  delete prob_parm;
  amrex::The_Arena()->free(prob_parm_d);
  m_initial_ba.clear();
  m_regrid_ba.clear();
#ifdef PELE_USE_SPRAY
  SprayParticleContainer::SprayCleanUp();
  SprayPC.reset();
#endif
#ifdef PELE_USE_SOOT
  cleanupSootModel();
#endif
}

PeleLM::LevelData*
PeleLM::getLevelDataPtr(
  const int lev, const TimeStamp a_time, const int /*useUMac*/)
{
  AMREX_ASSERT(
    a_time == AmrOldTime || a_time == AmrNewTime || a_time == AmrHalfTime);
  if (a_time == AmrOldTime) {
    return m_leveldata_old[lev].get();
  }
  if (a_time == AmrNewTime) {
    return m_leveldata_new[lev].get();
  }
  m_leveldata_floating = std::make_unique<LevelData>(
    grids[lev], dmap[lev], *m_factory[lev], m_incompressible, m_has_divu,
    m_nAux, m_nGrowState, m_use_soret, static_cast<int>(m_do_les));
  const amrex::Real time = getTime(lev, a_time);
  fillpatch_state(lev, time, m_leveldata_floating->state, m_nGrowState);
  if (m_nAux > 0) {
    fillpatch_aux(lev, time, m_leveldata_floating->auxiliaries, m_nGrowState);
  }
  return m_leveldata_floating.get();
}

PeleLM::LevelDataReact*
PeleLM::getLevelDataReactPtr(const int lev)
{
  if (m_do_react != 0) {
    return m_leveldatareact[lev].get();
  }
  return nullptr;
}

amrex::Vector<std::unique_ptr<amrex::MultiFab>>
PeleLM::getStateVect(const TimeStamp a_time)
{
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    if (m_incompressible != 0) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(
          std::make_unique<amrex::MultiFab>(
            m_leveldata_old[lev]->state, amrex::make_alias, 0, AMREX_SPACEDIM));
      }
    } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(
          std::make_unique<amrex::MultiFab>(
            m_leveldata_old[lev]->state, amrex::make_alias, 0, NVAR));
      }
    }
  } else {
    if (m_incompressible != 0) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(
          std::make_unique<amrex::MultiFab>(
            m_leveldata_new[lev]->state, amrex::make_alias, 0, AMREX_SPACEDIM));
      }
    } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(
          std::make_unique<amrex::MultiFab>(
            m_leveldata_new[lev]->state, amrex::make_alias, 0, NVAR));
      }
    }
  }
  return r;
}

amrex::Vector<std::unique_ptr<amrex::MultiFab>>
PeleLM::getVelocityVect(const TimeStamp a_time)
{
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_old[lev]->state, amrex::make_alias, VELX,
          AMREX_SPACEDIM));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_new[lev]->state, amrex::make_alias, VELX,
          AMREX_SPACEDIM));
    }
  }
  return r;
}

amrex::Vector<std::unique_ptr<amrex::MultiFab>>
PeleLM::getSpeciesVect(const TimeStamp a_time)
{
  AMREX_ASSERT(m_incompressible == 0);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_old[lev]->state, amrex::make_alias, FIRSTSPEC,
          NUM_SPECIES));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_new[lev]->state, amrex::make_alias, FIRSTSPEC,
          NUM_SPECIES));
    }
  }
  return r;
}

amrex::Vector<std::unique_ptr<amrex::MultiFab>>
PeleLM::getDensityVect(const TimeStamp a_time)
{
  AMREX_ASSERT(m_incompressible == 0);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_old[lev]->state, amrex::make_alias, DENSITY, 1));
    }
  } else if (a_time == AmrNewTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_new[lev]->state, amrex::make_alias, DENSITY, 1));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      const amrex::Real time = getTime(lev, a_time);
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          grids[lev], dmap[lev], 1, m_nGrowState));
      fillpatch_density(lev, time, *(r[lev]), 0, m_nGrowState);
    }
  }
  return r;
}

amrex::Vector<std::unique_ptr<amrex::MultiFab>>
PeleLM::getTempVect(const TimeStamp a_time)
{
  AMREX_ASSERT(m_incompressible == 0);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_old[lev]->state, amrex::make_alias, TEMP, 1));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_new[lev]->state, amrex::make_alias, TEMP, 1));
    }
  }
  return r;
}

amrex::Vector<std::unique_ptr<amrex::MultiFab>>
PeleLM::getRhoHVect(const TimeStamp a_time)
{
  AMREX_ASSERT(m_incompressible == 0);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_old[lev]->state, amrex::make_alias, RHOH, 1));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_new[lev]->state, amrex::make_alias, RHOH, 1));
    }
  }
  return r;
}

amrex::Vector<amrex::MultiFab*>
PeleLM::getDivUVect(const TimeStamp a_time)
{
  AMREX_ASSERT(m_incompressible == 0);
  amrex::Vector<amrex::MultiFab*> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
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

amrex::Vector<amrex::MultiFab*>
PeleLM::getDiffusivityVect(const TimeStamp a_time)
{
  AMREX_ASSERT(m_incompressible == 0);
  amrex::Vector<amrex::MultiFab*> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
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

amrex::Vector<amrex::MultiFab*>
PeleLM::getViscosityVect(const TimeStamp a_time)
{
  amrex::Vector<amrex::MultiFab*> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
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

amrex::Vector<amrex::MultiFab*>
PeleLM::getIRVect()
{
  amrex::Vector<amrex::MultiFab*> r;
  r.reserve(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    r.push_back(&(m_leveldatareact[lev]->I_R));
  }
  return r;
}

amrex::Vector<std::unique_ptr<amrex::MultiFab>>
PeleLM::getAuxVect(const TimeStamp a_time)
{
  AMREX_ASSERT(m_nAux > 0);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_old[lev]->auxiliaries, amrex::make_alias, 0, m_nAux));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_new[lev]->auxiliaries, amrex::make_alias, 0, m_nAux));
    }
  }
  return r;
}

amrex::Vector<amrex::MultiFab*>
PeleLM::getAuxDiffusivityVect(const TimeStamp a_time)
{
  AMREX_ASSERT(m_nAux > 0);
  amrex::Vector<amrex::MultiFab*> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(&(m_leveldata_old[lev]->diff_aux_cc));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(&(m_leveldata_new[lev]->diff_aux_cc));
    }
  }
  return r;
}

void
PeleLM::averageDownState(const TimeStamp a_time)
{
  const int nCompState = (m_incompressible != 0) ? AMREX_SPACEDIM : NVAR;
  for (int lev = finest_level; lev > 0; --lev) {
    auto* ldataFine_p = getLevelDataPtr(lev, a_time);
    auto* ldataCrse_p = getLevelDataPtr(lev - 1, a_time);
#ifdef AMREX_USE_EB
    EB_average_down(
      ldataFine_p->state, ldataCrse_p->state, 0, nCompState, refRatio(lev - 1));
#else
    average_down(
      ldataFine_p->state, ldataCrse_p->state, 0, nCompState, refRatio(lev - 1));
#endif
  }
}

void
PeleLM::averageDownScalars(const TimeStamp a_time)
{
#ifdef PELE_USE_PLASMA
  constexpr int nScal = NUM_SPECIES + 5;
#else
  constexpr int nScal = NUM_SPECIES + 3;
#endif
  for (int lev = finest_level; lev > 0; --lev) {
    auto* ldataFine_p = getLevelDataPtr(lev, a_time);
    auto* ldataCrse_p = getLevelDataPtr(lev - 1, a_time);
#ifdef AMREX_USE_EB
    EB_average_down(
      ldataFine_p->state, ldataCrse_p->state, DENSITY, nScal,
      refRatio(lev - 1));
#else
    average_down(
      ldataFine_p->state, ldataCrse_p->state, DENSITY, nScal,
      refRatio(lev - 1));
#endif
  }
}

void
PeleLM::averageDownAux(const TimeStamp a_time)
{
  for (int lev = finest_level; lev > 0; --lev) {
    auto* ldataFine_p = getLevelDataPtr(lev, a_time);
    auto* ldataCrse_p = getLevelDataPtr(lev - 1, a_time);
#ifdef AMREX_USE_EB
    EB_average_down(
      ldataFine_p->auxiliaries, ldataCrse_p->auxiliaries, 0, m_nAux,
      refRatio(lev - 1));
#else
    average_down(
      ldataFine_p->auxiliaries, ldataCrse_p->auxiliaries, 0, m_nAux,
      refRatio(lev - 1));
#endif
  }
}

void
PeleLM::averageDown(
  const TimeStamp a_time, const int state_comp, const int ncomp)
{
  for (int lev = finest_level; lev > 0; --lev) {
    auto* ldataFine_p = getLevelDataPtr(lev, a_time);
    auto* ldataCrse_p = getLevelDataPtr(lev - 1, a_time);
#ifdef AMREX_USE_EB
    EB_average_down(
      ldataFine_p->state, ldataCrse_p->state, state_comp, ncomp,
      refRatio(lev - 1));
#else
    average_down(
      ldataFine_p->state, ldataCrse_p->state, state_comp, ncomp,
      refRatio(lev - 1));
#endif
  }
}

void
PeleLM::averageDownVelocity(const TimeStamp a_time)
{
  for (int lev = finest_level; lev > 0; --lev) {
    auto* ldataFine_p = getLevelDataPtr(lev, a_time);
    auto* ldataCrse_p = getLevelDataPtr(lev - 1, a_time);
#ifdef AMREX_USE_EB
    EB_average_down(
      ldataFine_p->state, ldataCrse_p->state, VELX, AMREX_SPACEDIM,
      refRatio(lev - 1));
#else
    average_down(
      ldataFine_p->state, ldataCrse_p->state, VELX, AMREX_SPACEDIM,
      refRatio(lev - 1));
#endif
  }
}

void
PeleLM::averageDownReaction()
{
  for (int lev = finest_level; lev > 0; --lev) {
    auto* ldataRFine_p = getLevelDataReactPtr(lev);
    auto* ldataRCrse_p = getLevelDataReactPtr(lev - 1);
#ifdef AMREX_USE_EB
    EB_average_down(
      ldataRFine_p->I_R, ldataRCrse_p->I_R, 0, nCompIR(), refRatio(lev - 1));
#else
    average_down(
      ldataRFine_p->I_R, ldataRCrse_p->I_R, 0, nCompIR(), refRatio(lev - 1));
#endif
  }
}

#ifdef PELE_USE_PLASMA
amrex::Vector<std::unique_ptr<amrex::MultiFab>>
PeleLM::getPhiVVect(const TimeStamp a_time)
{
  AMREX_ASSERT(m_incompressible == 0);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_old[lev]->state, amrex::make_alias, PHIV, 1));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_new[lev]->state, amrex::make_alias, PHIV, 1));
    }
  }
  return r;
}

amrex::Vector<std::unique_ptr<amrex::MultiFab>>
PeleLM::getnEVect(const TimeStamp a_time)
{
  AMREX_ASSERT(m_incompressible == 0);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_old[lev]->state, amrex::make_alias, NE, 1));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(
        std::make_unique<amrex::MultiFab>(
          m_leveldata_new[lev]->state, amrex::make_alias, NE, 1));
    }
  }
  return r;
}

amrex::Vector<amrex::MultiFab*>
PeleLM::getnEDiffusivityVect(const TimeStamp a_time)
{
  amrex::Vector<amrex::MultiFab*> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
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
