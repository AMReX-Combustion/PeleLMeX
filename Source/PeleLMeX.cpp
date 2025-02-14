#include <PeleLMeX.H>
#include <memory>

#ifdef PELE_USE_SPRAY
#include "SprayParticles.H"
#endif

using namespace amrex;

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

  delete prob_parm;
  The_Arena()->free(prob_parm_d);
  m_initial_ba.clear();
  m_regrid_ba.clear();
#ifdef PELE_USE_SPRAY
  SprayParticleContainer::SprayCleanUp();
#endif
#ifdef PELE_USE_SOOT
  cleanupSootModel();
#endif
}

PeleLM::LevelData*
PeleLM::getLevelDataPtr(
  int lev, const PeleLM::TimeStamp& a_time, int /*useUMac*/)
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
  Real time = getTime(lev, a_time);
  fillpatch_state(lev, time, m_leveldata_floating->state, m_nGrowState);
  return m_leveldata_floating.get();
}

PeleLM::LevelDataReact*
PeleLM::getLevelDataReactPtr(int lev)
{
  if (m_do_react != 0) {
    return m_leveldatareact[lev].get();
  }
  return nullptr;
}

Vector<std::unique_ptr<MultiFab>>
PeleLM::getStateVect(const TimeStamp& a_time)
{
  Vector<std::unique_ptr<MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    if (m_incompressible != 0) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(std::make_unique<MultiFab>(
          m_leveldata_old[lev]->state, amrex::make_alias, 0, AMREX_SPACEDIM));
      }
    } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(std::make_unique<MultiFab>(
          m_leveldata_old[lev]->state, amrex::make_alias, 0, NVAR));
      }
    }
  } else {
    if (m_incompressible != 0) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(std::make_unique<MultiFab>(
          m_leveldata_new[lev]->state, amrex::make_alias, 0, AMREX_SPACEDIM));
      }
    } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
        r.push_back(std::make_unique<MultiFab>(
          m_leveldata_new[lev]->state, amrex::make_alias, 0, NVAR));
      }
    }
  }
  return r;
}

Vector<std::unique_ptr<MultiFab>>
PeleLM::getVelocityVect(const TimeStamp& a_time)
{
  Vector<std::unique_ptr<MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_old[lev]->state, amrex::make_alias, VELX, AMREX_SPACEDIM));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_new[lev]->state, amrex::make_alias, VELX, AMREX_SPACEDIM));
    }
  }
  return r;
}

Vector<std::unique_ptr<MultiFab>>
PeleLM::getSpeciesVect(const TimeStamp& a_time)
{
  AMREX_ASSERT(!m_incompressible);
  Vector<std::unique_ptr<MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_old[lev]->state, amrex::make_alias, FIRSTSPEC,
        NUM_SPECIES));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_new[lev]->state, amrex::make_alias, FIRSTSPEC,
        NUM_SPECIES));
    }
  }
  return r;
}

Vector<std::unique_ptr<MultiFab>>
PeleLM::getDensityVect(const TimeStamp& a_time)
{
  AMREX_ASSERT(!m_incompressible);
  Vector<std::unique_ptr<MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_old[lev]->state, amrex::make_alias, DENSITY, 1));
    }
  } else if (a_time == AmrNewTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_new[lev]->state, amrex::make_alias, DENSITY, 1));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      Real time = getTime(lev, a_time);
      r.push_back(
        std::make_unique<MultiFab>(grids[lev], dmap[lev], 1, m_nGrowState));
      fillpatch_density(lev, time, *(r[lev]), 0, m_nGrowState);
    }
  }
  return r;
}

Vector<std::unique_ptr<MultiFab>>
PeleLM::getTempVect(const TimeStamp& a_time)
{
  AMREX_ASSERT(!m_incompressible);
  Vector<std::unique_ptr<MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_old[lev]->state, amrex::make_alias, TEMP, 1));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_new[lev]->state, amrex::make_alias, TEMP, 1));
    }
  }
  return r;
}

Vector<std::unique_ptr<MultiFab>>
PeleLM::getRhoHVect(const TimeStamp& a_time)
{
  AMREX_ASSERT(!m_incompressible);
  Vector<std::unique_ptr<MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_old[lev]->state, amrex::make_alias, RHOH, 1));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_new[lev]->state, amrex::make_alias, RHOH, 1));
    }
  }
  return r;
}

Vector<MultiFab*>
PeleLM::getDivUVect(const TimeStamp& a_time)
{
  AMREX_ASSERT(!m_incompressible);
  Vector<MultiFab*> r;
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

Vector<MultiFab*>
PeleLM::getDiffusivityVect(const TimeStamp& a_time)
{
  AMREX_ASSERT(!m_incompressible);
  Vector<MultiFab*> r;
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

Vector<MultiFab*>
PeleLM::getViscosityVect(const TimeStamp& a_time)
{
  Vector<MultiFab*> r;
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

Vector<MultiFab*>
PeleLM::getIRVect()
{
  Vector<MultiFab*> r;
  r.reserve(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    r.push_back(&(m_leveldatareact[lev]->I_R));
  }
  return r;
}

void
PeleLM::averageDownState(const PeleLM::TimeStamp& a_time)
{
  int nCompState = (m_incompressible) != 0 ? AMREX_SPACEDIM : NVAR;
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
PeleLM::averageDownScalars(const PeleLM::TimeStamp& a_time)
{
  int nScal = NUM_SPECIES + 3; // rho, rhoYs, rhoH, Temp
#ifdef PELE_USE_PLASMA
  nScal += 2; // rhoRT, nE
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
PeleLM::averageDown(
  const PeleLM::TimeStamp& a_time, const int state_comp, const int ncomp)
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
PeleLM::averageDownVelocity(const PeleLM::TimeStamp& a_time)
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
Vector<std::unique_ptr<MultiFab>>
PeleLM::getPhiVVect(const TimeStamp& a_time)
{
  AMREX_ASSERT(!m_incompressible);
  Vector<std::unique_ptr<MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_old[lev]->state, amrex::make_alias, PHIV, 1));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_new[lev]->state, amrex::make_alias, PHIV, 1));
    }
  }
  return r;
}

Vector<std::unique_ptr<MultiFab>>
PeleLM::getnEVect(const TimeStamp& a_time)
{
  AMREX_ASSERT(!m_incompressible);
  Vector<std::unique_ptr<MultiFab>> r;
  r.reserve(finest_level + 1);
  if (a_time == AmrOldTime) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_old[lev]->state, amrex::make_alias, NE, 1));
    }
  } else {
    for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(std::make_unique<MultiFab>(
        m_leveldata_new[lev]->state, amrex::make_alias, NE, 1));
    }
  }
  return r;
}

Vector<MultiFab*>
PeleLM::getnEDiffusivityVect(const TimeStamp& a_time)
{
  Vector<MultiFab*> r;
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
