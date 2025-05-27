
#include "SprayInjection.H"
#include "pelelmex_prob.H"

bool
SprayParticleContainer::injectParticles(
  amrex::Real /*time*/,
  amrex::Real /*dt*/,
  int /*nstep*/,
  int /*lev*/,
  int /*finest_level*/)
{
  return false;
}

void
SprayParticleContainer::InitSprayParticles(const bool init_parts)
{
  if (!init_parts) {
    return;
  }
  const int level = 0;
  amrex::ParmParse pp("prob");
  int numRedist = -1;
  amrex::IntVect partNum(AMREX_D_DECL(100, 100, 100));
  // Find the number of redistributions during particle initialization
  pp.query("init_redist", numRedist);

  // Get particle specs
  pp.query("num_particles", partNum);
  amrex::RealVect partVel;
  std::array<amrex::Real, AMREX_SPACEDIM> pvel = {0.0};
  pp.query<amrex::Real>("part_vel", pvel);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    partVel[dir] = pvel[dir];
  }
  amrex::Real partDia, partTemp;
  pp.get("part_dia", partDia);
  pp.get("part_temp", partTemp);

  // Get composition of liquid
  std::array<amrex::Real, SPRAY_FUEL_NUM> partY = {0.0};
  if (SPRAY_FUEL_NUM > 1) {
    pp.query<amrex::Real>("Y_drop", partY);
    amrex::Real sumtest = 0.;
    for (int n = 0; n < SPRAY_FUEL_NUM; ++n) {
      sumtest += partY[n];
    }
    if (std::abs(1. - sumtest) > 1.E-6) {
      amrex::Abort("Liquid mass fractions must sum to 1!");
    }
  } else {
    partY[0] = 1.;
  }
  // Initialize particles using uniform distribution
  uniformSprayInit(
    partNum, partVel, partDia, partTemp, partY.begin(), level, numRedist);
}
