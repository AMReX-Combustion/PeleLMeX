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
  amrex::ParmParse pp_part("particles");
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
    pp_part.query<amrex::Real>("Y_0", partY);
    amrex::Real sumtest = 0.;
    amrex::Print() << "\nY_0 = ";
    for (int n = 0; n < SPRAY_FUEL_NUM; ++n) {
      sumtest += partY[n];
      amrex::Print() << partY[n] << " ";
    }
    amrex::Print() << "\n\n";
    if (std::abs(1. - sumtest) > 1.E-6) {
      amrex::Abort("Liquid mass fractions must sum to 1!");
    }
  } else {
    partY[0] = 1.;
  }

  // =========================================================================
  // MESH MAPPING: Read the ConstantMap scaling factors so particle
  // initialization is performed in physical (not index) space.  This
  // ensures consistency with the gas-phase IC when mesh mapping is active.
  // =========================================================================
  amrex::RealVect fac(AMREX_D_DECL(1.0, 1.0, 1.0));
  {
    amrex::ParmParse ppcm("ConstantMap");
    amrex::Vector<amrex::Real> fac_vec(AMREX_SPACEDIM, 1.0);
    ppcm.queryarr("scaling_factor", fac_vec, 0, AMREX_SPACEDIM);
    AMREX_D_TERM(
      fac[0] = fac_vec[0];, fac[1] = fac_vec[1];, fac[2] = fac_vec[2];);
  }

  // Scale partNum by the reciprocal of the mesh-mapping factors.
  // Under mesh mapping, the AMReX index space spans prob_lo to
  // prob_hi / fac_i.  To distribute particles uniformly across
  // the *physical* domain (prob_lo to prob_hi), we need to place them
  // at AMReX indices that account for the stretched grid.
  //   AMReX index range for direction i: [0, N_i / fac_i)
  //   To fill uniform spacing over [prob_lo_i, prob_hi_i], divide by fac_i
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    partNum[dir] =
      amrex::max(1, static_cast<int>(std::ceil(partNum[dir] / fac[dir])));
  }

  // Initialize particles using uniform distribution
  uniformSprayInit(
    partNum, partVel, partDia, partTemp, partY.begin(), level, numRedist);
}
