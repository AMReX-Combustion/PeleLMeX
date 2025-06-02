#ifndef PELELM_INDEX_H
#define PELELM_INDEX_H

#define VELX 0
#define VELY 1
#define VELZ 2
#define DENSITY AMREX_SPACEDIM
#define FIRSTSPEC (DENSITY + 1)
#define RHOH (FIRSTSPEC + NUM_SPECIES)
#define TEMP (RHOH + 1)
#define RHORT (TEMP + 1)
#ifdef PELE_USE_PLASMA
#define NE (RHORT + 1)
#define PHIV (NE + 1)
#define NUMPLASMA 2
#else
#define NUMPLASMA 0
#endif
#ifdef PELE_USE_SOOT
#define FIRSTSOOT (RHORT + NUMPLASMA + 1)
#define NUMSOOTVAR (NUM_SOOT_MOMENTS + 1)
#else
#define NUMSOOTVAR 0
#endif
#ifndef NUM_ODE
#define NUM_ODE 0
#endif
#define FIRSTODE (RHORT + NUMPLASMA + NUMSOOTVAR + 1)
#define FIRSTAUX (RHORT + NUMPLASMA + NUMSOOTVAR + NUM_ODE + 1)
#define NVAR (RHORT + NUMPLASMA + NUMSOOTVAR + NUM_ODE + 1)
// TODO: handle NAUX

namespace pelelmex {
enum BCTypeEB {
  wall_adiab = 0,
  wall_isot,
  wallslip_adiab,
  wallslip_isot,
  inflow
};
}

#endif
