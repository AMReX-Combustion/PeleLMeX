#ifndef PELELM_INDEX_H_
#define PELELM_INDEX_H_

#define VELX 0
#define VELY 1
#define VELZ 2
#define DENSITY AMREX_SPACEDIM
#define FIRSTSPEC DENSITY + 1
#define RHOH FIRSTSPEC + NUM_SPECIES
#define TEMP RHOH + 1
#define RHORT TEMP + 1
#ifdef PELE_USE_EFIELD
#define NE RHORT + 1
#define PHIV NE + 1
#define NUMEFIELD 2
#else
#define NUMEFIELD 0
#endif
#ifdef PELELM_USE_SOOT
#define FIRSTSOOT (RHORT + NUMEFIELD + 1)
#define NUMSOOTVAR (NUM_SOOT_MOMENTS + 1)
#else
#define NUMSOOTVAR 0
#endif
#define FIRSTAUX RHORT + NUMEFIELD + NUMSOOTVAR + 1
#define NVAR FIRSTAUX + 0
//TODO: handle NAUX

#endif
