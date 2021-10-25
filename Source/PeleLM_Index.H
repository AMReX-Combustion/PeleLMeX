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
#define FIRSTAUX PHIV + 1 
#else
#define FIRSTAUX RHORT + 1 
#endif
#define NVAR FIRSTAUX + 0
//TODO: handle NAUX

#endif
