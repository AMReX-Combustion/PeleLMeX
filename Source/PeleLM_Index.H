#ifndef PELELM_INDEX_H_
#define PELELM_INDEX_H_

#define VELX 0  
#define VELY 1  
#define VELZ 2  
#define DENSITY AMREX_SPACEDIM
#define FIRSTSPEC DENSITY + 1
#define RHOH FIRSTSPEC + NUM_SPECIES
#define TEMP RHOH + 1
#define FIRSTAUX TEMP + 1 
#define NVAR FIRSTAUX + m_nAux

#endif
