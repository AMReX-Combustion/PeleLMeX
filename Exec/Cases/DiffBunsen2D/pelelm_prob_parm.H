#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{  
   
   //Definition of the geometry (in meters)
   amrex::Real Xmax = 0.015; //Maximum size in x direction of the calcuation zone
   amrex::Real Xf = 0.0011; //Radius of the fuel pipe
   amrex::Real Xe = 0.0001; //Thickness of the extruded tube
   amrex::Real Xc = 0.0113; //Size of the oxidizer fuel
   amrex::Real Ymax = 0.037; //Maximum size in y direction of the calculation zone == Position of the electrode ine the experiment
   amrex::Real Zmax = 0.015; //Maximum size in z direction of the calculation zone : set equal to Xmax
   
   //Definition of the inflow conditions
   amrex::Real Yin = 0.004; //Height of the inflow entrance at t=0s
   
   amrex::Real V_fu = 0.15; //Speed of the fuel in inflow (m/s)
   amrex::Real V_ox = 0.031; //Speed of the oxidizer in inflow
   amrex::Real V_air = 0.0; //Speed of the air flow around the burner
   
   amrex::Real T_fu = 300; //Temperature of the fuel (in Kelvin)
   amrex::Real T_ox = 300; //Temperature of the oxidizer
   amrex::Real T_air = 300; //Temperature of the air
   amrex::Real T_obst = 450; //Temperature of the fuel pipe (=obstacle)
   
   //Definition of the initial conditions
   amrex::Real P_mean   = 100000.0_rt; //Ambient pressure in the chamber
  
   //Definition of the ignition zone
   int do_ignition = 1; //Do you want to do ignition ?
   amrex::Real ign_rad = 0.001; //Size of the ignition blob
   int ign_T = 1800; //Temperature inside the ignition blob (in Kelvin)
   
   //Definition of the fuel parameters
   amrex::Real dilution = 1.0; //Proportion of CH4 in the fuel mix. Ch4 is diluted in N2
  
   int  meanFlowDir = 1;
   amrex::Real phiV_loy = 0.0_rt;
   amrex::Real phiV_hiy = 0.0_rt;
};
#endif
