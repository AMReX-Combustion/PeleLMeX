#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

#ifdef PELE_USE_EFIELD
#include <EOS_Extension.H>
#include <PeleLMEF_Constants.H>
#endif


// Find the zone of the cell to identify the inflow, initial conditions and BC
// There are 4 zones : 1 = Fuel inflow ; 2 = Obstacle ; 3 = Oxidizer inflow ; 4 = Air inflow
// Schema : [ 1 | 2 | 3 | 4 ]
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
int get_zone(const amrex::Real x, 
             amrex::GeometryData const& geomdata,
             ProbParm const& prob_parm)
{
   const amrex::Real* prob_lo = geomdata.ProbLo();
   const amrex::Real* prob_hi = geomdata.ProbHi();

   amrex::Real split12 = prob_lo[0] + prob_parm.Xf;
   amrex::Real split23 = prob_lo[0] + prob_parm.Xf + prob_parm.Xe;
   amrex::Real split34 = prob_lo[0] + prob_parm.Xf + prob_parm.Xe + prob_parm.Xc;

   int zone = 0;

   if (x < split12) {
      zone = 1;
   } else if (x < split23) {
      zone = 2;
   } else if (x < split34) {
      zone = 3;
   } else {
      zone = 4;
   }

   return zone;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void pelelm_initdata(int i, int j, int k,
                     int is_incompressible,
                     amrex::Array4<amrex::Real> const& state,
                     amrex::Array4<amrex::Real> const& aux,
                     amrex::GeometryData const& geomdata,
                     ProbParm const& prob_parm,
                     pele::physics::PMF::PmfData::DataContainer const * pmf_data)
{
    // Initialisation of some varibales for the cell, the geometry and the initial values
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                 const amrex::Real Lz = prob_hi[2] - prob_lo[2];);

    auto eos = pele::physics::PhysicsType::eos();
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};

    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    // Set the pressure
    amrex::Real P_cgs = prob_parm.P_mean*10.0; // Conversion to CGS

    // We set temp, velocity, mass fractions and electric density for all the domain; then density, enthalpy and mass.
    if (y >= prob_parm.Yin) // Zone outside the inflow for initial conditions (see geometry)
    {
        // Mass fractions
        molefrac[O2_ID] = 0.21;
        molefrac[N2_ID] = 1 - molefrac[O2_ID];
        eos.X2Y(molefrac,massfrac);

        // Temperature
        state(i,j,k,TEMP) = prob_parm.T_air;

        // Velocity
        AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                     state(i,j,k,VELY) = prob_parm.V_air;,
                     state(i,j,k,VELZ) = 0.0;);
    } else {
        int zone = get_zone(x, geomdata, prob_parm);

        if (zone == 1) {
           // Mass fractions
           massfrac[CH4_ID] = prob_parm.dilution;
           massfrac[N2_ID] = 1 - massfrac[CH4_ID];

           // Temperature
           state(i,j,k,TEMP) = prob_parm.T_fu;

           // Velocity
           AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                        state(i,j,k,VELY) = prob_parm.V_fu;,
                        state(i,j,k,VELZ) = 0.0;);

        } else { // zone = 4
           // Mass fractions
           molefrac[O2_ID] = 0.21;
           molefrac[N2_ID] = 1 - molefrac[O2_ID];
           eos.X2Y(molefrac,massfrac);

           // Temperature
           state(i,j,k,TEMP) = prob_parm.T_air;

           // Velocity
           AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                        state(i,j,k,VELY) = prob_parm.V_air;,
                        state(i,j,k,VELZ) = 0.0;);
         }
    } 

    // Ignition zone for the flame - This zone is set to the stoechimetry for the methane combustion - The temperature profile inside the blob is hyperbolic
    if (prob_parm.do_ignition == 1) {
        amrex::Real xc = prob_parm.Xf;
        amrex::Real yc = prob_parm.Yin;
        amrex::Real rad = std::sqrt((y-yc)*(y-yc)+(x-xc)*(x-xc));
        
        if (rad <= prob_parm.ign_rad)
        {
           amrex::Real thick = 0.005*prob_parm.ign_rad; //This is set so that the blob will loose temperature only in the last percent of radius
           amrex::Real eta = 0.5*(1.0 - tanh((rad-prob_parm.ign_rad)/(thick))); // in the center, T = ign_T ; at the edge, T = average of T_fu et ign_T 
           state(i,j,k,TEMP) = prob_parm.ign_T*eta + prob_parm.T_fu*(1.0-eta);

           molefrac[CH4_ID] = 0.095;
           molefrac[O2_ID] = 0.190;
           molefrac[N2_ID] = 1 - molefrac[CH4_ID] - molefrac[O2_ID];
           eos.X2Y(molefrac,massfrac);
        }
    }

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i,j,k,TEMP), massfrac, h_cgs);
    state(i,j,k,RHOH) = h_cgs * 1.0e-4 * state(i,j,k,DENSITY);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) 
    {
       state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
}


// Set the BC for inflows
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal( const amrex::Real x[AMREX_SPACEDIM],
          const int m_nAux,
          amrex::Real s_ext[NVAR],
          const int idir,
          const int sgn,
          const amrex::Real time,
          amrex::GeometryData const& geomdata,
          ProbParm const& prob_parm,
          pele::physics::PMF::PmfData::DataContainer const * pmf_data)
{
   const amrex::Real* prob_lo = geomdata.ProbLo();

   auto eos = pele::physics::PhysicsType::eos();
   amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};

   amrex::Real molefrac[NUM_SPECIES] = {0.0};
   amrex::Real massfrac[NUM_SPECIES] = {0.0};
   amrex::Real rho_cgs, P_cgs, RhoH_temp;

   for (int n = 0; n < NVAR; n++)
   {
      s_ext[n] = 0.0;
   }

   // BC Inflow (idir = 1 and sgn = 1 always as inflow BC is set at the bottom of the calculation zone)
   if ((idir == 1) && (sgn == 1)) {
       s_ext[VELX] = 0.0;
       s_ext[VELZ] = 0.0;
       P_cgs = prob_parm.P_mean * 10.0; //CGS->MKS conversion

       int zone = get_zone(x[0], geomdata, prob_parm);

       amrex::Real a = 0;
       amrex::Real c = 0;
       amrex::Real eta = 1;
       amrex::Real delta = 5e-5;

       switch (zone)
       {
          case 1: // Fuel inflow
          {
             // The speed in the fuel inflow has a parabolic profile - the fuel tube is supposed ideal and the speed is set at 0 at the limits of the tube. 
             a = (-1.5)*prob_parm.V_fu/(prob_parm.Xf*prob_parm.Xf);
             c = 1.5*prob_parm.V_fu;
             s_ext[VELY] = a*x[0]*x[0] + c;
             
             // Temperature and mass fractions are set as uniform
             s_ext[TEMP] = prob_parm.T_fu;
             massfrac[CH4_ID] = prob_parm.dilution;
             massfrac[N2_ID] = 1 - massfrac[CH4_ID];
          }
          break;

          case 2: // Obstacle
          {
             // Speed, Temperature and mass fractions are set as uniforms
             s_ext[VELY] = 0.0;
             s_ext[TEMP] = prob_parm.T_obst;
             massfrac[O2_ID] = 0.233;
             massfrac[N2_ID] = 0.767;
          }
          break;

          case 3: // Oxidizer inflow
          {
             // The speed in tha oxidizer inflow has an hyperbolic profile. This is due to the boundary layer of the extruded tube a=in the real experiment. 
             delta = 5e-5;
             eta = tanh((x[0]-prob_parm.Xf-prob_parm.Xe)/delta);
             s_ext[VELY] = eta*prob_parm.V_ox;

             // Temperature and mass fractions are set as uniform
             s_ext[TEMP] = prob_parm.T_ox;
             molefrac[O2_ID] = 0.21; 
             molefrac[N2_ID] = 1 - molefrac[O2_ID];
             eos.X2Y(molefrac,massfrac);
          }
          break;

          case 4: // Air around
          {
             // Speed, Temperature and mass fractions are set as uniforms
             s_ext[VELY] = prob_parm.V_air;
             s_ext[TEMP] = prob_parm.T_air;
             molefrac[O2_ID] = 0.21;
             molefrac[N2_ID] = 1 - molefrac[O2_ID];
             eos.X2Y(molefrac,massfrac);
          }
          break;

          default:
             break;
       }

      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

      for (int n = 0; n < NUM_SPECIES; n++) 
      {
         s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }
   } 
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);

    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
    if ( (beta_comp == 0) && (nComp == NUM_SPECIES) ) {
        if (dir == 1 && j <= domainBox.smallEnd(dir)) {
            //In the zone 2, we have an obstacle (fuel pipe). Thus, the flow is cut.
            const amrex::Real* prob_lo = geomdata.ProbLo();
            const amrex::Real* dx      = geomdata.CellSize();
            const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];

            if ((x > 0.0011) and (x < 0.0012)) // Hard coded size of the burner lip (we can't use void get_zone as Prob_parm is not called in this function
            {
                for (int n = beta_comp; n < beta_comp + nComp; n++) {
                   beta(i,j,k,n) = 0.0;
                }
            }
        }
    }
}
#endif
