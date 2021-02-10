#include <PeleLM.H>
#include <AMReX_ParmParse.H>
#include <EOS.H>
#include <Transport.H>

using namespace amrex;

void PeleLM::Setup() {
   BL_PROFILE_VAR("PeleLM::Setup()", Setup);

   // Read PeleLM parameters
   ReadParameters();

   // Setup the state variables
   VariablesSetup();

   // Initialize Level Hierarchy data
   ResizeArray();
 
   // Initialize EOS and others
   amrex::Print() << " Initialization of EOS ... \n";
   EOS::init();
   amrex::Print() << " Initialization of Transport ... \n";
   transport_init();

   // Problem parameters
   prob_parm.reset(new ProbParm{});

}

void PeleLM::ReadParameters() {
   BL_PROFILE_VAR("PeleLM::ReadParameters()", ReadParameters);

   ParmParse pp("peleLM");

   // -----------------------------------------
   // Misc
   // -----------------------------------------

   pp.query("v", m_verbose);

   // -----------------------------------------
   // Boundary conditions
   // -----------------------------------------

   Vector<std::string> lo_bc_char(AMREX_SPACEDIM);
   Vector<std::string> hi_bc_char(AMREX_SPACEDIM);
   pp.getarr("lo_bc",lo_bc_char,0,AMREX_SPACEDIM);
   pp.getarr("hi_bc",hi_bc_char,0,AMREX_SPACEDIM);

   Vector<int> lo_bc(AMREX_SPACEDIM), hi_bc(AMREX_SPACEDIM);
   for (int dir = 0; dir<AMREX_SPACEDIM; dir++){
      if (lo_bc_char[dir] == "Interior") {
         lo_bc[dir] = 0;
      } else if (lo_bc_char[dir] == "Inflow") {
         lo_bc[dir] = 1;
      } else if (lo_bc_char[dir] == "Outflow") {
         lo_bc[dir] = 2;
      } else if (lo_bc_char[dir] == "Symmetry") {
         lo_bc[dir] = 3;
      } else if (lo_bc_char[dir] == "SlipWallAdiab") {
         lo_bc[dir] = 4;
      } else if (lo_bc_char[dir] == "NoSlipWallAdiab") {
         lo_bc[dir] = 5;
      } else if (lo_bc_char[dir] == "SlipWallIsotherm") {
         lo_bc[dir] = 6;
      } else if (lo_bc_char[dir] == "NoSlipWallIsotherm") {
         lo_bc[dir] = 7;
      } else {
         amrex::Abort("Wrong boundary condition word in lo_bc, please use: Interior, Inflow, Outflow, "
                      "Symmetry, SlipWallAdiab, NoSlipWallAdiab, SlipWallIsotherm, NoSlipWallIsotherm");
      }

      if (hi_bc_char[dir] == "Interior") {
         hi_bc[dir] = 0;
      } else if (hi_bc_char[dir] == "Inflow") {
         hi_bc[dir] = 1;
      } else if (hi_bc_char[dir] == "Outflow") {
         hi_bc[dir] = 2;
      } else if (hi_bc_char[dir] == "Symmetry") {
         hi_bc[dir] = 3;
      } else if (hi_bc_char[dir] == "SlipWallAdiab") {
         hi_bc[dir] = 4;
      } else if (hi_bc_char[dir] == "NoSlipWallAdiab") {
         hi_bc[dir] = 5;
      } else if (hi_bc_char[dir] == "SlipWallIsotherm") {
         hi_bc[dir] = 6;
      } else if (hi_bc_char[dir] == "NoSlipWallIsotherm") {
         hi_bc[dir] = 7;
      } else {
         amrex::Abort("Wrong boundary condition word in hi_bc, please use: Interior, Inflow, Outflow, "
                      "Symmetry, SlipWallAdiab, NoSlipWallAdiab, SlipWallIsotherm, NoSlipWallIsotherm");
      }
   }

   //for (int i = 0; i < AMREX_SPACEDIM; i++) {
   //   phys_bc.setLo(i,lo_bc[i]);
   //   phys_bc.setHi(i,hi_bc[i]);
   //}
}

void PeleLM::VariablesSetup() {
   BL_PROFILE_VAR("PeleLM::VariablesSetup()", ReadParameters);

   if (m_verbose<1) return;

   std::string PrettyLine = std::string(78, '=') + "\n";

   // Variables ordering is defined through macro in PeleLM_Index.H 
   // Simply print on screen the state layout
   Print() << PrettyLine;
   Print() << " State components \n"; 
   Print() << PrettyLine;

   Print() << " Velocity X: " << VELX;
#if AMREX_SPACEDIM > 1 
   Print() << ", Velocity Y: " << VELY;
#if AMREX_SPACEDIM > 2 
   Print() << ", Velocity Z: " << VELZ;
#endif
#endif
   Print() << " \n";

   Print() << " Density: " << DENSITY << "\n";
   Print() << " First species: " << FIRSTSPEC << "\n";
   Print() << " Enthalpy: " << RHOH << "\n";
   Print() << " Temperature: " << TEMP << "\n";
   if (m_nAux > 0) {
      Print() << " First passive scalar: " << FIRSTAUX << "\n";
   }

   Print() << " => Total number of state variables: " << NVAR << "\n";

}

void PeleLM::ResizeArray() {

   if (m_verbose) {
      Print() << " Initializing data on " << max_level+1 << " levels \n";
   }

   // State data
   m_leveldata_old.resize(max_level+1);
   m_leveldata_new.resize(max_level+1);

   // Factory
   m_factory.resize(max_level+1);

}

PeleLM::LevelData::LevelData(amrex::BoxArray const& ba,
                             amrex::DistributionMapping const& dm,
                             amrex::FabFactory<FArrayBox> const& factory,
                             int nAux, int nGrowState, int nGrowMAC) :
   velocity(ba, dm, AMREX_SPACEDIM, nGrowState, MFInfo(), factory),
   density (ba, dm, 1             , nGrowState, MFInfo(), factory),
   species (ba, dm, NUM_SPECIES   , nGrowState, MFInfo(), factory),
   rhoH    (ba, dm, 1             , nGrowState, MFInfo(), factory),
   temp    (ba, dm, 1             , nGrowState, MFInfo(), factory),
   gp      (ba, dm, AMREX_SPACEDIM, 0         , MFInfo(), factory),
   press   (amrex::convert(ba,IntVect::TheNodeVector()), 
                dm, 1             , 1         , MFInfo(), factory),   
   diff_cc (ba, dm, NUM_SPECIES+2 , 1    , MFInfo(), factory)
{
   if ( nAux > 0 ) {
      auxiliaries.define(ba, dm, nAux, nGrowState, MFInfo(), factory);
   }

   velocity_mac = new MultiFab[AMREX_SPACEDIM];
   for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      const BoxArray& faceba = amrex::convert(ba,IntVect::TheDimensionVector(dir));
      velocity_mac[dir].define(faceba,dm, 1, nGrowMAC, MFInfo(), factory);
   }

}
