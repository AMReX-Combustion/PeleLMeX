#ifndef _PeleLM_H_
#define _PeleLM_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <PeleLM_Index.H>
#include <mechanism.h>
#include <pelelm_prob_parm.H>
#include <AMReX_VisMF.H>

class PeleLM : public amrex::AmrCore {

   public:

   enum TimeStamp { AmrOldTime, AmrHalfTime, AmrNewTime };

   // constructor
   PeleLM();

   //destructor
   virtual ~PeleLM();

   // Setup function
   void Setup();

   // Init function
   void Init();

   // Advance function
   void Evolve();

   // Unit Test function
   void Evaluate();

   //-----------------------------------------------------------------------------
   // Virtual AmrCore functions
   virtual void MakeNewLevelFromScratch( int lev, 
                                               amrex::Real time,
                                         const amrex::BoxArray& ba,
                                         const amrex::DistributionMapping& dm) override;

   virtual void ErrorEst( int lev, 
                         amrex::TagBoxArray& tags,
                         amrex::Real time,
                         int ng) override;

   virtual void MakeNewLevelFromCoarse( int lev, 
                                              amrex::Real time,
                                        const amrex::BoxArray& ba,
                                        const amrex::DistributionMapping& dm) override;

   virtual void RemakeLevel( int lev,
                                   amrex::Real time,
                             const amrex::BoxArray& ba,
                             const amrex::DistributionMapping& dm) override;

   virtual void ClearLevel(int lev) override;
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // SETUP

   // ReadParameters
   void ReadParameters();

   // VariablesSetup
   void VariablesSetup();

   // ResizeArray
   void ResizeArray();
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // INIT

   // initData
   void initData();

   // initLevelData
   void initLevelData(int lev);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // DIFFUSION

   // compute cell-centered diffusivity
   void calcDiffusivity(PeleLM::TimeStamp time);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // I/O
   void WritePlotFile();
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Gather data of a level into a struct
   struct LevelData {
      LevelData () = default;
      LevelData (amrex::BoxArray const& ba,
                 amrex::DistributionMapping const& dm,
                 amrex::FabFactory<amrex::FArrayBox> const& factory,
                 int nAux, int nGrowState, int nGrowMAC);

      // cell-centered state multifabs
      amrex::MultiFab velocity;        // Velocity (dim:AMREX_SPACEDIM)   
      amrex::MultiFab density;         // Density (dim:1)
      amrex::MultiFab species;         // Species mass fraction (dim:NUM_SPECIES)
      amrex::MultiFab rhoH;            // Rho * enthalpy (dim:1)
      amrex::MultiFab temp;            // Temperature (dim:1) 
      amrex::MultiFab auxiliaries;     // Auxiliary variables (passive scalars and others) (dim:m_nAux) 
      amrex::MultiFab gp;              // pressure gradient (dim:AMREX_SPACEDIM)

      // node-centered state multifabs
      amrex::MultiFab press;           // nodal pressure (dim:1)

      // cell-centered advance multifabs
      amrex::MultiFab diff_cc;         // Diffusivity (dim:NUM_SPECIES+2)

      // face-centered multifabs
      amrex::MultiFab* velocity_mac;
   };
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Accessors
   std::string runMode() {
      return m_run_mode;
   }

   amrex::FabFactory<amrex::FArrayBox> const&
   Factory (int lev) const noexcept { return *m_factory[lev]; }

   LevelData* getLevelDataPtr(int lev, PeleLM::TimeStamp time) const {
      AMREX_ASSERT(time==AmrOldTime || time==AmrNewTime);
      if ( time == AmrOldTime ) { 
         return m_leveldata_old[lev].get();
      } else {
         return m_leveldata_new[lev].get();
      }
   }
   //-----------------------------------------------------------------------------

   protected:

   // Class data
   amrex::Vector<std::unique_ptr<LevelData> > m_leveldata_old;
   amrex::Vector<std::unique_ptr<LevelData> > m_leveldata_new;

   // FAB factory
   amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > > m_factory;

   // Problem parameters
   std::unique_ptr<ProbParm> prob_parm;

   int m_verbose = 0;

   std::string m_restart_file{""};
   std::string m_plot_file{"plt"};

   amrex::Real m_cur_time     = 0.0;

   //-----------------------------------------------------------------------------
   // ALGORITHM

   // Switch Evolve/Evaluate
   std::string m_run_mode;

   // Number if passive scalars
   int m_nAux = 0;

   // Number of ghost cells
   int m_nGrowState = 0;
   int m_nGrowMAC = 0;

   int m_nstep = -1; 

   // Init proj flag
   int m_do_init_proj = 1;
   //-----------------------------------------------------------------------------

};

void writeBuildInfo();

#endif
