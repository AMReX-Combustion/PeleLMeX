#ifndef _PeleLM_H_
#define _PeleLM_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_MacProjector.H>
#include <AMReX_NodalProjector.H>
#include <AMReX_VisMF.H>

#include <mechanism.h>
#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <DiffusionOp.H>

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
   void readParameters();

   // ReadIOParameters
   void readIOParameters();

   // VariablesSetup
   void variablesSetup();

   // ResizeArray
   void resizeArray();
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // INIT

   // initData
   void initData();

   // initLevelData
   void initLevelData(int lev);

   // initialIterations
   void initialIterations();
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // PROJECTION

   void initialProjection();

   void doNodalProject(amrex::Vector<amrex::MultiFab*> &a_vel,
                       const amrex::Vector<amrex::MultiFab*> &a_sigma,
                       const amrex::Vector<amrex::MultiFab*> &rhs_cc,
                       const amrex::Vector<const amrex::MultiFab*> &rhs_nd,
                       int increment_gp,
                       amrex::Real scaling_factor);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // DIFFUSION

   // compute cell-centered diffusivity (stored in LevelData)
   void calcDiffusivity(PeleLM::TimeStamp a_time);
   
   // get edge-centered diffusivity on a per level / per comp basis
   amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>
   getDiffusivity(int lev, int beta_comp, int ncomp,
                  amrex::Vector<amrex::BCRec> a_bcrec,
                  const amrex::MultiFab &a_diff_cc);

   // get species/enthalpy differential diffusion fluxes
   void computeDifferentialDiffusionFluxes(PeleLM::TimeStamp a_time,
                                           const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_fluxes);

   // get species/enthalpy differential diffusion terms
   void computeDifferentialDiffusionTerms(PeleLM::TimeStamp a_time,
                                          amrex::Vector<amrex::MultiFab> &a_viscTerm);

   // add wbar term to the species fluxes
   void addWbarTerm(const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_fluxes,
                    amrex::Vector<amrex::MultiFab const*> const &a_spec,
                    amrex::Vector<amrex::MultiFab const*> const &a_rho,
                    amrex::Vector<amrex::MultiFab const*> const &a_beta);

   void adjustSpeciesFluxes(const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_fluxes,
                            amrex::Vector<amrex::MultiFab const*> const &a_spec);

   void computeSpeciesEnthalpyFlux(const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_fluxes,
                                   amrex::Vector<amrex::MultiFab const*> const &a_temp);
                    
   // diffusion Linear operator BC
   amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>   
   getDiffusionLinOpBC(amrex::Orientation::Side a_side,
                       const amrex::BCRec &a_bc);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // EOS

   // Thermodynamic pressure
   void setThermoPress(int lev, PeleLM::TimeStamp a_time);

   // Divergence constraint
   void calcDivU(int is_init, PeleLM::TimeStamp a_time);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // I/O
   void WritePlotFile();
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // BOUNDARY CONDITIONS
   // Convert PhysBC into field BCs
   void setBoundaryConditions();

   // FillPatch the entire state
   void fillPatchState(int lev, PeleLM::TimeStamp a_time);

   // FillPatch state components
   void fillpatch_velocity(int lev, amrex::Real a_time, amrex::MultiFab& a_vel, int nGhost);
   void fillpatch_mass(int lev, amrex::Real a_time, amrex::MultiFab& a_rho, amrex::MultiFab& a_rhoY, int nGhost);
   void fillpatch_energy(int lev, amrex::Real a_time, amrex::MultiFab& a_rhoh, amrex::MultiFab& a_temp, int nGhost);

   // Fill physical boundaries
   void setPhysBoundaryVel (amrex::MultiFab &a_vel, int lev, PeleLM::TimeStamp a_time);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // UTILS

   // flux divergence
   void fluxDivergence(amrex::Vector<amrex::MultiFab> &a_divergence,
                       int div_comp,
                       const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> > &a_fluxes, 
                       int flux_comp,
                       int ncomp,
                       amrex::Real scale);

   void fluxDivergenceLevel(int lev,
                            amrex::MultiFab &a_divergence,
                            int div_comp,
                            const amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                            int flux_comp,
                            int ncomp,
                            amrex::Real scale);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Gather data of a level into a struct
   struct LevelData {
      LevelData () = default;
      LevelData (amrex::BoxArray const& ba,
                 amrex::DistributionMapping const& dm,
                 amrex::FabFactory<amrex::FArrayBox> const& factory,
                 int a_incompressible, int a_has_divu, 
                 int a_nAux, int a_nGrowState, int a_nGrowMAC);

      // cell-centered state multifabs
      amrex::MultiFab velocity;        // Velocity (dim:AMREX_SPACEDIM)   
      amrex::MultiFab density;         // Density (dim:1)
      amrex::MultiFab species;         // Species mass fraction (dim:NUM_SPECIES)
      amrex::MultiFab rhoh;            // Rho * enthalpy (dim:1)
      amrex::MultiFab rhoRT;           // Thermodynamic pressure (dim:1)
      amrex::MultiFab temp;            // Temperature (dim:1) 
      amrex::MultiFab auxiliaries;     // Auxiliary variables (passive scalars and others) (dim:m_nAux) 
      amrex::MultiFab gp;              // pressure gradient (dim:AMREX_SPACEDIM)

      // node-centered state multifabs
      amrex::MultiFab press;           // nodal pressure (dim:1)

      // cell-centered advance multifabs
      amrex::MultiFab divu;            // Velocity divergence constraint
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

   amrex::Real getTime(int lev, PeleLM::TimeStamp time) const {
      AMREX_ASSERT(time==AmrOldTime || time==AmrNewTime);
      if ( time == AmrOldTime ) {
         return m_t_old[lev];
      } else {
         return m_t_new[lev];
      }
   }

   amrex::Vector<amrex::MultiFab* > getSpeciesVect(PeleLM::TimeStamp time);
   amrex::Vector<amrex::MultiFab* > getDensityVect(PeleLM::TimeStamp time);
   amrex::Vector<amrex::MultiFab* > getTempVect(PeleLM::TimeStamp time);
   amrex::Vector<amrex::MultiFab* > getDiffusivityVect(PeleLM::TimeStamp time);

   amrex::Vector<amrex::BCRec> fetchBCRecArray(int scomp, int ncomp);
   DiffusionOp* getDiffusionOp ();
   //-----------------------------------------------------------------------------

   //protected:

   // Class data
   amrex::Vector<std::unique_ptr<LevelData> > m_leveldata_old;
   amrex::Vector<std::unique_ptr<LevelData> > m_leveldata_new;

   // FAB factory
   amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > > m_factory;

   // Times
   amrex::Vector<amrex::Real> m_t_old;
   amrex::Vector<amrex::Real> m_t_new;

   // Boundary Conditions
   amrex::BCRec  m_phys_bc;                      // Store physical BC from input
   amrex::Vector<amrex::BCRec> m_bcrec_state;

   // Problem parameters
   std::unique_ptr<ProbParm> prob_parm;

   // Linear solvers
   std::unique_ptr<DiffusionOp> m_diffusion_op;

   int m_verbose = 0;

   // IO options
   std::string m_restart_file{""};
   std::string m_plot_file{"plt"};
   int m_plot_int = 0;

   amrex::Real m_cur_time     = 0.0;
   amrex::Real m_dt           = -1.0;

   //-----------------------------------------------------------------------------
   // ALGORITHM

   // Switch Evolve/Evaluate
   std::string m_run_mode;

   // Incompressible ?
   int m_incompressible = 0;
   amrex::Real m_rho = 1.175;

   // DivU ? default ON
   int m_has_divu = 1;

   // Number if passive scalars
   int m_nAux = 0;

   // R-z
   int m_rz_correction = 0;

   // cc->ec average
   int m_harm_avg_cen2edge = 0;

   // wbar diffusion term
   int m_use_wbar = 1;

   // Number of ghost cells
   int m_nGrowState = 1;
   int m_nGrowMAC = 0;

   int m_nstep = -1; 

   // Init proj flag
   int m_do_init_proj = 1;
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Linear Solvers

   // Nodal projection
   int m_nodal_mg_max_coarsening_level = 100;
   amrex::Real m_nodal_mg_rtol = 1.0e-11;   
   amrex::Real m_nodal_mg_atol = 1.0e-14;
   //-----------------------------------------------------------------------------
};

void writeBuildInfo();

#endif
