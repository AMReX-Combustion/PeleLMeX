#ifndef _PeleLM_H_
#define _PeleLM_H_

#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BC_TYPES.H>
#include <hydro_MacProjector.H>
#include <hydro_NodalProjector.H>
#include <AMReX_ErrorList.H>
#include <AMReX_VisMF.H>

#include <mechanism.H>
#include "PelePhysics.H"
#include "PMFData.H"
#include "turbinflow.H"
#include "ReactorBase.H"
#include <PeleLM_Index.H>
#include <PeleLMDerive.H>
#include <pelelm_prob_parm.H>
#include <DiffusionOp.H>
#include <DiagBase.H>

#ifdef PELE_USE_EFIELD
#include <PrecondOp.H>
#endif

// Forward declarations
#ifdef PELELM_USE_SPRAY
class SprayParticleContainer;
#endif

#ifdef PELELM_USE_SOOT
class SootModel;
#endif

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

   // Simulation function
   void Evolve();

   // Unit Test function
   void Evaluate();

   // Time advance function
   void Advance(int is_init);

   //-----------------------------------------------------------------------------
   // Virtual AmrCore functions
   virtual void regrid(int lbase,
                       amrex::Real time,
                       bool initial = false) override;

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

   // ReadProblemParameters
   void readProbParm();

   // ReadGridFile
  void readGridFile(std::string grid_file,
                    amrex::Vector<amrex::BoxArray>& input_ba);

#ifdef PELELM_USE_SPRAY
   // ReadSprayParameters
   void readSprayParameters();

   // SprayParticleSetup
   void sprayParticleSetup();

   // CreateSprayData
   void createSprayData();

   // InitSprays
   void initSprays();

   // SprayRestart
   void sprayRestart();
#endif

   // VariablesSetup
   void variablesSetup();

   // DerivedSetup
   void derivedSetup();

   // EvaluateSetup
   void evaluateSetup();

   // TaggingSetup
   void taggingSetup();

   // ResizeArray
   void resizeArray();
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Gather State Data of a level into a struct
   struct LevelData {
      LevelData () = default;
      LevelData (amrex::BoxArray const& ba,
                 amrex::DistributionMapping const& dm,
                 amrex::FabFactory<amrex::FArrayBox> const& factory,
                 int a_incompressible, int a_has_divu, 
                 int a_nAux, int a_nGrowState);

      // cell-centered state multifabs
      amrex::MultiFab state;           // Stqte data: velocity, density, rhoYs, rhoH, Temp, RhoRT
      amrex::MultiFab auxiliaries;     // Auxiliary variables (passive scalars and others) (dim:m_nAux)
      amrex::MultiFab gp;              // pressure gradient (dim:AMREX_SPACEDIM)
      amrex::MultiFab divu;            // Velocity divergence constraint

      // node-centered state multifabs
      amrex::MultiFab press;           // nodal pressure (dim:1)

      // cell-centered transport multifabs
      amrex::MultiFab visc_cc;         // Viscosity (dim:1)
      amrex::MultiFab diff_cc;         // Diffusivity (dim:NUM_SPECIES+2)
#ifdef PELE_USE_EFIELD
      amrex::MultiFab diffE_cc;        // Electron diffusivity (dim:1)
      amrex::MultiFab mobE_cc;         // Electron mobility (dim:1)
      amrex::MultiFab mob_cc;          // Species mobility (dim:NUM_IONS)
#endif
   };

   struct LevelDataReact {
      LevelDataReact () = default;
      LevelDataReact (const amrex::BoxArray &ba,
                      const amrex::DistributionMapping &dm,
                      const amrex::FabFactory<amrex::FArrayBox> &factory);
      amrex::MultiFab  I_R;            // Species reaction rates
      amrex::MultiFab  functC;         // Implicit integrator function call count
#ifdef PELE_USE_EFIELD
      amrex::MultiFab  I_RnE;          // Electron number density reaction term
#endif
   };
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Gather Advance Datas into a structs
   struct AdvanceDiffData {
      AdvanceDiffData () = default;
      AdvanceDiffData (int a_finestLevel,
                       const amrex::Vector<amrex::BoxArray> &ba,
                       const amrex::Vector<amrex::DistributionMapping> &dm,
                       const amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>>> &factory,
                       int a_nGrowAdv,
                       int a_use_wbar,
                       int is_init = 0);
      amrex::Vector<amrex::MultiFab> Dn;           // Diffusion term t^n
      amrex::Vector<amrex::MultiFab> Dnp1;         // Diffusion term t^(n+1,k)
      amrex::Vector<amrex::MultiFab> Dhat;         // Diffusion term t^(n+1,k+1)

      amrex::Vector<amrex::MultiFab> Dwbar;        // Diffusion term of Wbar
      amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> wbar_fluxes;     // Wbar flux correction
   };
   struct AdvanceAdvData {
      AdvanceAdvData () = default;
      AdvanceAdvData (int a_finestLevel,
                      const amrex::Vector<amrex::BoxArray> &ba,
                      const amrex::Vector<amrex::DistributionMapping> &dm,
                      const amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>>> &factory,
                      int a_incompressible,
                      int a_nGrowAdv,
                      int a_nGrowMAC);
      amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> umac;        // MAC face velocity
      amrex::Vector<amrex::MultiFab>  AofS;        // Scalar advection term at t^(n+1/2)
      amrex::Vector<amrex::MultiFab>  chi;         // Thermodynamic contraint
      amrex::Vector<amrex::MultiFab>  Forcing;     // Scalar forcing for both advection and diffusion
      amrex::Vector<amrex::MultiFab>  mac_divu;    // divu used in MAC projection
#ifdef PELE_USE_EFIELD
      amrex::Vector<amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>> uDrift;        // ions drift face velocity
#endif

   };
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // INIT

   // initData
   void initData();

   // initLevelData
   void initLevelData(int lev);

   // initialIterations
   void initialIterations();

   // initFromGridFile
   void InitFromGridFile(amrex::Real time);

   // regridFromFile
   void regridFromGridFile(int lbase,
                           amrex::Real time,
                           bool initial);

   // checkRunParams
   void checkRunParams();

   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // PROJECTION

   // Initial projections
   void initialProjection();

   void initialPressProjection();

   // Velocity projection at Time
   void velocityProjection(int is_init,
                           const PeleLM::TimeStamp &a_rhoTime,
                           const amrex::Real &a_dt);

   // Actual projection function
   void doNodalProject(const amrex::Vector<amrex::MultiFab*> &a_vel,
                       const amrex::Vector<amrex::MultiFab*> &a_sigma,
                       const amrex::Vector<amrex::MultiFab*> &rhs_cc,
                       const amrex::Vector<const amrex::MultiFab*> &rhs_nd,
                       int increment_gp,
                       amrex::Real scaling_factor);

   // Scaling functions for 2D-RZ calculations
   void scaleProj_RZ(int a_lev,
                     amrex::MultiFab &a_rhs);
   void unscaleProj_RZ(int a_lev,
                       amrex::MultiFab &a_rhs);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // MAC PROJECTION

   void resetMacProjector();
   
   // Predict face velocity using Godunov
   void predictVelocity(std::unique_ptr<AdvanceAdvData> &advData);

   // Create the MAC projection RHS
   void createMACRHS(std::unique_ptr<AdvanceAdvData> &advData);

   // Append the chi increment
   void addChiIncrement(int a_sdcIter,
                        const PeleLM::TimeStamp &a_time,
                        std::unique_ptr<AdvanceAdvData> &advData);

   // Project face velocity
   void macProject(const PeleLM::TimeStamp &a_time,
                   std::unique_ptr<AdvanceAdvData> &advData,
                   const amrex::Vector<amrex::MultiFab*> &a_divu);

   // Umac ghost cells
   void create_constrained_umac_grown(int a_lev, int a_nGrow,
                                      const amrex::Geometry* crse_geom,
                                      const amrex::Geometry* fine_geom,
                                      amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> u_mac_crse,
                                      amrex::Array<amrex::MultiFab *, AMREX_SPACEDIM> u_mac_fine,
                                      const amrex::MultiFab *divu,
                                      const amrex::IntVect& crse_ratio);

   // MAC proj BC
   amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
   getMACProjectionBC(amrex::Orientation::Side a_side);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // DIFFUSION

   // compute cell-centered diffusivity (stored in LevelData)
   void calcViscosity(const PeleLM::TimeStamp &a_time);
   void calcDiffusivity(const PeleLM::TimeStamp &a_time);
   
   // get edge-centered diffusivity on a per level / per comp basis
   amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>
   getDiffusivity(int lev, int beta_comp, int ncomp, int doZeroVisc,
                  amrex::Vector<amrex::BCRec> a_bcrec,
                  const amrex::MultiFab &a_diff_cc);

   // get species/enthalpy differential diffusion fluxes
   void computeDifferentialDiffusionFluxes(const PeleLM::TimeStamp &a_time,
                                           const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_fluxes,
                                           const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_wbarfluxes);

   // get species/enthalpy differential diffusion terms
   void computeDifferentialDiffusionTerms(const PeleLM::TimeStamp &a_time,
                                          std::unique_ptr<AdvanceDiffData> &diffData,
                                          int is_init = 0);

   // add wbar term to the species fluxes
   void addWbarTerm(const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_fluxes,
                    const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_wbarfluxes, 
                    amrex::Vector<amrex::MultiFab const*> const &a_spec,
                    amrex::Vector<amrex::MultiFab const*> const &a_rho,
                    amrex::Vector<amrex::MultiFab const*> const &a_beta);

   // Ensure sum of specues diffusion fluxes sum to 0
   void adjustSpeciesFluxes(const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_fluxes,
                            amrex::Vector<amrex::MultiFab const*> const &a_spec);

   // Compute enthalpy flux due to species diffusion
   void computeSpeciesEnthalpyFlux(const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_fluxes,
                                   amrex::Vector<amrex::MultiFab const*> const &a_temp);
                    
   void differentialDiffusionUpdate(std::unique_ptr<AdvanceAdvData> &advData,
                                    std::unique_ptr<AdvanceDiffData> &diffData);

   void deltaTIter_prepare(const amrex::Vector<amrex::MultiFab*> &a_rhs,
                           const amrex::Vector<amrex::MultiFab*> &a_Tsave,
                           const amrex::Vector<amrex::MultiFab*> &a_rhoCp,
                           std::unique_ptr<AdvanceAdvData> &advData,
                           std::unique_ptr<AdvanceDiffData> &diffData);

   void deltaTIter_update(int a_dtiter,
                          const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_fluxes,
                          const amrex::Vector<amrex::MultiFab const*> &a_Tsave,
                          std::unique_ptr<AdvanceDiffData> &diffData,
                          amrex::Real &a_deltaT_norm);


   // Viscous forces
   void computeDivTau(const PeleLM::TimeStamp &a_time,
                      const amrex::Vector<amrex::MultiFab*> &a_divtau,
                      int use_density, amrex::Real scale = 1.0);

   // Diffuse velocity
   void diffuseVelocity();

   // diffusion Linear operator BC
   amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>
   getDiffusionLinOpBC(amrex::Orientation::Side a_side,
                       const amrex::BCRec &a_bc);

   // diffusion Tensor operator BC
   amrex::Vector<amrex::Array<amrex::LinOpBCType,AMREX_SPACEDIM>>
   getDiffusionTensorOpBC(amrex::Orientation::Side a_side,
                          const amrex::Vector<amrex::BCRec> a_bc);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // REACTION

   void advanceChemistry(std::unique_ptr<AdvanceAdvData> &advData);

   void advanceChemistry(int lev,
                         const amrex::Real &a_dt,
                         amrex::MultiFab &a_extForcing);

   void advanceChemistry(int lev,
                         const amrex::Real &a_dt,
                         amrex::MultiFab &a_extForcing,
                         amrex::MultiFab *avgDownIR);

   void computeInstantaneousReactionRate(const amrex::Vector<amrex::MultiFab*> &a_I_R,
                                         const PeleLM::TimeStamp &a_time);

   void computeInstantaneousReactionRate(int lev,
                                         const PeleLM::TimeStamp &a_time,
                                         amrex::MultiFab* a_I_R);
   void getHeatRelease(int lev, amrex::MultiFab *a_HR);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Time step size

   amrex::Real computeDt(int is_init, const PeleLM::TimeStamp &a_time);
   amrex::Real estConvectiveDt(const PeleLM::TimeStamp &a_time);
   amrex::Real estDivUDt(const PeleLM::TimeStamp &a_time);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // FORCES

   void getVelForces(const TimeStamp &a_time,
                     const amrex::Vector<amrex::MultiFab*> &a_divTau,
                     const amrex::Vector<amrex::MultiFab*> &a_velForce,
                     int nGrowForce,
                     int add_gradP);
   
   void getVelForces(const TimeStamp &a_time,
                     int lev,
                     amrex::MultiFab* a_divTau,
                     amrex::MultiFab* a_velForce,
                     int add_gradP);

   void getVelForces(int lev,
                     const amrex::Box &bx,
                     const amrex::Real &a_time,
                     amrex::Array4<amrex::Real> const& force,
                     amrex::Array4<const amrex::Real> const& vel,
                     amrex::Array4<const amrex::Real> const& rho,
                     amrex::Array4<const amrex::Real> const& rhoY,
                     amrex::Array4<const amrex::Real> const& rhoh,
                     amrex::Array4<const amrex::Real> const& temp,
                     amrex::Array4<const amrex::Real> const& extMom,
                     amrex::Array4<const amrex::Real> const& extRho);
   //-----------------------------------------------------------------------------

#ifdef PELELM_USE_SPRAY
   //-----------------------------------------------------------------------------
   // Sprays
   void setupVirtualParticles(const int level);
   void removeVirtualParticles(const int level);
   void setupGhostParticles(const int ngrow, const int level);
   void removeGhostParticles(const int level);
   amrex::Real estSprayDt();
   void sprayMKD(const amrex::Real time, const amrex::Real dt);
   void sprayMKDLevel(const int level,
                      const amrex::Real time,
                      const amrex::Real dt);
   void sprayInjectRedist();
   void sprayPostRegrid();
   void setSprayState(const amrex::Real& a_flow_dt);
   void addSpraySource(const int level);
   static int write_spray_ascii_files;
   static int num_spray_src;
   static bool do_spray_particles;
   static int plot_spray_src;
   static std::string spray_fuel_names[SPRAY_FUEL_NUM];
   static amrex::Vector<std::string> spray_derive_vars;
   static SprayParticleContainer* theSprayPC();
   static SprayParticleContainer* theVirtPC();
   static SprayParticleContainer* theGhostPC();
   // State used for spray interpolation
   amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_spraysource;
   amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_spraystate;
   //-----------------------------------------------------------------------------
#endif
#ifdef PELELM_USE_SOOT
   //-----------------------------------------------------------------------------
   // Soot

   // SetSootIndx
   void setSootIndx();

   // CleanupSootModel
   void cleanupSootModel();

   // AddSootDerivePlotVars
   //void addSootDerivePlotVars(PeleLMDeriveList& derive_lst);

   // Compute the soot source terms
   void computeSootSource(const PeleLM::TimeStamp &a_timestamp,
                          const amrex::Real a_dt);

   // Clip soot moments
   void clipSootMoments();

   //-----------------------------------------------------------------------------
#endif
   //-----------------------------------------------------------------------------
   // EOS

   // Thermodynamic pressure
   void setThermoPress(const PeleLM::TimeStamp &a_time);
   void setThermoPress(int lev, const PeleLM::TimeStamp &a_time);

   // update Temperature from rhoY and rhoH
   void setTemperature(const PeleLM::TimeStamp &a_time);
   void setTemperature(int lev, const PeleLM::TimeStamp &a_time);

   // Divergence constraint
   void calcDivU(int is_init,
                 int computeDiff,
                 int do_avgDown,
                 const PeleLM::TimeStamp &a_time,
                 std::unique_ptr<AdvanceDiffData> &diffData);

   // Pressure drift
   void calc_dPdt(const PeleLM::TimeStamp &a_time,
                  const amrex::Vector<amrex::MultiFab*> &a_dPdt);
   void calc_dPdt(int lev, const PeleLM::TimeStamp &a_time,
                  amrex::MultiFab* a_dPdt);

   // Closed Chamber
   amrex::Real
   adjustPandDivU(std::unique_ptr<AdvanceAdvData> &advData);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // I/O
   void WritePlotFile();
   bool writePlotNow();
   bool checkMessage(const std::string &a_action);
   void WriteCheckPointFile();
   void ReadCheckPointFile();
   bool writeCheckNow();
   void WriteJobInfo(const std::string &path) const;
   void WriteHeader(const std::string &name, bool is_checkpoint) const;
   void WriteDebugPlotFile(const amrex::Vector<const amrex::MultiFab*> &a_MF,
                           const std::string &pltname);
   void initLevelDataFromPlt(int a_lev, const std::string &a_pltFile);
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Diagnostics
   void createDiagnostics();
   void updateDiagnostics();
   void doDiagnostics();

   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Temporal
   void massBalance();
   void speciesBalance();
   void rhoHBalance();
   void initTemporals(const PeleLM::TimeStamp &a_time = AmrOldTime);
   void writeTemporals();
   void addMassFluxes(const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> &a_fluxes, const amrex::Geometry &a_geom);
   void addRhoHFluxes(const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> &a_fluxes, const amrex::Geometry &a_geom);
   void addRhoYFluxes(const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                      const amrex::Geometry &a_geom,
                      const amrex::Real &a_factor = 1.0);
   void addUmacFluxes(std::unique_ptr<AdvanceAdvData> &advData, const amrex::Geometry &a_geom);
   void openTempFile();
   void closeTempFile();
   bool doTemporalsNow();
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // BOUNDARY CONDITIONS / FILLPATCH

   amrex::InterpBase* getInterpolator();
   // Convert PhysBC into field BCs
   void setBoundaryConditions();

   // FillPatch the entire class state
   void fillPatchState(const PeleLM::TimeStamp &a_time);
   void fillPatchState(int lev, const PeleLM::TimeStamp &a_time);

   // FillPatch the various state components
   void fillPatchDensity(const PeleLM::TimeStamp &a_time);
   void fillPatchSpecies(const PeleLM::TimeStamp &a_time);
   void fillPatchTemp(const PeleLM::TimeStamp &a_time);
#ifdef PELE_USE_EFIELD
   void fillPatchPhiV(const PeleLM::TimeStamp &a_time);
#endif

   // FillPatch and return a state MF
   std::unique_ptr<amrex::MultiFab>
   fillPatchState(int lev, amrex::Real a_time, int nGrow);

   // FillPatch state components
   void fillpatch_state(int lev, amrex::Real a_time, amrex::MultiFab& a_state, int nGhost);
   void fillpatch_density(int lev, amrex::Real a_time, amrex::MultiFab& a_rhoY, int rho_comp, int nGhost);
   void fillpatch_species(int lev, amrex::Real a_time, amrex::MultiFab& a_rhoY, int rhoY_comp, int nGhost);
   void fillpatch_temp(int lev, amrex::Real a_time,
                       amrex::MultiFab& a_temp, int temp_comp, int nGhost);
   void fillpatch_divu(int lev, amrex::Real a_time, amrex::MultiFab& a_divu, int nGhost);
   void fillpatch_gradp(int lev, amrex::Real a_time, amrex::MultiFab& a_gp, int nGhost);
   void fillpatch_reaction(int lev, amrex::Real a_time, amrex::MultiFab& a_I_R, int nGhost);
   void fillpatch_forces(amrex::Real a_time, amrex::Vector<amrex::MultiFab*> const &a_velForces, int nGrowForce);
#ifdef PELE_USE_EFIELD
   void fillpatch_phiV(int lev, amrex::Real a_time, amrex::MultiFab& a_phiV, int phiV_comp, int nGhost);
   void fillPatchExtrap(amrex::Real a_time, amrex::Vector<amrex::MultiFab*> const &a_MF, int a_nGrow);
   void fillPatchNLnE(amrex::Real a_time, amrex::Vector<amrex::MultiFab*> const &a_nE, int a_nGrow);
   void fillPatchNLphiV(amrex::Real a_time, amrex::Vector<amrex::MultiFab*> const &a_phiV, int a_nGrow);
#endif

   // FillCoarsePatch state components
   void fillcoarsepatch_state(int lev, amrex::Real a_time, amrex::MultiFab& a_state, int nGhost);
   void fillcoarsepatch_divu(int lev, amrex::Real a_time, amrex::MultiFab& a_divu, int nGhost);
   void fillcoarsepatch_gradp(int lev, amrex::Real a_time, amrex::MultiFab& a_gp, int nGhost);
   void fillcoarsepatch_reaction(int lev, amrex::Real a_time, amrex::MultiFab& a_I_R, int nGhost);

   // Fill physical boundaries
   void setInflowBoundaryVel (amrex::MultiFab &a_vel, int lev, PeleLM::TimeStamp a_time);

   void fillTurbInflow(amrex::MultiFab &a_vel, int vel_comp, int lev, const amrex::Real a_time);

   // Average down operations
   void averageDownState(const PeleLM::TimeStamp &a_time);
   void averageDownScalars(const PeleLM::TimeStamp &a_time);
   void averageDownVelocity(const PeleLM::TimeStamp &a_time);
   void averageDownDensity(const PeleLM::TimeStamp &a_time);
   void averageDownRhoRT(const PeleLM::TimeStamp &a_time);
   void averageDown(const PeleLM::TimeStamp &a_time, const int state_comp, const int ncomp);
#ifdef PELE_USE_EFIELD
   void averageDownnE(const PeleLM::TimeStamp &a_time);
   void averageDownPhiV(const PeleLM::TimeStamp &a_time);
#endif
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // UTILS

   // flux divergence
   void fluxDivergence(const amrex::Vector<amrex::MultiFab*> &a_divergence,
                       int div_comp,
                       const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> > &a_fluxes,
                       int flux_comp,
                       int ncomp,
                       int intensiveFluxes,
                       amrex::Real scale);

   void fluxDivergenceRD(const amrex::Vector<const amrex::MultiFab*> &a_state,
                         int state_comp,
                         const amrex::Vector<amrex::MultiFab*> &a_divergence,
                         int div_comp,
                         const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> > &a_fluxes,
                         int flux_comp,
                         int ncomp,
                         int intensiveFluxes,
                         const amrex::BCRec *state_bc_d,
                         const amrex::Real &scale,
                         const amrex::Real &a_dt);

   void advFluxDivergence(int a_lev,
                          amrex::MultiFab &a_divergence, int div_comp,
                          amrex::MultiFab &a_divu,
                          const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                          int flux_comp,
                          const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> &a_faceState,
                          int face_comp,
                          int ncomp,
                          int const* l_conserv_d,
                          const amrex::Geometry &a_geom,
                          amrex::Real scale,
                          bool fluxes_are_area_weighted);

   void extFluxDivergenceLevel(int lev,
                               amrex::MultiFab &a_divergence,
                               int div_comp,
                               const amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                               int flux_comp,
                               int ncomp,
                               amrex::Real scale);

   void intFluxDivergenceLevel(int lev,
                               amrex::MultiFab &a_divergence,
                               int div_comp,
                               const amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                               int flux_comp,
                               int ncomp,
                               amrex::Real scale);

   void floorSpecies(const PeleLM::TimeStamp &a_time);

   amrex::Real MLNorm0(const amrex::Vector<const amrex::MultiFab*> &a_MF);
   amrex::Vector<amrex::Real>
   MLNorm0(const amrex::Vector<const amrex::MultiFab*> &a_MF,
           int startcomp, int ncomp);

   amrex::Real MFSum(const amrex::Vector<const amrex::MultiFab*> &a_MF, int comp);
   amrex::Real MFmax(const amrex::MultiFab *a_MF,
                     const amrex::iMultiFab &a_mask,
                     int comp);
   amrex::Real MFmin(const amrex::MultiFab *a_MF,
                     const amrex::iMultiFab &a_mask,
                     int comp);

   amrex::Vector<amrex::Real>
   MLmax(const amrex::Vector<const amrex::MultiFab*> &a_MF,
         int scomp, int ncomp);

   amrex::Vector<amrex::Real>
   MLmin(const amrex::Vector<const amrex::MultiFab*> &a_MF,
         int scomp, int ncomp);

   void resetCoveredMask();

   // Derived variable evaluation
   std::unique_ptr<amrex::MultiFab> derive(const std::string &a_name,
                                           amrex::Real a_time, int lev, int nGrow);
   std::unique_ptr<amrex::MultiFab> deriveComp(const std::string &a_name,
                                               amrex::Real a_time, int lev, int nGrow);

   // Evaluate function
   void MLevaluate(const amrex::Vector<amrex::MultiFab*> &a_MFVec,
                   int a_comp,
                   int &nComp,
                   const std::string &a_var);

   // Lightweight state description access
   bool isStateVariable(const std::string &a_name);
   int stateVariableIndex(const std::string &a_name);

   // Typical values
   void setTypicalValues(const PeleLM::TimeStamp &a_time, int is_init = 0);
   void updateTypicalValuesChem();
    
   void checkMemory(const std::string &a_message);

   // Mixture fraction & Progress variable
   void initMixtureFraction();
   void parseComposition(amrex::Vector<std::string> compositionIn,
                         std::string                compositionType,
                         amrex::Real               *massFrac);
   void parseVars(const amrex::Vector<std::string> &a_varNames,
                  const amrex::Vector<std::string> &a_stringIn,
                  amrex::Vector<amrex::Real> &a_rVars);
   void initProgressVariable();
   //-----------------------------------------------------------------------------

#ifdef AMREX_USE_EB
   //-----------------------------------------------------------------------------
   // EB

   void makeEBGeometry();

   int getRestartEBMaxLevel();

   void initialRedistribution();

   void initCoveredState();
   void setCoveredState(const PeleLM::TimeStamp &a_time);
   void setCoveredState(int lev, const PeleLM::TimeStamp &a_time);

   void redistributeAofS(int a_lev,
                         amrex::Real &a_dt,
                         amrex::MultiFab &a_tmpDiv, int div_comp,
                         amrex::MultiFab &a_AofS, int aofs_comp,
                         amrex::MultiFab &a_state, int state_comp,
                         int ncomp,
                         const amrex::BCRec * d_bc,
                         const amrex::Geometry &a_geom);

   void redistributeDiff(int a_lev,
                         const amrex::Real &a_dt,
                         amrex::MultiFab &a_tmpDiv, int div_comp,
                         amrex::MultiFab &a_AofS, int aofs_comp,
                         const amrex::MultiFab &a_state, int state_comp,
                         int ncomp,
                         const amrex::BCRec * d_bc,
                         const amrex::Geometry &a_geom);


   void getEBDistance(int a_lev,
                      amrex::MultiFab &a_signDistLev); 

   void extendSignedDistance( amrex::MultiFab *a_signDist,
                              amrex::Real a_extendFactor);

   // static_cast factory into ebfactory
   amrex::EBFArrayBoxFactory const& 
   EBFactory (int lev) const noexcept { 
      return static_cast<amrex::EBFArrayBoxFactory const&>(*m_factory[lev]);
   }

   void correct_vel_small_cells(amrex::Vector<amrex::MultiFab*> const &a_vel,
                                amrex::Vector<amrex::Array<amrex::MultiFab const*,AMREX_SPACEDIM> > const &a_umac);
#endif


   //-----------------------------------------------------------------------------
   // ADVANCE

   void advanceSetup();

   void oneSDC(int sdcIter, 
               std::unique_ptr<AdvanceAdvData> &advData,
               std::unique_ptr<AdvanceDiffData> &diffData);
   
   void copyTransportOldToNew();
   void copyStateNewToOld(int nGhost = 0);
   void copyPressNewToOld();
   void copyStateOldToNew(int nGhost = 0);
   void copyDiffusionOldToNew(std::unique_ptr<AdvanceDiffData> &diffData);

   void computeVelocityAdvTerm(std::unique_ptr<AdvanceAdvData> &advData);
   void updateVelocity(std::unique_ptr<AdvanceAdvData> &advData);

   void getScalarAdvForce(std::unique_ptr<AdvanceAdvData> &advData,
                          std::unique_ptr<AdvanceDiffData> &diffData);
   void computeScalarAdvTerms(std::unique_ptr<AdvanceAdvData> &advData);
   void computePassiveAdvTerms(std::unique_ptr<AdvanceAdvData> &advData, int state_comp, int ncomp);
   void updateDensity(std::unique_ptr<AdvanceAdvData> &advData);
   void updateScalarComp(std::unique_ptr<AdvanceAdvData> &advData, int state_comp, int ncomp);

   void getScalarDiffForce(std::unique_ptr<AdvanceAdvData> &advData,
                           std::unique_ptr<AdvanceDiffData> &diffData);

   void getScalarReactForce(std::unique_ptr<AdvanceAdvData> &advData);
   //-----------------------------------------------------------------------------

#ifdef PELE_USE_EFIELD
#include <PeleLMEF.H>
   //-----------------------------------------------------------------------------
   // EFIELD

   void poissonSolveEF(const PeleLM::TimeStamp &a_time);

   void calcEFTransport(const PeleLM::TimeStamp &a_time);

   void initializeElectronNeutral();
   void initializeElectronFromMassFraction();

   void implicitNonLinearSolve(int sdcIter,
                               const amrex::Real &a_dt,
                               std::unique_ptr<AdvanceDiffData> &diffData,
                               std::unique_ptr<AdvanceAdvData> &advData);
   void computeBGcharge(const amrex::Real &a_time,
                        std::unique_ptr<AdvanceDiffData> &diffData,
                        std::unique_ptr<AdvanceAdvData> &advData);

   void nonLinearResidual(const amrex::Real &a_dt,
                          const amrex::Vector<amrex::MultiFab *> &a_nlstate,
                          const amrex::Vector<amrex::MultiFab *> &a_nlresid,
                          int updateScaling,
                          int updatePrecond);

   void getAdvectionTerm(const amrex::Vector<const amrex::MultiFab*> &a_nE,
                         const amrex::Vector<amrex::MultiFab*> &a_advTerm,
                         const amrex::Vector<amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM>> &a_gPhiVCur);

   void getAdvectionFluxes(int lev,
                           const amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                           const amrex::MultiFab &a_nE,
                           const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> &a_ueff,
                           amrex::BCRec bcrec);

   void getAdvectionFluxesMOL(int lev,
                              const amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                              const amrex::MultiFab &a_nE,
                              const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> &a_ueff,
                              amrex::BCRec bcrec);

   void ionDriftVelocity(std::unique_ptr<AdvanceAdvData> &advData);
   void ionDriftAddUmac(int lev, std::unique_ptr<AdvanceAdvData> &advData);

   void computeInstantaneousReactionRateEF(int lev,
                                           const PeleLM::TimeStamp &a_time,
                                           amrex::MultiFab* a_I_R);

   amrex::BCRec
   hackBCChargedParticle(const amrex::Real &zk,
                         const amrex::BCRec &bc_in);

   amrex::Real estEFIonsDt(const PeleLM::TimeStamp &a_time);

   void addLorentzVelForces(int lev,
                            const amrex::Box &bx,
                            const amrex::Real &a_time,
                            amrex::Array4<amrex::Real> const& force,
                            amrex::Array4<const amrex::Real> const& rhoY,
                            amrex::Array4<const amrex::Real> const& phiV,
                            amrex::Array4<const amrex::Real> const& nE);

   void nlSolveNorm(const amrex::Vector<amrex::MultiFab*> &a_MF,
                    amrex::Real &norm);

   int testExitNewton(int newtonIter,
                      const amrex::Real &max_res,
                      const amrex::Real &norm_newtonDir); 

   void updateNLState(const amrex::Vector<amrex::MultiFab*> &a_update);

   void incrementElectronForcing(int a_sstep,
                                 std::unique_ptr<AdvanceAdvData> &advData);

   void getNLStateScaling(amrex::Real &nEScale, amrex::Real &phiVScale);
   void getNLResidScaling(amrex::Real &nEScale, amrex::Real &phiVScale);
   void scaleNLState(const amrex::Real &nEScale, const amrex::Real &phiVScale);
   void scaleNLResid(const amrex::Vector<amrex::MultiFab *> &a_resid, const amrex::Real &nEScale, const amrex::Real &phiVScale);

   // get upwinded edge-centered MF from cc MF and ec Vel
   amrex::Array<amrex::MultiFab,AMREX_SPACEDIM>
   getUpwindedEdge(int lev, int edge_comp, int ncomp,
                   amrex::Vector<amrex::BCRec> a_bcrec,
                   const amrex::MultiFab &ccMF,
                   const amrex::Array<const amrex::MultiFab*,AMREX_SPACEDIM> &ecVel);

   void jTimesV(const amrex::Vector<amrex::MultiFab*> &a_x,
                const amrex::Vector<amrex::MultiFab*> &a_Ax);

   void setUpPrecond(const amrex::Real &a_dt,
                     const amrex::Vector<const amrex::MultiFab*> &a_nE);
   void applyPrecond(const amrex::Vector<amrex::MultiFab*> &a_v,
                     const amrex::Vector<amrex::MultiFab*> &a_Pv);

   void ionsBalance();
   //-----------------------------------------------------------------------------
#endif


   //-----------------------------------------------------------------------------
   // Accessors
   std::string runMode() {
      return m_run_mode;
   }

   int nCompIR() {
#ifdef PELE_USE_EFIELD
      return NUM_SPECIES+1;
#else
      return NUM_SPECIES;
#endif
   }

   int nCompForcing() {
#ifdef PELE_USE_EFIELD
      return NUM_SPECIES+2;
#else
      return NUM_SPECIES+1;
#endif
   }

   amrex::FabFactory<amrex::FArrayBox> const&
   Factory (int lev) const noexcept { return *m_factory[lev]; }

   LevelData* getLevelDataPtr(int lev, const PeleLM::TimeStamp &a_time, int useUmac = 0);
   LevelDataReact* getLevelDataReactPtr(int lev);

   amrex::Real getTime(int lev, const PeleLM::TimeStamp &a_time) const {
      AMREX_ASSERT(a_time==AmrOldTime || a_time==AmrNewTime || a_time==AmrHalfTime);
      if ( a_time == AmrOldTime ) {
         return m_t_old[lev];
      } else if ( a_time == AmrNewTime ) {
         return m_t_new[lev];
      } else {
         return 0.5 * ( m_t_old[lev] + m_t_new[lev] );
      }
   }

   amrex::Vector<std::unique_ptr<amrex::MultiFab> > getStateVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<std::unique_ptr<amrex::MultiFab> > getVelocityVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<std::unique_ptr<amrex::MultiFab> > getDensityVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<std::unique_ptr<amrex::MultiFab> > getSpeciesVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<std::unique_ptr<amrex::MultiFab> > getRhoHVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<std::unique_ptr<amrex::MultiFab> > getTempVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<amrex::MultiFab* > getDivUVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<amrex::MultiFab* > getDiffusivityVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<amrex::MultiFab* > getViscosityVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<amrex::MultiFab* > getIRVect();
#ifdef PELE_USE_EFIELD
   amrex::Vector<std::unique_ptr<amrex::MultiFab> > getPhiVVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<std::unique_ptr<amrex::MultiFab> > getnEVect(const PeleLM::TimeStamp &a_time);
   amrex::Vector<amrex::MultiFab* > getnEDiffusivityVect(const PeleLM::TimeStamp &a_time);
#endif

   amrex::Vector<amrex::BCRec> fetchBCRecArray(int scomp, int ncomp);
   amrex::Vector<int>          fetchAdvTypeArray(int scomp, int ncomp);
   amrex::Vector<int>          fetchDiffTypeArray(int scomp, int ncomp);
   DiffusionOp* getDiffusionOp ();
   DiffusionOp* getMCDiffusionOp (int ncomp = 1);
   DiffusionTensorOp* getDiffusionTensorOp ();
   //-----------------------------------------------------------------------------

   //protected:
   // Brief state descriptor
   std::list<std::tuple<int,std::string>> stateComponents;

   // Class level data
   amrex::Vector<std::unique_ptr<LevelData> > m_leveldata_old;
   amrex::Vector<std::unique_ptr<LevelData> > m_leveldata_new;
   amrex::Vector<std::unique_ptr<LevelDataReact> > m_leveldatareact;

   // Temporary holders for halftime data
   std::unique_ptr<LevelData> m_leveldata_floating;
   amrex::Vector<std::unique_ptr<amrex::MultiFab> > m_halfTimeDensity;

   // FAB factory
   amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > > m_factory;

   // External source terms
   amrex::Vector<std::unique_ptr<amrex::MultiFab> >  m_extSource;

   // Covered cells mask
   amrex::Vector<std::unique_ptr<amrex::iMultiFab> > m_coveredMask;
   amrex::Vector<std::unique_ptr<amrex::BoxArray> > m_baChem;
   amrex::Vector<std::unique_ptr<amrex::DistributionMapping> > m_dmapChem;
   amrex::Vector<amrex::Vector<int > > m_baChemFlag;
   int m_resetCoveredMask;

   // Times
   amrex::Vector<amrex::Real> m_t_old;
   amrex::Vector<amrex::Real> m_t_new;

   // Boundary Conditions
   amrex::BCRec  m_phys_bc;                      // Store physical BC from input
#ifdef PELE_USE_EFIELD
   amrex::BCRec  m_phiV_bc;                      // Store separate phiV BC from input
   amrex::BCRec  m_phiV_bcpol;                   // Store separate phiV BC polarity from input
#endif
   amrex::Vector<amrex::BCRec> m_bcrec_state;
   amrex::BCRec                m_bcrec_divu;
   amrex::Vector<amrex::BCRec> m_bcrec_force;

   // Advection/Diffusion types
   amrex::Vector<int> m_AdvTypeState;            // Tag Conservative/NonConservative state components
   amrex::Vector<int> m_DiffTypeState;           // Tag Diffusive/NonDiffusive state components

   // Derived quantities
   PeleLMDeriveList     derive_lst;

   // Derived quantities
   PeleLMDeriveList     evaluate_lst;

   // Tagging/Refinement
   amrex::Vector<amrex::AMRErrorTag> errTags;
   std::string m_initial_grid_file;
   std::string m_regrid_file;
   amrex::Vector<amrex::BoxArray> m_initial_ba;
   amrex::Vector<amrex::BoxArray> m_regrid_ba;

   // Problem parameters
   ProbParm* prob_parm;
   ProbParm* prob_parm_d;

#ifdef PELELM_USE_SOOT
   // Soot parameters
   SootModel* soot_model;
   bool do_soot_solve;
#endif

   // PMF parameters
   pele::physics::PMF::PmfData pmf_data;

   // Transport pointer
   static pele::physics::transport::TransportParams<
             pele::physics::PhysicsType::transport_type> trans_parms;

   // Reactor pointer
   std::string m_chem_integrator;
   std::unique_ptr<pele::physics::reactions::ReactorBase> m_reactor;

   // Turbulence injection
   pele::physics::turbinflow::TurbInflow turb_inflow;

   // Diagnostics
   amrex::Vector<std::unique_ptr<DiagBase>> m_diagnostics; 

   // Linear solvers
   std::unique_ptr<DiffusionOp> m_diffusion_op;
   std::unique_ptr<DiffusionOp> m_mcdiffusion_op;
   std::unique_ptr<DiffusionTensorOp> m_diffusionTensor_op;
   std::unique_ptr<Hydro::MacProjector> macproj;
   int m_macProjNeedReset = 0;
   int m_macProjOldSize = 0;

   int m_verbose = 0;

   // IO options
   std::string m_restart_chkfile{""};
   std::string m_restart_pltfile{""};
   std::string pltfileSource{"LM"};
   std::string m_plot_file{"plt"};
   std::string m_check_file{"chk"};
   int m_derivePlotVarCount = 0;
   amrex::Vector<std::string> m_derivePlotVars;
   int m_plotChemDiag = 0;
   int m_plotHeatRelease = 0;
   int m_plotStateSpec = 0;
   int m_plot_int = 0;
   amrex::Real m_plot_per_approx = -1.;
   amrex::Real m_plot_per_exact = -1.;
   int m_check_int = 0;
   int m_message_int = 10;
   int m_evaluatePlotVarCount = 0;
   int m_plot_grad_p = 1;
   int m_ioDigits = 5;
   amrex::Vector<std::string> m_evaluatePlotVars;
   bool m_write_hdf5_pltfile = false;

   //-----------------------------------------------------------------------------
   // ALGORITHM

   // Time & step
   amrex::Real m_cur_time     = 0.0;
   amrex::Real m_dt           = -1.0;
   amrex::Real m_cfl          = 0.7;
   amrex::Real m_fixed_dt     = -1.0;
   amrex::Real m_init_dt      = -1.0;
   amrex::Real m_min_dt       = 1.0e-12;
   amrex::Real m_dtshrink     = 1.0;
   amrex::Real m_prev_dt      = -1.0;
   amrex::Real m_dtChangeMax  = 1.1;
   amrex::Real m_stop_time    = -1.0;
   amrex::Real m_max_dt       = 1e12;
   int m_nstep = -1; 
   int m_max_step = -1;

   int m_regrid_int = -1;

   // Switch Evolve/Evaluate
   std::string m_run_mode = "normal";

   // Incompressible ?
   int m_incompressible = 0;
   amrex::Real m_rho = -1.0;
   amrex::Real m_mu  = -1.0;

   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> m_gravity{AMREX_D_DECL(0.0,0.0,0.0)};
   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> m_background_gp{AMREX_D_DECL(0.0,0.0,0.0)};

   // DivU ? default ON when not incompressible
   int m_has_divu = 1;

   // DivU dt
   int m_divu_checkFlag = 1;
   amrex::Real m_divu_dtFactor = 0.5;
   amrex::Real m_divu_rhoMin = 0.1;  

   // Aux data size
   int m_nAux = 0;

   // R-z
   int m_rz_correction = 0;

   // cc->ec average
   int m_harm_avg_cen2edge = 0;

   // wbar diffusion term
   int m_use_wbar = 1;

   // switch on/off different processes
   // TODO: actually implement that
   int m_do_react = 1;
   int m_do_diff = 1;
   int m_do_adv = 1;

   // Misc.
   int m_floor_species = 0;
   int m_checkMem = 0;

   // SDC
   int m_nSDCmax = 1;
   int m_sdcIter = 0;

   // DeltaT iterations
   int m_deltaT_verbose = 0;
   int m_deltaTIterMax = 10;
   amrex::Real m_deltaT_norm_max = 1.0e-10;
   int m_crashOnDeltaTFail = 1;

   // Pressure
   amrex::Real m_dpdtFactor = 1.0;
   amrex::Real m_pOld = -1.0;
   amrex::Real m_pNew = -1.0;
   amrex::Real m_dp0dt = 0.0;
   amrex::Real m_uncoveredVol = -1.0;
   int m_closed_chamber = 0;

   // Godunov is default
   std::string m_predict_advection_key = "Godunov_PLM";
   std::string m_advection_key = "Godunov_PLM";
   std::string m_predict_advection_type;
   std::string m_advection_type;
   int m_Godunov_ppm = 1;
   int m_Godunov_ForceInTrans = 0;

   // Temporals
   int m_do_temporals = 0;
   int m_temp_int = 5;
   int m_do_extremas = 0;
   int m_do_massBalance = 0;
   int m_do_energyBalance = 0;
   int m_do_speciesBalance = 0;
   amrex::Real m_massOld;
   amrex::Real m_massNew;
   amrex::Real m_RhoHOld;
   amrex::Real m_RhoHNew;
   amrex::Array<amrex::Real,NUM_SPECIES> m_RhoYOld;
   amrex::Array<amrex::Real,NUM_SPECIES> m_RhoYNew;
   amrex::Array<amrex::Real,2*AMREX_SPACEDIM> m_domainMassFlux;
   amrex::Array<amrex::Real,2*AMREX_SPACEDIM> m_domainRhoHFlux;
   amrex::Array<amrex::Real,2*AMREX_SPACEDIM*NUM_SPECIES> m_domainRhoYFlux;
   amrex::Array<amrex::Real,2*AMREX_SPACEDIM> m_domainUmacFlux;
   
   std::ofstream tmpStateFile;
   std::ofstream tmpExtremasFile;
   std::ofstream tmpMassFile;
   std::ofstream tmpSpecFile;

   // Number of ghost cells
#ifdef AMREX_USE_EB
   int m_nGrowState = 4;
#else
   int m_nGrowState = 3;
#endif
   int m_nGrowAdv = 1;
   int m_nGrowMAC = 1;

   // Initialization
   int m_do_init_proj = 1;
   int m_numDivuIter = 1;
   int m_init_iter = 3;

   // Chemistry
   int m_skipInstantRR = 0;
   int m_plot_react = 1;

   // Typical values
   int m_resetTypValInt = 10;
   int m_useTypValChem = 1;
   amrex::Real m_typicalYvalMin = 1.0e-10;
   amrex::Vector<amrex::Real> typical_values;

#ifdef AMREX_USE_EB
   // EB
   std::string m_EB_refine_type = "Static";
   int m_EB_refine_LevMax = -1;
   int m_EB_refine_LevMin = -1;
   int m_EB_refine_LevAdapt = -1;
   int m_signDistNeeded = 0;
   int m_EB_generate_max_level = -1;
   amrex::Real m_derefineEBBuffer = 3.0;
   std::unique_ptr<amrex::MultiFab> m_signedDist0 = nullptr;
   amrex::Vector<amrex::Real> coveredState_h;
   amrex::Gpu::DeviceVector<amrex::Real> coveredState_d;
   std::string m_adv_redist_type = "StateRedist";
   std::string m_diff_redist_type = "FluxRedist";
#endif
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Combustion physics
   int fuelID = -1;

   // Mixture fraction
   amrex::Array<amrex::Real, 4> Beta_mix;
   amrex::Array<amrex::Real, NUM_SPECIES> spec_Bilger_fact;
   amrex::Real Zfu = -1.0;
   amrex::Real Zox = -1.0;

   // Progress variable
   // C is a linear combination of Ys and T
   amrex::GpuArray<amrex::Real,NUM_SPECIES+1> m_Cweights;
   amrex::Real m_C0 = -1.0;
   amrex::Real m_C1 = -1.0;
   int m_Crevert = 0;           // Flip the definition with 1 - ∑_k wgt*v_k
   //-----------------------------------------------------------------------------

   //-----------------------------------------------------------------------------
   // Linear Solvers

   // Nodal projection
   int m_nodal_mg_max_coarsening_level = 100;
   amrex::Real m_nodal_mg_rtol = 1.0e-11;   
   amrex::Real m_nodal_mg_atol = 1.0e-14;
   std::string m_hypre_namespace_nodal = "nodal_hypre";

   // MAC projection
   int m_mac_mg_verbose = 0;
   int m_mac_mg_max_coarsening_level = 100;
   int m_mac_max_order = 4;
   amrex::Real m_mac_mg_rtol = 1.0e-11;   
   amrex::Real m_mac_mg_atol = 1.0e-14;
   std::string m_hypre_namespace_mac = "mac_hypre";
   //-----------------------------------------------------------------------------
};

void writeBuildInfo();

#endif
