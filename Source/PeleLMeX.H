#ifndef PeleLMeX_H
#define PeleLMeX_H

// PeleLMeX includes
#include "PeleLMeX_UserKeys.H"
#include "PeleLMeX_Index.H"
#include "PeleLMeX_Derive.H"
#include "pelelmex_prob_parm.H"
#include "PeleLMeX_DiffusionOp.H"
#include "DiagBase.H"
#include "PeleLMeX_FlowControllerData.H"
#include "PeleLMeX_BPatch.H"

#ifdef PELE_USE_EFIELD
#include "PrecondOp.H"
#endif

// PelePhysics lib
#include <mechanism.H>
#include <PelePhysics.H>
#include <PMFData.H>
#include <turbinflow.H>
#include <ReactorBase.H>

// AMReX-Hydro lib
#include <hydro_MacProjector.H>
#include <hydro_NodalProjector.H>
#include <hydro_godunov_ppm.H>

// AMReX lib
#include <AMReX_AmrCore.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_ErrorList.H>
#include <AMReX_VisMF.H>

#ifdef PELE_USE_RADIATION
#include <PeleLMRad.H>
#endif

// Forward declarations
#ifdef PELE_USE_SPRAY
class SprayParticleContainer;
#endif

#ifdef PELE_USE_SOOT
class SootModel;
#endif

class BPatch;

const std::string PrettyLine = " " + std::string(78, '=') + "\n";

class PeleLM : public amrex::AmrCore
{

public:
  enum TimeStamp { AmrOldTime, AmrHalfTime, AmrNewTime };

  // constructor
  PeleLM();

  // destructor
  ~PeleLM() override;

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
  void regrid(int lbase, amrex::Real time, bool initial = false) override;

  void MakeNewLevelFromScratch(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm) override;

  void ErrorEst(
    int lev, amrex::TagBoxArray& tags, amrex::Real time, int ng) override;

  void MakeNewLevelFromCoarse(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm) override;

  void RemakeLevel(
    int lev,
    amrex::Real time,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm) override;

  void ClearLevel(int lev) override;
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
  void
  readGridFile(std::string grid_file, amrex::Vector<amrex::BoxArray>& input_ba);

#ifdef PELE_USE_SPRAY
  // ReadSprayParameters
  void SprayReadParameters();

  // SpraySetup
  void SpraySetup();

  // SprayCreateData
  void SprayCreateData();

  // SprayInit
  void SprayInit();
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
  /**
   * \brief Data container for the state data and long lived derived data
   */
  struct LevelData
  {
    LevelData() = default;
    LevelData(
      amrex::BoxArray const& ba,
      amrex::DistributionMapping const& dm,
      amrex::FabFactory<amrex::FArrayBox> const& factory,
      int a_incompressible,
      int a_has_divu,
      int a_nAux,
      int a_nGrowState,
      int a_use_soret,
      int a_do_les);

    // cell-centered state multifabs
    amrex::MultiFab
      state; // Stqte data: velocity, density, rhoYs, rhoH, Temp, RhoRT
    amrex::MultiFab auxiliaries; // Auxiliary variables (passive scalars and
                                 // others) (dim:m_nAux)
    amrex::MultiFab gp;          // pressure gradient (dim:AMREX_SPACEDIM)
    amrex::MultiFab divu;        // Velocity divergence constraint

    // node-centered state multifabs
    amrex::MultiFab press; // nodal pressure (dim:1)

    // cell-centered transport multifabs
    amrex::MultiFab visc_cc; // Viscosity (dim:1)
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>
      visc_turb_fc; // Turbulent Viscosity (dim:)
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>
      lambda_turb_fc;        // Turbulent thermal conductivity (dim:)
    amrex::MultiFab diff_cc; // Diffusivity (dim:NUM_SPECIES+2)
#ifdef PELE_USE_EFIELD
    amrex::MultiFab diffE_cc; // Electron diffusivity (dim:1)
    amrex::MultiFab mobE_cc;  // Electron mobility (dim:1)
    amrex::MultiFab mob_cc;   // Species mobility (dim:NUM_IONS)
#endif
  };

  /**
   * \brief Data container for reaction data
   */
  struct LevelDataReact
  {
    LevelDataReact() = default;
    LevelDataReact(
      const amrex::BoxArray& ba,
      const amrex::DistributionMapping& dm,
      const amrex::FabFactory<amrex::FArrayBox>& factory);
    amrex::MultiFab I_R;    // Species reaction rates
    amrex::MultiFab functC; // Implicit integrator function call count
#ifdef PELE_USE_EFIELD
    amrex::MultiFab I_RnE; // Electron number density reaction term
#endif
  };
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  /**
   * \brief Advance function data container for diffusion term
   */
  struct AdvanceDiffData
  {
    AdvanceDiffData() = default;
    AdvanceDiffData(
      int a_finestLevel,
      const amrex::Vector<amrex::BoxArray>& ba,
      const amrex::Vector<amrex::DistributionMapping>& dm,
      const amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>>>&
        factory,
      int a_nGrowAdv,
      int a_use_wbar,
      int a_use_soret,
      int is_init = 0);
    amrex::Vector<amrex::MultiFab> Dn;   // Diffusion term t^n
    amrex::Vector<amrex::MultiFab> Dnp1; // Diffusion term t^(n+1,k)
    amrex::Vector<amrex::MultiFab> Dhat; // Diffusion term t^(n+1,k+1)

    amrex::Vector<amrex::MultiFab> Dwbar; // Diffusion term of Wbar
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
      wbar_fluxes;                     // Wbar flux correction
    amrex::Vector<amrex::MultiFab> DT; // Diffusion term of T (soret flux)
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
      soret_fluxes; // Soret flux correction
  };
  /**
   * \brief Advance function data container for advection term
   */
  struct AdvanceAdvData
  {
    AdvanceAdvData() = default;
    AdvanceAdvData(
      int a_finestLevel,
      const amrex::Vector<amrex::BoxArray>& ba,
      const amrex::Vector<amrex::DistributionMapping>& dm,
      const amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>>>&
        factory,
      int a_incompressible,
      int a_nGrowAdv,
      int a_nGrowMAC);
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
      umac;                              // MAC face velocity
    amrex::Vector<amrex::MultiFab> AofS; // Scalar advection term at t^(n+1/2)
    amrex::Vector<amrex::MultiFab> chi;  // Thermodynamic constraint
    amrex::Vector<amrex::MultiFab>
      Forcing; // Scalar forcing for both advection and diffusion
    amrex::Vector<amrex::MultiFab> mac_divu; // divu used in MAC projection
#ifdef PELE_USE_EFIELD
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>
      uDrift; // ions drift face velocity
#endif
  };
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // INIT

  /**
   * \brief Top level initialization function
   */
  void initData();

  /**
   * \brief Fill the initial level data container
   * using the user defined init function
   */
  void initLevelData(int lev);

  /**
   * \brief Fill the initial level data container
   * using the data from a pltfile
   * \param a_pltFile path to PeleLMeX plot file
   */
  void initLevelDataFromPlt(int a_lev, const std::string& a_dataPltFile);

  /**
   * \brief Project the initial solution velocity field,
   * initializing the reaction term and providing a first
   * estimate of the time step size
   */
  void projectInitSolution();

  /**
   * \brief Iterate the advance function to provide an
   * initial pressure field.
   */
  void initialIterations();

  // initFromGridFile
  void InitFromGridFile(amrex::Real time);

  // regridFromFile
  void regridFromGridFile(int lbase, amrex::Real time, bool initial);

  /**
   * \brief Check the consistency of the run
   * parameters.
   */
  void checkRunParams();

  //-----------------------------------------------------------------------------
  // REGRID / LOAD BALANCE

  /**
   * \brief Compute load balancing cost estimate on a given level
   * filling the PeleLM class LayoutData with default method
   * \param a_lev target level
   */
  void computeCosts(int a_lev);

  /**
   * \brief Compute load balancing cost estimate on a given level
   * \param a_lev target level
   * \param a_costs LayoutData holding the cost for each box
   * \param a_costMethod cost method
   */
  void computeCosts(
    int a_lev, amrex::LayoutData<amrex::Real>& a_costs, int a_costMethod);

  /**
   * \brief Create/update the DMap used for chemistry on all levels
   */
  void loadBalanceChem();

  /**
   * \brief Create/update the DMap used for chemistry on a given level
   * \param a_lev level of interest
   */
  void loadBalanceChemLev(int a_lev);

  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // NODAL PROJECTION

  /**
   * \brief Initial velocity projection. The pressure and
   * pressure gradient are reset to zero after the projection.
   * Sigma is 1/rho instead of dt/rho.
   */
  void initialProjection();

  /**
   * \brief Initial velocity projection considering zero
   * velocity divergence (no cell-centered RHS) to establish
   * the isostatic pressure when using gravity
   * Sigma is 1/rho instead of dt/rho.
   */
  void initialPressProjection();

  /**
   * \brief Advance function velocity projection.
   * \param is_init flag initial iteration calls
   * \param a_rhoTime time used to fill sigma
   * \param a_dt time time step size
   */
  void velocityProjection(
    int is_init, const PeleLM::TimeStamp& a_rhoTime, const amrex::Real& a_dt);

  /**
   * \brief Actual nodal projection call function.
   * \param a_vel vector of velocity MultiFabs
   * \param a_sigma vector of variable coefficient
   * \param rhs_cc vector of cell-centered projection RHS (can be empty)
   * \param rhs_cc vector of node-centered projection RHS (can be empty)
   * \param increment_gp flag incremental projection (where vel is U^{np1*} -
   * U^{n}) \param scaling_factor used for constant coefficient projection
   */
  void doNodalProject(
    const amrex::Vector<amrex::MultiFab*>& a_vel,
    const amrex::Vector<amrex::MultiFab*>& a_sigma,
    const amrex::Vector<amrex::MultiFab*>& rhs_cc,
    const amrex::Vector<const amrex::MultiFab*>& rhs_nd,
    int incremental,
    amrex::Real scaling_factor);

  /**
   * \brief For 2D-RZ, scale multifab components by radius
   * including 1 ghost cell on Dirichlet BC
   * \param a_lev level index
   * \param a_mf  MultiFab to act upon
   */
  void scaleProj_RZ(int a_lev, amrex::MultiFab& a_mf);

  /**
   * \brief For 2D-RZ, unscale multifab components by radius
   * \param a_lev level index
   * \param a_mf  MultiFab to act upon
   */
  void unscaleProj_RZ(int a_lev, amrex::MultiFab& a_mf);
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // MAC PROJECTION

  /**
   * \brief Re-create the MAC projector (reset the std::unique_ptr)
   * if the circonstances require it
   */
  void resetMacProjector();

  /**
   * \brief Use Godunov to obtain face-centered velocities
   * from cell-centered ones
   * \param advData container for the face-centered velocity umac
   */
  void predictVelocity(std::unique_ptr<AdvanceAdvData>& advData);

  /**
   * \brief Create MAC projection RHS by fillpatching divU
   *  between old and new time
   * \param advData container for mac_divu
   */
  void createMACRHS(std::unique_ptr<AdvanceAdvData>& advData);

  /**
   * \brief Create MAC projection RHS by fillpatching divU
   *  between old and new time
   * \param a_sdcIter current SDC iteration index
   * \param a_time time stamp to access thermodynamic pressure
   * \param advData container for chi and incremented mac_divu
   */
  void addChiIncrement(
    int a_sdcIter,
    const PeleLM::TimeStamp& a_time,
    std::unique_ptr<AdvanceAdvData>& advData);

  /**
   * \brief Perform MAC projection to obtained divergence constrained
   * face-centered (MAC) velocities
   * \param a_time time stamp to access the density used in the projection
   * \param advData container for the face-centered velocities umac
   * \param a_divu velocity divergence constraint (can be empty)
   */
  void macProject(
    const PeleLM::TimeStamp& a_time,
    std::unique_ptr<AdvanceAdvData>& advData,
    const amrex::Vector<amrex::MultiFab*>& a_divu);

  /**
   * \brief Fill the projected face velocities ghost faces, ensuring
   * that the ghost cell velocities are constrained by the ghost cell
   * divergence constraint
   * \param a_lev target level
   * \param a_nGrow number of ghost faces
   * \param crse_geom pointer to the geometry of the next coarser level
   * \param fine_geom pointer to the geometry of the current level
   * \param u_mac_crse face-centered umac of the next coarser level
   * \param u_mac_fine face-centered umac of the current level (those updated by
   * the function) \param divu pointer to the current level divU \param
   * crse_ratio refinement ratio with the next coarser level
   */
  void create_constrained_umac_grown(
    int a_lev,
    int a_nGrow,
    const amrex::Geometry* crse_geom,
    const amrex::Geometry* fine_geom,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> u_mac_crse,
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> u_mac_fine,
    const amrex::MultiFab* divu,
    const amrex::IntVect& crse_ratio);

  /**
   * \brief Setup the linear operator boundary condition for the
   * MAC projection in a given direction
   * \param a_side an AMReX orientation object
   */
  amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>
  getMACProjectionBC(amrex::Orientation::Side a_side);
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // DIFFUSION

  /**
   * \brief Compute the cell-centered viscosity on all levels
   * filling the specified container of Old or New LevelData
   * \param a_time either Old or New stamp
   */
  void calcViscosity(const PeleLM::TimeStamp& a_time);

  /**
   * \brief Compute the face-centered turbulent viscosity on all levels
   * filling the specified container of Old or New LevelData
   * \param a_time either Old or New stamp
   */
  void calcTurbViscosity(const PeleLM::TimeStamp& a_time);

  /**
   * \brief Compute the cell-centereds diffusivity on all levels
   * filling the specified container of Old or New LevelData
   * \param a_time either Old or New stamp
   */
  void calcDiffusivity(const PeleLM::TimeStamp& a_time);

  // get edge-centered diffusivity on a per level / per comp basis
  /**
   * \brief Compute face-averaged diffusivity (or else) from
   * cell-centered data, accounting for Dirichlet BCs, allowing
   * to call user-defined zero_visc() function and possibly adding
   * a turbulent contribution
   * \param lev target level
   * \param beta_comp first component to use in a_diff_cc MultiFab
   * \param ncomp number of components of the outgoing data
   * \param doZeroVisc zero_visc() call trigger
   * \param a_bcrec boundary conditions for the required components
   * \param a_diff_cc cell-centered data
   * \param addTurbContrib turbulent contribution increment trigger
   */
  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> getDiffusivity(
    int lev,
    int beta_comp,
    int ncomp,
    int doZeroVisc,
    amrex::Vector<amrex::BCRec> bcrec,
    const amrex::MultiFab& beta_cc,
    int addTurbContrib = 0);

  /**
   * \brief Compute the explicit face-centered diffusion fluxes on all levels
   * for species and enthalpy (NUM_SPECIES+2 components) using Old or
   * New state/coefficient data from the LevelData container
   * \param a_time either Old or New stamp
   * \param a_fluxes outgoing fluxes container, the first NUM_SPECIES components
   * filled \param a_EBfluxes optional EB fluxes container \param a_wbarfluxes
   * optional Wbar fluxes container \param a_soretfluxes optional Soret fluxes
   * container
   */
  void computeDifferentialDiffusionFluxes(
    const PeleLM::TimeStamp& a_time,
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_fluxes,
    const amrex::Vector<amrex::MultiFab*>& a_EBfluxes,
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_wbarfluxes,
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_soretfluxes);

  /**
   * \brief Compute the explicit cell-centered diffusion term on all levels
   * by computing the divergence of the diffusion fluxes at Old or New time,
   * filling the proper container in diffData
   * \param a_time either Old or New stamp
   * \param diffData container for the returning diffusion terms
   * \param is_init initial computation flag
   */
  void computeDifferentialDiffusionTerms(
    const PeleLM::TimeStamp& a_time,
    std::unique_ptr<AdvanceDiffData>& diffData,
    int is_init = 0);

  /**
   * \brief Add the Wbar contribution to the face-centered species diffusion
   * fluxes on all levels optionally returning that contribution in a separate
   * container \param a_fluxes diffusion fluxes to be updated \param
   * a_wbarfluxes optional container to return the Wbar contribution \param
   * a_spec species rhoYs state data on all levels \param a_rho density state
   * data on all levels \param a_beta cell-centered species diffusivity on all
   * levels
   */
  void addWbarTerm(
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_spfluxes,
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_spwbarfluxes,
    amrex::Vector<amrex::MultiFab const*> const& a_spec,
    amrex::Vector<amrex::MultiFab const*> const& a_rho,
    amrex::Vector<amrex::MultiFab const*> const& a_beta);

  /**
   * \brief Add the Soret contribution to the face-centered species diffusion
   * fluxes on all levels optionally returning that contribution in a separate
   * container \param a_fluxes diffusion fluxes to be updated \param
   * a_soretfluxes optional container to return the Soret contribution \param
   * a_spec species rhoYs state data on all levels \param a_rho density state
   * data on all levels \param a_temp temperature state data on all levels
   * \param a_beta cell-centered species diffusivity on all levels
   */
  void addSoretTerm(
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_spfluxes,
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_spsoretfluxes,
    amrex::Vector<amrex::MultiFab const*> const& a_spec,
    amrex::Vector<amrex::MultiFab const*> const& a_temp,
    amrex::Vector<amrex::MultiFab const*> const& a_beta);

  /**
   * \brief Modify the species diffusion fluxes by computing a correction
   * velocity, ensuring they sum up to zero. \param a_fluxes diffusion fluxes to
   * be updated \param a_spec species rhoYs state data on all levels
   */
  void adjustSpeciesFluxes(
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_spfluxes,
    amrex::Vector<amrex::MultiFab const*> const& a_spec);

  /**
   * \brief Compute the enthalpy flux due to species differential diffusion,
   * filling the NUM_SPECIES+1 component of a_fluxes
   * \param a_fluxes outgoing fluxes container, filled on the NUM_SPECIES+1
   * component \param a_temp temperature state data on all levels
   */
  void computeSpeciesEnthalpyFlux(
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_fluxes,
    amrex::Vector<amrex::MultiFab const*> const& a_temp);

  /**
   * \brief Implicit diffusion solves
   * \param advData container for the external forcing of the linear solve
   * \param diffData container for the returning implicit diffusion term
   */
  void differentialDiffusionUpdate(
    std::unique_ptr<AdvanceAdvData>& advData,
    std::unique_ptr<AdvanceDiffData>& diffData);

  /**
   * \brief Setup the implicit deltaT linear solve
   * \param a_rhs outgoing RHS for the linear solve
   * \param a_Tsave outgoing temporary container for temperature while solving
   * for deltaT \param a_rhoCp outgoing linear solve coefficient \param advData
   * container for external forcing used to build RHS \param diffData container
   * for diffusion forcing used to build RHS
   */
  void deltaTIter_prepare(
    const amrex::Vector<amrex::MultiFab*>& a_rhs,
    const amrex::Vector<amrex::MultiFab*>& a_Tsave,
    const amrex::Vector<amrex::MultiFab*>& a_rhoCp,
    std::unique_ptr<AdvanceAdvData>& advData,
    std::unique_ptr<AdvanceDiffData>& diffData);

  /**
   * \brief Restore state, update enthalpy and compute fluxes after deltaT
   * linear solve within the deltaT iteration scheme \param a_dtiter deltaT
   * iteration index \param a_fluxes diffusion fluxes, with updated Fourier and
   * differential species fluxes \param a_ebfluxes optional temperature fluxes
   * at the EB \param a_Tsave temporary container to restore state temperature
   * \param diffData container for diffusion forcing used to build RHS
   * \param a_deltaT_norm norm of deltaT residual
   */
  void deltaTIter_update(
    int a_dtiter,
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_fluxes,
    const amrex::Vector<amrex::MultiFab*>& a_ebfluxes,
    const amrex::Vector<amrex::MultiFab const*>& a_Tsave,
    std::unique_ptr<AdvanceDiffData>& diffData,
    amrex::Real& a_deltaT_norm);

  /**
   * \brief Compute the viscous force term in the momentum equation using AMReX
   * TensorOp
   * \param a_time state time for velocity, density and viscosity, either Old or
   * New \param a_divtau outgoing divTau term \param use_density flag to scale
   * divTau by 1/rho \param scale prefactor of divTau, used to weight by num.
   * scheme property (e.g. 0.5 for Crank Nicholson)
   */
  void computeDivTau(
    const PeleLM::TimeStamp& a_time,
    const amrex::Vector<amrex::MultiFab*>& a_divtau,
    int use_density,
    amrex::Real scale = 1.0);

  /**
   * \brief Implicit diffusion of the of velocity field using AMReX
   * TensorOp. Always act on New time data and include Crank Nicholson factor
   */
  void diffuseVelocity();

  /**
   * \brief Return LinOpBC for the scalar diffusion operator
   * \param a_side AMReX orientation
   * \param a_bc BCRec for the scalar of interest
   */
  amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> getDiffusionLinOpBC(
    amrex::Orientation::Side a_side, const amrex::BCRec& a_bc);

  /**
   * \brief Return AMReX_SPACEDIM LinOpBC vector for the velocity tensor
   * operator \param a_side AMReX orientation \param a_bc BCRec for all the
   * velocity components
   */
  amrex::Vector<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>>
  getDiffusionTensorOpBC(
    amrex::Orientation::Side a_side, const amrex::Vector<amrex::BCRec> a_bc);
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // REACTION

  /**
   * \brief Top-level reaction solve function, performing the chemistry
   * integration on all the AMR levels. \param advData container for Forcing,
   * the advection/diffusion component of the advance function in the SDC
   * context
   */
  void advanceChemistry(std::unique_ptr<AdvanceAdvData>& advData);

  /**
   * \brief Performing the chemistry integration on a given level, for a fixed
   * step length and relying on AmrCore BoxArray/DMap, not masking fine-covered
   * regions \param lev level of interest \param a_dt integration length \param
   * a_extForcing advection/diffusion forcing
   */
  void advanceChemistry(
    int lev, const amrex::Real& a_dt, amrex::MultiFab& a_extForcing);

  /**
   * \brief Performing the chemistry integration on a given level, for a fixed
   * step length and using the BoxArray/DMap masking fine-covered region and
   * specifically load balanced for reaction
   * \param lev level of interest
   * \param a_dt integration length
   * \param a_extForcing advection/diffusion forcing
   */
  void advanceChemistryBAChem(
    int lev, const amrex::Real& a_dt, amrex::MultiFab& a_extForcing);

  /**
   * \brief Top-level instantaneous reaction rate function, acting on all levels
   * \param a_I_R outgoing multi-level container inst. RR container
   * \param a_time state time used to compute RR, either Old or New
   */
  void computeInstantaneousReactionRate(
    const amrex::Vector<amrex::MultiFab*>& a_I_R,
    const PeleLM::TimeStamp& a_time);

  /**
   * \brief Compute instantaneous reaction rate on a given level
   * \param lev level of interest
   * \param a_time state time used to compute RR, either Old or New
   * \param a_I_R outgoing MultiFab for inst. RR container
   */
  void computeInstantaneousReactionRate(
    int lev, const PeleLM::TimeStamp& a_time, amrex::MultiFab* a_I_R);

  /**
   * \brief Compute the heat release rate on a given level, using the
   * LevelDataReact I_R data
   * \param lev level of interest
   * \param a_HR outgoing MultiFab
   */
  void getHeatRelease(int lev, amrex::MultiFab* a_HR);
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // Time step size

  /**
   * \brief Compute the time step across all levels and processes.
   * \param is_init flag to trigger init-only dt factors
   * \param a_time state time to work on, either Old or New
   */
  amrex::Real computeDt(int is_init, const PeleLM::TimeStamp& a_time);

  /**
   * \brief Compute the CFL step size
   * \param a_time state time to work on, either Old or New
   */
  amrex::Real estConvectiveDt(const PeleLM::TimeStamp& a_time);

  /**
   * \brief Compute the divU step size
   * \param a_time state time to work on, either Old or New
   */
  amrex::Real estDivUDt(const PeleLM::TimeStamp& a_time);
  void checkDt(const PeleLM::TimeStamp& a_time, const amrex::Real& a_dt);
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // FORCES

  void getVelForces(
    const TimeStamp& a_time,
    const amrex::Vector<amrex::MultiFab*>& a_divTau,
    const amrex::Vector<amrex::MultiFab*>& a_velForce,
    int nGrowForce,
    int add_gradP);

  void getVelForces(
    const TimeStamp& a_time,
    int lev,
    amrex::MultiFab* a_divTau,
    amrex::MultiFab* a_velForce,
    int add_gradP);

  void getVelForces(
    int lev,
    const amrex::Box& bx,
    const amrex::Real& a_time,
    amrex::Array4<amrex::Real> const& force,
    amrex::Array4<const amrex::Real> const& vel,
    amrex::Array4<const amrex::Real> const& rho,
    amrex::Array4<const amrex::Real> const& rhoY,
    amrex::Array4<const amrex::Real> const& rhoh,
    amrex::Array4<const amrex::Real> const& temp,
    amrex::Array4<const amrex::Real> const& extMom,
    amrex::Array4<const amrex::Real> const& extRho);
  //-----------------------------------------------------------------------------

#ifdef PELE_USE_SPRAY
  //-----------------------------------------------------------------------------
  // Sprays
  void setupVirtualParticles(const int level);
  void removeVirtualParticles(const int level);
  void setupGhostParticles(const int ngrow, const int level);
  void removeGhostParticles(const int level);
  amrex::Real SprayEstDt();
  void SprayMKD(const amrex::Real time, const amrex::Real dt);
  void
  SprayMKDLevel(const int level, const amrex::Real time, const amrex::Real dt);
  void SprayInjectRedist();
  void SprayPostRegrid();
  void SpraySetState(const amrex::Real& a_flow_dt);
  void SprayAddSource(const int level);
  static int num_spray_src;
  static bool do_spray_particles;
  static std::unique_ptr<SprayParticleContainer> SprayPC;
  static std::unique_ptr<SprayParticleContainer> VirtPC;
  static std::unique_ptr<SprayParticleContainer> GhostPC;
  // State used for spray interpolation
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_spraysource;
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_spraystate;
  //-----------------------------------------------------------------------------
#endif
#ifdef PELE_USE_SOOT
  //-----------------------------------------------------------------------------
  // Soot

  // SetSootIndx
  void setSootIndx();

  // CleanupSootModel
  void cleanupSootModel();

  // AddSootDerivePlotVars
  // void addSootDerivePlotVars(PeleLMDeriveList& derive_lst);

  // Compute the soot source terms
  void computeSootSource(
    const PeleLM::TimeStamp& a_timestamp, const amrex::Real a_dt);

  // Clip soot moments
  void clipSootMoments();

  //-----------------------------------------------------------------------------
#endif

#ifdef PELE_USE_RADIATION
  // Initialize the radiation module
  void RadInit();

  // Compute the radiation source term
  void computeRadSource(const PeleLM::TimeStamp& a_timestamp);
#endif

  //-----------------------------------------------------------------------------
  // EOS

  /**
   * \brief Compute thermodynamic pressure (RhoRT) from state on all levels
   * \param a_time state time to work on, either Old or New
   */
  void setThermoPress(const PeleLM::TimeStamp& a_time);

  /**
   * \brief Compute thermodynamic pressure (RhoRT) from state on a single level
   * \param lev level of interest
   * \param a_time state time to work on, either Old or New
   */
  void setThermoPress(int lev, const PeleLM::TimeStamp& a_time);

  /**
   * \brief Enforce rho = \sum rhoYs on a given level
   * \param lev level of interest
   * \param a_time state time to work on, either Old or New
   */
  void setRhoToSumRhoY(int lev, const PeleLM::TimeStamp& a_time);

  /**
   * \brief Compute temperature from rhoY and rhoH on all levels
   * \param a_time state time to work on, either Old or New
   */
  void setTemperature(const PeleLM::TimeStamp& a_time);

  /**
   * \brief Compute temperature from rhoY and rhoH on a given level
   * \param lev level of interest
   * \param a_time state time to work on, either Old or New
   */
  void setTemperature(int lev, const PeleLM::TimeStamp& a_time);

  /**
   * \brief Compute the divergence constraint from reaction/diffusion
   * terms
   * \param is_init flag for init-only checks
   * \param computeDiff flag to trigger computing the diffusion terms or
   * using the one already precompute and available in diffData
   * \param a_time state time to work with
   * \param diffData container for the diffusion terms
   */
  void calcDivU(
    int is_init,
    int computeDiff,
    int do_avgDown,
    const PeleLM::TimeStamp& a_time,
    std::unique_ptr<AdvanceDiffData>& diffData);

  /**
   * \brief Compute the pressure drift term using RhoRT and
   * state data at a_time, acting on all levels
   * \param a_time state time to work on, either Old or New
   * \param a_dPdt outgoing container
   */
  void calc_dPdt(
    const PeleLM::TimeStamp& a_time,
    const amrex::Vector<amrex::MultiFab*>& a_dPdt);

  /**
   * \brief Compute the pressure drift term using RhoRT and
   * state data at a_time on a given level
   * \param lev level of interest
   * \param a_time state time to work on, either Old or New
   * \param a_dPdt outgoing container
   */
  void
  calc_dPdt(int lev, const PeleLM::TimeStamp& a_time, amrex::MultiFab* a_dPdt);

  /**
   * \brief For the closed chamber algorithm, compute averaged
   * mac_divU and theta, update divU and the spatially averaged
   * thermodynamic pressure. Return the mean divU.
   * \param advData container including mac_divu
   */
  amrex::Real adjustPandDivU(std::unique_ptr<AdvanceAdvData>& advData);
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // I/O
  void WritePlotFile();
  bool writePlotNow() const;
  bool checkMessage(const std::string& a_action) const;
  void WriteCheckPointFile();
  void ReadCheckPointFile();
  bool writeCheckNow() const;
  void WriteJobInfo(const std::string& path) const;
  void WriteHeader(const std::string& name, bool is_checkpoint) const;
  void WriteDebugPlotFile(
    const amrex::Vector<const amrex::MultiFab*>& a_MF,
    const std::string& pltname);
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
  void speciesBalancePatch();
  void initBPatches(amrex::Geometry& a_geom);

  void rhoHBalance();
  void initTemporals(const PeleLM::TimeStamp& a_time = AmrOldTime);
  void writeTemporals();
  void addMassFluxes(
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    const amrex::Geometry& a_geom);
  void addRhoHFluxes(
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    const amrex::Geometry& a_geom);
  void addRhoYFluxes(
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    const amrex::Geometry& a_geom,
    const amrex::Real& a_factor = 1.0);
  void addRhoYFluxesPatch(
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    const amrex::Geometry& a_geom,
    const amrex::Real& a_factor = 1.0);
  void addUmacFluxes(
    std::unique_ptr<AdvanceAdvData>& advData, const amrex::Geometry& a_geom);
  void openTempFile();
  void closeTempFile();
  bool doTemporalsNow() const;
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // BOUNDARY CONDITIONS / FILLPATCH
  amrex::InterpBase* getInterpolator(int a_method = 1) const;

  // Interpolation method for regrid (default to cell cons interp)
  int m_regrid_interp_method = 1;

  // Convert PhysBC into field BCs
  void setBoundaryConditions();

  // FillPatch the entire class state
  void fillPatchState(const PeleLM::TimeStamp& a_time);
  void fillPatchState(int lev, const PeleLM::TimeStamp& a_time);

  // FillPatch the various state components
  void fillPatchDensity(const PeleLM::TimeStamp& a_time);
  void fillPatchSpecies(const PeleLM::TimeStamp& a_time);
  void fillPatchTemp(const PeleLM::TimeStamp& a_time);
#ifdef PELE_USE_EFIELD
  void fillPatchPhiV(const PeleLM::TimeStamp& a_time);
#endif

  // FillPatch and return a MF
  std::unique_ptr<amrex::MultiFab>
  fillPatchState(int lev, amrex::Real a_time, int nGrow);
  std::unique_ptr<amrex::MultiFab>
  fillPatchReact(int lev, amrex::Real a_time, int nGrow);

  // FillPatch state components
  void fillpatch_state(
    int lev, amrex::Real a_time, amrex::MultiFab& a_state, int nGhost);
  void fillpatch_density(
    int lev,
    amrex::Real a_time,
    amrex::MultiFab& a_density,
    int rho_comp,
    int nGhost);
  void fillpatch_species(
    int lev,
    amrex::Real a_time,
    amrex::MultiFab& a_species,
    int rhoY_comp,
    int nGhost);
  void fillpatch_temp(
    int lev,
    amrex::Real a_time,
    amrex::MultiFab& a_temp,
    int temp_comp,
    int nGhost);
  void fillpatch_divu(
    int lev, amrex::Real a_time, amrex::MultiFab& a_divu, int nGhost);
  void fillpatch_gradp(
    int lev, amrex::Real a_time, amrex::MultiFab& a_gp, int nGhost);
  void fillpatch_reaction(
    int lev, amrex::Real a_time, amrex::MultiFab& a_I_R, int nGhost);
  void fillpatch_forces(
    amrex::Real a_time,
    amrex::Vector<amrex::MultiFab*> const& a_force,
    int nGrowForce);
  void fillpatch_chemFunctCall(
    int lev, amrex::Real a_time, amrex::MultiFab& a_fctC, int nGhost);
#ifdef PELE_USE_EFIELD
  void fillpatch_phiV(
    int lev,
    amrex::Real a_time,
    amrex::MultiFab& a_phiV,
    int phiV_comp,
    int nGhost);
  void fillPatchExtrap(
    amrex::Real a_time,
    amrex::Vector<amrex::MultiFab*> const& a_MF,
    int a_nGrow);
  void fillPatchNLnE(
    amrex::Real a_time,
    amrex::Vector<amrex::MultiFab*> const& a_nE,
    int a_nGrow);
  void fillPatchNLphiV(
    amrex::Real a_time,
    amrex::Vector<amrex::MultiFab*> const& a_phiV,
    int a_nGrow);
#endif

  // FillCoarsePatch state components
  void fillcoarsepatch_state(
    int lev, amrex::Real a_time, amrex::MultiFab& a_state, int nGhost);
  void fillcoarsepatch_divu(
    int lev, amrex::Real a_time, amrex::MultiFab& a_divu, int nGhost);
  void fillcoarsepatch_gradp(
    int lev, amrex::Real a_time, amrex::MultiFab& a_gp, int nGhost);
  void fillcoarsepatch_reaction(
    int lev, amrex::Real a_time, amrex::MultiFab& a_I_R, int nGhost);
  void fillcoarsepatch_chemFunctCall(
    int lev, amrex::Real a_time, amrex::MultiFab& a_fctC, int nGhost);

  // Fill physical boundaries
  void setInflowBoundaryVel(
    amrex::MultiFab& a_vel, int lev, PeleLM::TimeStamp a_time);

  void fillTurbInflow(
    amrex::MultiFab& a_vel, int vel_comp, int lev, const amrex::Real a_time);

  // Average down operations
  void averageDownState(const PeleLM::TimeStamp& a_time);
  void averageDownScalars(const PeleLM::TimeStamp& a_time);
  void averageDownVelocity(const PeleLM::TimeStamp& a_time);
  void averageDown(
    const PeleLM::TimeStamp& a_time, const int state_comp, const int ncomp);
  void averageDownReaction();
#ifdef PELE_USE_EFIELD
  void averageDownnE(const PeleLM::TimeStamp& a_time);
  void averageDownPhiV(const PeleLM::TimeStamp& a_time);
#endif
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // UTILS

  template <class K>
  void parseUserKey(
    amrex::ParmParse& a_pp,
    const std::string& a_search_key,
    const K& keyT,
    int& key,
    int keyIndex = -1)
  {
    std::string userkey = "default";
    if (keyIndex < 0) {
      a_pp.query(a_search_key.c_str(), userkey);
    } else {
      int search_key_size = a_pp.countval(a_search_key.c_str());
      if (search_key_size > 0) {
        AMREX_ASSERT(keyIndex < search_key_size);
        amrex::Vector<std::string> userkeyVec(search_key_size);
        a_pp.getarr(a_search_key.c_str(), userkeyVec, keyIndex, 1);
        userkey = userkeyVec[keyIndex];
      }
    }
    std::transform(userkey.begin(), userkey.end(), userkey.begin(), ::tolower);
    if (
      std::find(
        std::begin(keyT.searchKey), std::end(keyT.searchKey), a_search_key) ==
      std::end(keyT.searchKey)) {
      std::string errmsg = "Improper use of parseUserKey: search_key "
                           "incompatible with provided userKey. \n";
      errmsg += "The valid search keys are:\n";
      for (const auto& validKey : keyT.searchKey) {
        errmsg += " - " + validKey + " \n";
      }
      amrex::Abort(errmsg);
    }

    if (keyT.str2int.count(userkey) == 0) {
      std::string errmsg =
        "Invalid key for search_key " + a_search_key + ". \n";
      errmsg += "Valid keys are:\n";
      for (const auto& valid_pair : keyT.str2int) {
        if (valid_pair.first != "default") {
          errmsg += " - " + valid_pair.first + " \n";
        }
      }
      amrex::Abort(errmsg);
    }

    key = keyT.str2int.at(userkey);
  }

  // flux divergence
  void fluxDivergence(
    const amrex::Vector<amrex::MultiFab*>& a_divergence,
    int div_comp,
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_fluxes,
    int flux_comp,
    int ncomp,
    int intensiveFluxes,
    amrex::Real scale);

  void fluxDivergence(
    const amrex::Vector<amrex::MultiFab*>& a_divergence,
    int div_comp,
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_fluxes,
    int flux_comp,
    const amrex::Vector<amrex::MultiFab*>& a_EBfluxes,
    int ebflux_comp,
    int ncomp,
    int intensiveFluxes,
    amrex::Real scale);

  void fluxDivergenceRD(
    const amrex::Vector<const amrex::MultiFab*>& a_state,
    int state_comp,
    const amrex::Vector<amrex::MultiFab*>& a_divergence,
    int div_comp,
    const amrex::Vector<amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_fluxes,
    int flux_comp,
    const amrex::Vector<amrex::MultiFab*>& a_EBfluxes,
    int ebflux_comp,
    int ncomp,
    int intensiveFluxes,
    const amrex::BCRec* state_bc_d,
    const amrex::Real& scale,
    const amrex::Real& a_dt);

  void advFluxDivergence(
    int a_lev,
    amrex::MultiFab& a_divergence,
    int div_comp,
    amrex::MultiFab& a_divu,
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    int flux_comp,
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_faceState,
    int face_comp,
    int ncomp,
    int const* l_conserv_d,
    const amrex::Geometry& a_geom,
    amrex::Real scale,
    bool fluxes_are_area_weighted) const;

  void extFluxDivergenceLevel(
    int lev,
    amrex::MultiFab& a_divergence,
    int div_comp,
    const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    int flux_comp,
    int ncomp,
    amrex::Real scale);

  void intFluxDivergenceLevel(
    int lev,
    amrex::MultiFab& a_divergence,
    int div_comp,
    const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    int flux_comp,
    int ncomp,
    amrex::Real scale);

  void intFluxDivergenceLevelEB(
    int lev,
    amrex::MultiFab& a_divergence,
    int div_comp,
    const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    int flux_comp,
    const amrex::MultiFab* a_EBfluxes,
    int ebflux_comp,
    int ncomp,
    amrex::Real scale);

  void floorSpecies(const PeleLM::TimeStamp& a_time);

  amrex::Real MLNorm0(const amrex::Vector<const amrex::MultiFab*>& a_MF);
  amrex::Vector<amrex::Real> MLNorm0(
    const amrex::Vector<const amrex::MultiFab*>& a_MF,
    int startcomp,
    int ncomp);

  amrex::Real
  MFSum(const amrex::Vector<const amrex::MultiFab*>& a_MF, int comp);
  static amrex::Real
  MFmax(const amrex::MultiFab* a_MF, const amrex::iMultiFab& a_mask, int comp);
  static amrex::Real
  MFmin(const amrex::MultiFab* a_MF, const amrex::iMultiFab& a_mask, int comp);

  amrex::Vector<amrex::Real> MLmax(
    const amrex::Vector<const amrex::MultiFab*>& a_MF, int scomp, int ncomp);

  amrex::Vector<amrex::Real> MLmin(
    const amrex::Vector<const amrex::MultiFab*>& a_MF, int scomp, int ncomp);

  void resetCoveredMask();

  // Active control
  void initActiveControl();
  void activeControl(int is_restart);
  void loadActiveControlHistory();
  void getActiveControlLowT(amrex::Real& a_coft);

  // Derived variable evaluation
  std::unique_ptr<amrex::MultiFab>
  derive(const std::string& a_name, amrex::Real a_time, int lev, int nGrow);
  std::unique_ptr<amrex::MultiFab>
  deriveComp(const std::string& a_name, amrex::Real a_time, int lev, int nGrow);

  // Evaluate function
  void MLevaluate(
    const amrex::Vector<amrex::MultiFab*>& a_MFVec,
    int a_comp,
    int& nComp,
    const std::string& a_var);

  void
  evaluateChemExtForces(const amrex::Vector<amrex::MultiFab*>& a_chemForces);
  void
  evaluateAdvectionTerms(const amrex::Vector<amrex::MultiFab*>& a_advTerms);

  // Lightweight descriptors access
  bool isStateVariable(std::string_view a_name);
  bool isReactVariable(std::string_view a_name);
  int stateVariableIndex(std::string_view a_name);
  int reactVariableIndex(std::string_view a_name);

  // Typical values
  void setTypicalValues(const PeleLM::TimeStamp& a_time, int is_init = 0);
  void updateTypicalValuesChem();

  void checkMemory(const std::string& a_message) const;

  // Mixture fraction & Progress variable
  void initMixtureFraction();
  static void parseComposition(
    amrex::Vector<std::string> compositionIn,
    std::string compositionType,
    amrex::Real* massFrac);
  static void parseVars(
    const amrex::Vector<std::string>& a_varsNames,
    const amrex::Vector<std::string>& a_stringIn,
    amrex::Vector<amrex::Real>& a_rVars);
  void initProgressVariable();
  //-----------------------------------------------------------------------------

#ifdef AMREX_USE_EB
  //-----------------------------------------------------------------------------
  // EB

  void makeEBGeometry();

  int getRestartEBMaxLevel() const;

  void initialRedistribution();

  void initCoveredState();
  void setCoveredState(const PeleLM::TimeStamp& a_time);
  void setCoveredState(int lev, const PeleLM::TimeStamp& a_time);

  void redistributeAofS(
    int a_lev,
    amrex::Real& a_dt,
    amrex::MultiFab& a_tmpDiv,
    int div_comp,
    amrex::MultiFab& a_AofS,
    int aofs_comp,
    amrex::MultiFab& a_state,
    int state_comp,
    int ncomp,
    const amrex::BCRec* d_bc,
    const amrex::Geometry& a_geom) const;

  void redistributeDiff(
    int a_lev,
    const amrex::Real& a_dt,
    amrex::MultiFab& a_tmpDiv,
    int div_comp,
    amrex::MultiFab& a_diff,
    int diff_comp,
    const amrex::MultiFab& a_state,
    int state_comp,
    int ncomp,
    const amrex::BCRec* d_bc,
    const amrex::Geometry& a_geom) const;

  void getCoveredIMask(int a_lev, amrex::iMultiFab& a_imask) const;

  void getEBDistance(int a_lev, amrex::MultiFab& a_signDistLev);

  void
  extendSignedDistance(amrex::MultiFab* a_signDist, amrex::Real a_extendFactor);

  void getEBState(
    int a_lev,
    const amrex::Real& a_time,
    amrex::MultiFab& a_EBstate,
    int stateComp,
    int nComp);

  void getEBDiff(
    int a_lev,
    const TimeStamp& a_time,
    amrex::MultiFab& a_EBDiff,
    int diffComp);

  // static_cast factory into ebfactory
  amrex::EBFArrayBoxFactory const& EBFactory(int lev) const noexcept
  {
    return static_cast<amrex::EBFArrayBoxFactory const&>(*m_factory[lev]);
  }

  void correct_vel_small_cells(
    amrex::Vector<amrex::MultiFab*> const& a_vel,
    amrex::Vector<amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM>> const&
      a_umac);
#endif

  //-----------------------------------------------------------------------------
  // ADVANCE

  void advanceSetup();

  void oneSDC(
    int sdcIter,
    std::unique_ptr<AdvanceAdvData>& advData,
    std::unique_ptr<AdvanceDiffData>& diffData);

  void copyTransportOldToNew();
  void copyStateNewToOld(int nGhost = 0);
  void copyPressNewToOld();
  void copyStateOldToNew(int nGhost = 0);
  void copyDiffusionOldToNew(std::unique_ptr<AdvanceDiffData>& diffData);

  void computeVelocityAdvTerm(std::unique_ptr<AdvanceAdvData>& advData);
  void updateVelocity(std::unique_ptr<AdvanceAdvData>& advData);

  void getScalarAdvForce(
    std::unique_ptr<AdvanceAdvData>& advData,
    std::unique_ptr<AdvanceDiffData>& diffData);
  void computeScalarAdvTerms(std::unique_ptr<AdvanceAdvData>& advData);
  void computePassiveAdvTerms(
    std::unique_ptr<AdvanceAdvData>& advData, int state_comp, int ncomp);
  void updateDensity(std::unique_ptr<AdvanceAdvData>& advData);
  void updateScalarComp(
    std::unique_ptr<AdvanceAdvData>& advData, int state_comp, int ncomp);

  void getScalarDiffForce(
    std::unique_ptr<AdvanceAdvData>& advData,
    std::unique_ptr<AdvanceDiffData>& diffData);

  void getScalarReactForce(std::unique_ptr<AdvanceAdvData>& advData);
  //-----------------------------------------------------------------------------

#ifdef PELE_USE_EFIELD
#include <PeleLMeX_EF.H>
  //-----------------------------------------------------------------------------
  // EFIELD

  void poissonSolveEF(const PeleLM::TimeStamp& a_time);

  void calcEFTransport(const PeleLM::TimeStamp& a_time);

  void initializeElectronNeutral();
  void initializeElectronFromMassFraction();

  void implicitNonLinearSolve(
    int sdcIter,
    const amrex::Real& a_dt,
    std::unique_ptr<AdvanceDiffData>& diffData,
    std::unique_ptr<AdvanceAdvData>& advData);
  void computeBGcharge(
    const amrex::Real& a_time,
    std::unique_ptr<AdvanceDiffData>& diffData,
    std::unique_ptr<AdvanceAdvData>& advData);

  void nonLinearResidual(
    const amrex::Real& a_dt,
    const amrex::Vector<amrex::MultiFab*>& a_nlstate,
    const amrex::Vector<amrex::MultiFab*>& a_nlresid,
    int updateScaling,
    int updatePrecond);

  void getAdvectionTerm(
    const amrex::Vector<const amrex::MultiFab*>& a_nE,
    const amrex::Vector<amrex::MultiFab*>& a_advTerm,
    const amrex::Vector<amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>>&
      a_gPhiVCur);

  void getAdvectionFluxes(
    int lev,
    const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    const amrex::MultiFab& a_nE,
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_ueff,
    amrex::BCRec bcrec);

  void getAdvectionFluxesMOL(
    int lev,
    const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    const amrex::MultiFab& a_nE,
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& a_ueff,
    amrex::BCRec bcrec);

  void ionDriftVelocity(std::unique_ptr<AdvanceAdvData>& advData);
  void ionDriftAddUmac(int lev, std::unique_ptr<AdvanceAdvData>& advData);

  void computeInstantaneousReactionRateEF(
    int lev, const PeleLM::TimeStamp& a_time, amrex::MultiFab* a_I_R);

  amrex::BCRec
  hackBCChargedParticle(const amrex::Real& zk, const amrex::BCRec& bc_in);

  amrex::Real estEFIonsDt(const PeleLM::TimeStamp& a_time);

  void addLorentzVelForces(
    int lev,
    const amrex::Box& bx,
    const amrex::Real& a_time,
    amrex::Array4<amrex::Real> const& force,
    amrex::Array4<const amrex::Real> const& rhoY,
    amrex::Array4<const amrex::Real> const& phiV,
    amrex::Array4<const amrex::Real> const& nE);

  void
  nlSolveNorm(const amrex::Vector<amrex::MultiFab*>& a_MF, amrex::Real& norm);

  int testExitNewton(
    int newtonIter,
    const amrex::Real& max_res,
    const amrex::Real& norm_newtonDir);

  void updateNLState(const amrex::Vector<amrex::MultiFab*>& a_update);

  void incrementElectronForcing(
    int a_sstep, std::unique_ptr<AdvanceAdvData>& advData);

  void getNLStateScaling(amrex::Real& nEScale, amrex::Real& phiVScale);
  void getNLResidScaling(amrex::Real& nEScale, amrex::Real& phiVScale);
  void scaleNLState(const amrex::Real& nEScale, const amrex::Real& phiVScale);
  void scaleNLResid(
    const amrex::Vector<amrex::MultiFab*>& a_resid,
    const amrex::Real& nEScale,
    const amrex::Real& phiVScale);

  // get upwinded edge-centered MF from cc MF and ec Vel
  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> getUpwindedEdge(
    int lev,
    int edge_comp,
    int ncomp,
    amrex::Vector<amrex::BCRec> a_bcrec,
    const amrex::MultiFab& ccMF,
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& ecVel);

  void jTimesV(
    const amrex::Vector<amrex::MultiFab*>& a_x,
    const amrex::Vector<amrex::MultiFab*>& a_Ax);

  void setUpPrecond(
    const amrex::Real& a_dt, const amrex::Vector<const amrex::MultiFab*>& a_nE);
  void applyPrecond(
    const amrex::Vector<amrex::MultiFab*>& a_v,
    const amrex::Vector<amrex::MultiFab*>& a_Pv);

  void ionsBalance();
  //-----------------------------------------------------------------------------
#endif

  //-----------------------------------------------------------------------------
  // Accessors
  std::string runMode() const { return m_run_mode; }

  static int nCompIR()
  {
#ifdef PELE_USE_EFIELD
    return NUM_SPECIES + 1;
#else
    return NUM_SPECIES;
#endif
  }

  static int nCompForcing()
  {
#ifdef PELE_USE_EFIELD
    return NUM_SPECIES + 2;
#else
    return NUM_SPECIES + 1;
#endif
  }

  amrex::FabFactory<amrex::FArrayBox> const& Factory(int lev) const noexcept
  {
    return *m_factory[lev];
  }

  LevelData*
  getLevelDataPtr(int lev, const PeleLM::TimeStamp& a_time, int useUmac = 0);
  LevelDataReact* getLevelDataReactPtr(int lev);

  amrex::Real getTime(int lev, const PeleLM::TimeStamp& a_time) const
  {
    AMREX_ASSERT(
      a_time == AmrOldTime || a_time == AmrNewTime || a_time == AmrHalfTime);
    if (a_time == AmrOldTime) {
      return m_t_old[lev];
    }
    if (a_time == AmrNewTime) {
      return m_t_new[lev];
    }
    return 0.5 * (m_t_old[lev] + m_t_new[lev]);
  }

  amrex::Vector<std::unique_ptr<amrex::MultiFab>>
  getStateVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>>
  getVelocityVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>>
  getDensityVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>>
  getSpeciesVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>>
  getRhoHVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>>
  getTempVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<amrex::MultiFab*> getDivUVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<amrex::MultiFab*>
  getDiffusivityVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<amrex::MultiFab*>
  getViscosityVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<amrex::MultiFab*> getIRVect();
#ifdef PELE_USE_EFIELD
  amrex::Vector<std::unique_ptr<amrex::MultiFab>>
  getPhiVVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<std::unique_ptr<amrex::MultiFab>>
  getnEVect(const PeleLM::TimeStamp& a_time);
  amrex::Vector<amrex::MultiFab*>
  getnEDiffusivityVect(const PeleLM::TimeStamp& a_time);
#endif

  amrex::Vector<amrex::BCRec> fetchBCRecArray(int scomp, int ncomp);
  amrex::Vector<int> fetchAdvTypeArray(int scomp, int ncomp);
  amrex::Vector<int> fetchDiffTypeArray(int scomp, int ncomp);
  DiffusionOp* getDiffusionOp();
  DiffusionOp* getMCDiffusionOp(int ncomp = 1);
  DiffusionTensorOp* getDiffusionTensorOp();
  //-----------------------------------------------------------------------------

  // protected:
  //  Brief state descriptor
  std::list<std::tuple<int, std::string>> stateComponents;
  std::list<std::tuple<int, std::string>> reactComponents;

  // Class level data
  amrex::Vector<std::unique_ptr<LevelData>> m_leveldata_old;
  amrex::Vector<std::unique_ptr<LevelData>> m_leveldata_new;
  amrex::Vector<std::unique_ptr<LevelDataReact>> m_leveldatareact;

  // Temporary holders for halftime data
  std::unique_ptr<LevelData> m_leveldata_floating;
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_halfTimeDensity;

  // FAB factory
  amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox>>> m_factory;

  // External source terms
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_extSource;

  // Fine-covered cells mask
  amrex::Vector<std::unique_ptr<amrex::iMultiFab>> m_coveredMask;
  int m_resetCoveredMask;

  // Chemistry BA & DM
  amrex::Vector<std::unique_ptr<amrex::BoxArray>> m_baChem;
  amrex::Vector<std::unique_ptr<amrex::DistributionMapping>> m_dmapChem;
  amrex::Vector<amrex::Vector<int>> m_baChemFlag;
  amrex::IntVect m_max_grid_size_chem{AMREX_D_DECL(-1, -1, -1)};

  // Times
  amrex::Vector<amrex::Real> m_t_old;
  amrex::Vector<amrex::Real> m_t_new;

  // Boundary Conditions
  amrex::BCRec m_phys_bc; // Store physical BC from input
#ifdef PELE_USE_EFIELD
  amrex::BCRec m_phiV_bc;    // Store separate phiV BC from input
  amrex::BCRec m_phiV_bcpol; // Store separate phiV BC polarity from input
#endif
  amrex::Vector<amrex::BCRec> m_bcrec_state;
  amrex::BCRec m_bcrec_divu;
  amrex::Vector<amrex::BCRec> m_bcrec_force;

  // Advection/Diffusion types
  amrex::Vector<int>
    m_AdvTypeState; // Tag Conservative/NonConservative state components
  amrex::Vector<int>
    m_DiffTypeState; // Tag Diffusive/NonDiffusive state components

  // Derived quantities
  PeleLMDeriveList derive_lst;

  // Derived quantities
  PeleLMDeriveList evaluate_lst;

  // Tagging/Refinement
  amrex::Vector<amrex::AMRErrorTag> errTags;
  std::string m_initial_grid_file;
  std::string m_regrid_file;
  amrex::Vector<amrex::BoxArray> m_initial_ba;
  amrex::Vector<amrex::BoxArray> m_regrid_ba;

  // Problem parameters
  ProbParm* prob_parm = nullptr;
  ProbParm* prob_parm_d = nullptr;

  // Patch parameters
  BPatch* boundary_patches = nullptr;
  BPatch* boundary_patches_d = nullptr;

#ifdef PELE_USE_SOOT
  // Soot parameters
  SootModel* soot_model = nullptr;
  bool do_soot_solve;
#endif

#ifdef PELE_USE_RADIATION
  // Radiation parameters
  std::unique_ptr<PeleRad::Radiation> rad_model;
  bool do_rad_solve;
#endif

  // Active controller
  int m_ctrl_active{0};
  int m_ctrl_useTemp{0};
  amrex::Real m_ctrl_tauControl{-1.0};
  amrex::Real m_ctrl_cfix{0.0};
  amrex::Real m_ctrl_coftOld{-1.0};
  amrex::Real m_ctrl_sest{0.0};
  amrex::Real m_ctrl_corr{1.0};
  amrex::Real m_ctrl_V_in{0.0};
  amrex::Real m_ctrl_V_in_old{0.0};
  amrex::Real m_ctrl_changeMax{1.0};
  amrex::Real m_ctrl_tBase{0.0};
  amrex::Real m_ctrl_dV{0.0};
  amrex::Real m_ctrl_scale{0.0};
  amrex::Real m_ctrl_zBase{0.0};
  amrex::Real m_ctrl_h{0.0};
  amrex::Real m_ctrl_velMax{-1.0};
  amrex::Real m_ctrl_temperature{-1.0};
  int m_ctrl_verbose{0};
  int m_ctrl_NavgPts{3};
  int m_ctrl_nfilled{-1};
  int m_ctrl_flameDir{AMREX_SPACEDIM - 1};
  int m_ctrl_pseudoGravity{0};
  int m_ctrl_method{3};
  std::string m_ctrl_AChistory{"AC_History.dat"};
  amrex::Vector<amrex::Real> m_ctrl_time_pts;
  amrex::Vector<amrex::Real> m_ctrl_velo_pts;
  amrex::Vector<amrex::Real> m_ctrl_cntl_pts;

  // PMF parameters
  pele::physics::PeleParams<pele::physics::PMF::PmfData::DataContainer>
    pmf_data;

  // Transport pointer
  static pele::physics::PeleParams<pele::physics::transport::TransParm<
    pele::physics::PhysicsType::eos_type,
    pele::physics::PhysicsType::transport_type>>
    trans_parms;

  // Reactor pointer
  std::string m_chem_integrator;
  std::unique_ptr<pele::physics::reactions::ReactorBase> m_reactor;

  // Turbulence injection
  pele::physics::turbinflow::TurbInflow turb_inflow;

  // Diagnostics
  amrex::Vector<std::unique_ptr<DiagBase>> m_diagnostics;
  amrex::Vector<std::string> m_diagVars;

  int m_verbose = 0;

  // IO options
  std::string m_restart_chkfile;
  std::string m_restart_pltfile;
  std::string pltfileSource{"LM"};
  std::string m_plot_file{"plt"};
  std::string m_check_file{"chk"};
  int m_derivePlotVarCount = 0;
  amrex::Vector<std::string> m_derivePlotVars;
  int m_plotChemDiag = 0;
  int m_plotHeatRelease = 0;
  int m_plotStateSpec = 0;
  int m_plot_int = 0;
  int m_plot_zeroEBcovered = 1;
  amrex::Real m_plot_per_approx = -1.;
  amrex::Real m_plot_per_exact = -1.;
  int m_check_int = 0;
  amrex::Real m_check_per = -1.;
  int m_message_int = 10;
  int m_evaluatePlotVarCount = 0;
  int m_plot_grad_p = 1;
  int m_ioDigits = 5;
  amrex::Vector<std::string> m_evaluatePlotVars;
  bool m_write_hdf5_pltfile = false;
  bool m_do_patch_flow_variables = false;

  //-----------------------------------------------------------------------------
  // ALGORITHM

  // Time & step
  amrex::Real m_cur_time = 0.0;
  amrex::Real m_dt = -1.0;
  amrex::Real m_cfl = 0.7;
  amrex::Real m_fixed_dt = -1.0;
  amrex::Real m_init_dt = -1.0;
  amrex::Real m_min_dt = 1.0e-12;
  amrex::Real m_dtshrink = 1.0;
  amrex::Real m_prev_dt = -1.0;
  amrex::Real m_dtChangeMax = 1.1;
  amrex::Real m_stop_time = -1.0;
  amrex::Real m_max_dt = 1e12;
  amrex::Real m_wall_start = 0.0;
  amrex::Real m_max_wall_time = -1.0;
  int m_nstep = -1;
  int m_max_step = -1;
  int m_nfiles = 1;

  int m_regrid_int = -1;
  int m_regrid_on_restart = 0;

  // Switch Evolve/Evaluate
  std::string m_run_mode = "normal";

  // Incompressible ?
  int m_incompressible = 0;
  amrex::Real m_rho = -1.0;
  amrex::Real m_mu = -1.0;

  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_gravity{
    AMREX_D_DECL(0.0, 0.0, 0.0)};
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> m_background_gp{
    AMREX_D_DECL(0.0, 0.0, 0.0)};

  int m_do_periodic_channel = 0;
  int m_periodic_channel_dir = -1;

  // DivU ? default ON when not incompressible
  int m_has_divu = 1;

  // DivU dt
  int m_divu_checkFlag{1};
  amrex::Real m_divu_dtFactor = 0.5;
  amrex::Real m_divu_rhoMin = 0.1;

  // Aux data size
  int m_nAux = 0;

  // R-z
  int m_rz_correction = 0;

  // cc->ec average
  int m_harm_avg_cen2edge = 0;

  // Diffusion
  int m_unity_Le = 0;
  int m_fixed_Pr = 0;
  int m_fixed_Le = 0;

  amrex::Real m_Lewis_inv = 1.0;
  amrex::Real m_Prandtl_inv = 1.0 / 0.7;
  amrex::Real m_Schmidt_inv = 1.0 / 0.7;

  int m_use_wbar = 1;
  int m_use_soret = 0;

  // LES Model
  bool m_do_les = false;
  bool m_plot_les = false;
  int m_les_verbose = 0;
  std::string m_les_model = "None";
  amrex::Real m_les_cs_smag = 0.18;
  amrex::Real m_les_cm_wale = 0.60;
  amrex::Real m_les_cs_sigma = 1.35;
  amrex::Vector<amrex::Real> m_turb_visc_time;

  // switch on/off different processes
  // TODO: actually implement that
  int m_do_react = 1;
  int m_do_diff = 1;
  int m_do_adv = 1;

  // Misc.
  int m_floor_species = 0;

  // Performances
  int m_checkMem{0};
  int m_doLoadBalance{0};
  int m_loadBalanceCost{LoadBalanceCost::Ncell};
  int m_loadBalanceMethod{LoadBalanceMethod::SFC};
  int m_loadBalanceCostChem{LoadBalanceCost::ChemFunctCallMax};
  int m_loadBalanceMethodChem{LoadBalanceMethod::Knapsack};
  amrex::Real m_loadBalanceKSfactor{1.2};
  amrex::Real m_loadBalanceEffRatioThreshold{1.1};
  amrex::Vector<std::unique_ptr<amrex::LayoutData<amrex::Real>>> m_costs;
  amrex::Vector<amrex::Real> m_loadBalanceEff;

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
  int m_Godunov_ppm_limiter = PPM::VanLeer;
  int m_Godunov_ForceInTrans = 0;

  // Temporals
  int m_do_temporals = 0;
  int m_do_patch_mfr = 0;
  amrex::Vector<std::unique_ptr<BPatch>> m_bPatches;
  int m_temp_int = 5;
  int m_do_extremas = 0;
  int m_do_massBalance = 0;
  int m_do_energyBalance = 0;
  int m_do_speciesBalance = 0;
  amrex::Real m_massOld;
  amrex::Real m_massNew;
  amrex::Real m_RhoHOld;
  amrex::Real m_RhoHNew;
  amrex::Array<amrex::Real, NUM_SPECIES> m_RhoYOld;
  amrex::Array<amrex::Real, NUM_SPECIES> m_RhoYNew;
  amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> m_domainMassFlux;
  amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> m_domainRhoHFlux;
  amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM * NUM_SPECIES> m_domainRhoYFlux;
  amrex::Array<amrex::Real, 2 * AMREX_SPACEDIM> m_domainUmacFlux;

  std::ofstream tmpStateFile;
  std::ofstream tmpExtremasFile;
  std::ofstream tmpMassFile;
  std::ofstream tmpSpecFile;
  std::ofstream tmppatchmfrFile;

  // Number of ghost cells
#ifdef AMREX_USE_EB
  int m_nGrowState = 4;
#else
  int m_nGrowState = 3;
#endif
  int m_nGrowAdv = 1;
  int m_nGrowdivu = 1;
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

  // EB bc/covered
  amrex::Vector<amrex::Real> coveredState_h;
  amrex::Gpu::DeviceVector<amrex::Real> coveredState_d;
  int m_isothermalEB = 0;

  // EB redistribution
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
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 1> m_Cweights;
  amrex::Real m_C0 = -1.0;
  amrex::Real m_C1 = -1.0;
  int m_Crevert = 0; // Flip the definition with 1 - ∑_k wgt*v_k
  //-----------------------------------------------------------------------------

  //-----------------------------------------------------------------------------
  // Linear Solvers
  std::unique_ptr<DiffusionOp> m_diffusion_op;
  std::unique_ptr<DiffusionOp> m_mcdiffusion_op;
  std::unique_ptr<DiffusionTensorOp> m_diffusionTensor_op;
  std::unique_ptr<Hydro::MacProjector> macproj;
  int m_macProjNeedReset{0};
  int m_macProjOldSize{0};

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
