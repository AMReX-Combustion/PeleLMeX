#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>
#include <PeleLMeX_DeriveFunc.H>
#include <PeleLMeX_BPatch.H>
#include "PelePhysics.H"
#include <AMReX_buildInfo.H>
#ifdef PELE_USE_EFIELD
#include "PeleLMeX_EOS_Extension.H"
#endif

#ifdef AMREX_USE_GPU
#include <AMReX_SUNMemory.H>
#endif

#ifdef PELE_USE_SOOT
#include "SootModel.H"
#endif
using namespace amrex;

static Box
the_same_box(const Box& b)
{
  return b;
}
#ifdef PELE_USE_EFIELD
static Box
grow_box_by_one(const Box& b)
{
  return amrex::grow(b, 1);
}
#endif
static Box
grow_box_by_two(const Box& b)
{
  return amrex::grow(b, 2);
}

void
PeleLM::Setup()
{
  BL_PROFILE("PeleLMeX::Setup()");

  m_wall_start = amrex::ParallelDescriptor::second();

  // Ensure grid is isotropic
  {
    auto const dx = geom[0].CellSizeArray();
    AMREX_ALWAYS_ASSERT(AMREX_D_TERM(
      , amrex::almostEqual(dx[0], dx[1], 10),
      &&amrex::almostEqual(dx[1], dx[2], 10)));
  }
  // Print build info to screen
  const char* githash1 = buildInfoGetGitHash(1);
  const char* githash2 = buildInfoGetGitHash(2);
  const char* githash3 = buildInfoGetGitHash(3);
  const char* githash4 = buildInfoGetGitHash(4);
  amrex::Print() << "\n ================= Build infos =================\n";
  amrex::Print() << " PeleLMeX    git hash: " << githash1 << "\n";
  amrex::Print() << " AMReX       git hash: " << githash2 << "\n";
  amrex::Print() << " PelePhysics git hash: " << githash3 << "\n";
  amrex::Print() << " AMReX-Hydro git hash: " << githash4 << "\n";
  amrex::Print() << " ===============================================\n";

#ifdef PELE_USE_SOOT
  soot_model = new SootModel{};
#endif

  // Read PeleLMeX parameters
  readParameters();

#ifdef AMREX_USE_EB
  makeEBGeometry();
#endif

  // Setup the state variables
  variablesSetup();

  // Derived variables
  derivedSetup();

  // Evaluate variables
  evaluateSetup();

  // Tagging setup
  taggingSetup();

#ifdef PELE_USE_SPRAY
  SpraySetup();
#endif
#ifdef PELE_USE_SOOT
  if (do_soot_solve) {
    soot_model->define();
  }
#endif
  // Diagnostics setup
  createDiagnostics();

  // Boundary Patch Setup
  if (m_do_patch_mfr != 0) {
    initBPatches(Geom(0));
  }

  // Initialize Level Hierarchy data
  resizeArray();

  // Initialize EOS and others
  if (m_incompressible == 0) {
    amrex::Print() << " Initialization of Transport ... \n";
    trans_parms.initialize();
    if ((m_les_verbose != 0) and m_do_les) { // Say what transport model we're
                                             // going to use
      amrex::Print() << "    Using LES in transport with Sc = "
                     << 1.0 / m_Schmidt_inv
                     << " and Pr = " << 1.0 / m_Prandtl_inv << std::endl;
    } else if (m_verbose != 0) {
      if (m_fixed_Le == 0 && m_fixed_Pr == 0) {
        if (m_use_soret == 0) {
          amrex::Print() << "    Using mixture-averaged transport" << std::endl;
        } else {
          amrex::Print()
            << "    Using mixture-averaged transport with Soret effects"
            << std::endl;
        }
      } else {
        if (m_fixed_Le != 0) {
          amrex::Print() << "    Using fixed Le = " << 1.0 / m_Lewis_inv
                         << std::endl;
        }
        if (m_fixed_Pr != 0) {
          amrex::Print() << "    Using fixed Pr = " << 1.0 / m_Prandtl_inv
                         << std::endl;
        }
      }
    }
    if (m_do_react != 0) {
      int reactor_type = 2;
      int ncells_chem = 1;
      amrex::Print() << " Initialization of chemical reactor ... \n";
      m_chem_integrator = "ReactorNull";
      ParmParse pp("peleLM");
      pp.query("chem_integrator", m_chem_integrator);
      m_reactor =
        pele::physics::reactions::ReactorBase::create(m_chem_integrator);
      m_reactor->init(reactor_type, ncells_chem);
      // For ReactorNull, we need to also skip instantaneous RR used in divU
      if (m_chem_integrator == "ReactorNull") {
        m_skipInstantRR = 1;
        m_plotChemDiag = 0;
        m_plotHeatRelease = 0;
        m_useTypValChem = 0;
        reactComponents.clear();
      }
      pp.query("plot_react", m_plot_react);
    }

#ifdef PELE_USE_EFIELD
    pele::physics::eos::charge_mass(zk.arr);
    for (int n = 0; n < NUM_SPECIES; n++) {
      zk[n] *= 1000.0; // CGS->MKS
    }
#endif
  }

  // Mixture fraction & Progress variable
  initMixtureFraction();
  initProgressVariable();

  // Initialize turbulence injection
  turb_inflow.init(Geom(0));

  // Initialize BCs
  setBoundaryConditions();

  // Problem parameters
  prob_parm = new ProbParm{};
  prob_parm_d = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm));

  // Problem parameters
  readProbParm();

  // Initialize ambient pressure
  // Will be overwritten on restart.
  m_pOld = prob_parm->P_mean;
  m_pNew = prob_parm->P_mean;

  // Copy problem parameters into device copy
  Gpu::copy(Gpu::hostToDevice, prob_parm, prob_parm + 1, prob_parm_d);

  // Initialize active control
  initActiveControl();
}

void
PeleLM::readParameters()
{
  BL_PROFILE("PeleLMeX::readParameters()");

  readIOParameters();

  ParmParse pp("peleLM");

  // -----------------------------------------
  // Misc
  // -----------------------------------------
  pp.query("run_mode", m_run_mode);
  pp.query("v", m_verbose);

  // -----------------------------------------
  // Boundary conditions
  // -----------------------------------------
  int isOpenDomain = 0;
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    int lo_bc = BoundaryCondition::BCInterior;
    int hi_bc = BoundaryCondition::BCInterior;
    parseUserKey(pp, "lo_bc", boundarycondition, lo_bc, idim);
    parseUserKey(pp, "hi_bc", boundarycondition, hi_bc, idim);
    m_phys_bc.setLo(idim, lo_bc);
    m_phys_bc.setHi(idim, hi_bc);
    if (
      lo_bc == BoundaryCondition::BCOutflow ||
      hi_bc == BoundaryCondition::BCOutflow) {
      isOpenDomain = 1;
    }
  }

  // Activate closed chamber if !isOpenDomain
  // enable overwriting
  m_closed_chamber = (isOpenDomain) != 0 ? 0 : 1;
  pp.query("closed_chamber", m_closed_chamber);
  if ((verbose != 0) && (m_closed_chamber != 0)) {
    Print() << " Simulation performed with the closed chamber algorithm \n";
  }

#ifdef PELE_USE_EFIELD
  ParmParse ppef("ef");

  // Get the phiV bc
  Vector<std::string> lo_bc_char(AMREX_SPACEDIM);
  Vector<std::string> hi_bc_char(AMREX_SPACEDIM);
  ppef.getarr("phiV_lo_bc", lo_bc_char, 0, AMREX_SPACEDIM);
  ppef.getarr("phiV_hi_bc", hi_bc_char, 0, AMREX_SPACEDIM);
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    if (lo_bc_char[idim] == "Interior") {
      m_phiV_bc.setLo(idim, 0);
    } else if (lo_bc_char[idim] == "Dirichlet") {
      m_phiV_bc.setLo(idim, 1);
    } else if (lo_bc_char[idim] == "Neumann") {
      m_phiV_bc.setLo(idim, 2);
    } else {
      amrex::Abort("Wrong PhiV bc. Should be : Interior, Dirichlet or Neumann");
    }
    if (hi_bc_char[idim] == "Interior") {
      m_phiV_bc.setHi(idim, 0);
    } else if (hi_bc_char[idim] == "Dirichlet") {
      m_phiV_bc.setHi(idim, 1);
    } else if (hi_bc_char[idim] == "Neumann") {
      m_phiV_bc.setHi(idim, 2);
    } else {
      amrex::Abort("Wrong PhiV bc. Should be : Interior, Dirichlet or Neumann");
    }
  }

  // Get the polarity of BCs
  ppef.getarr("phiV_polarity_lo", lo_bc_char, 0, AMREX_SPACEDIM);
  ppef.getarr("phiV_polarity_hi", hi_bc_char, 0, AMREX_SPACEDIM);
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    if (lo_bc_char[idim] == "Neutral") {
      m_phiV_bcpol.setLo(idim, 0);
    } else if (lo_bc_char[idim] == "Anode") { // Pos. elec = 1
      m_phiV_bcpol.setLo(idim, 1);
      AMREX_ASSERT(m_phiV_bc.lo(idim) == 1);
    } else if (lo_bc_char[idim] == "Cathode") { // Neg. elec = 2
      m_phiV_bcpol.setLo(idim, 2);
      AMREX_ASSERT(m_phiV_bc.lo(idim) == 1);
    } else {
      amrex::Abort(
        "Wrong PhiV polarity. Should be : Neutral, Anode or Cathode");
    }
    if (hi_bc_char[idim] == "Neutral") {
      m_phiV_bcpol.setHi(idim, 0);
    } else if (hi_bc_char[idim] == "Anode") { // Pos. elec = 1
      m_phiV_bcpol.setHi(idim, 1);
      AMREX_ASSERT(m_phiV_bc.hi(idim) == 1);
    } else if (hi_bc_char[idim] == "Cathode") { // Neg. elec = 2
      m_phiV_bcpol.setHi(idim, 2);
      AMREX_ASSERT(m_phiV_bc.hi(idim) == 1);
    } else {
      amrex::Abort(
        "Wrong PhiV polarity. Should be : Neutral, Anode or Cathode");
    }
  }
#endif

  // -----------------------------------------
  // Algorithm
  // -----------------------------------------

  // -----------------------------------------
  // incompressible vs. low Mach
  pp.query("use_divu", m_has_divu);
  pp.query("incompressible", m_incompressible);
  if (m_incompressible != 0) {
    m_has_divu = 0;
    m_do_react = 0;
    pp.query("rho", m_rho);
    pp.query("mu", m_mu);
    AMREX_ASSERT_WITH_MESSAGE(
      m_rho > 0.0, "peleLM.rho is needed when running incompressible");
    AMREX_ASSERT_WITH_MESSAGE(
      m_mu > 0.0, "peleLM.mu is needed when running incompressible");
  }
  Vector<Real> grav(AMREX_SPACEDIM, 0);
  pp.queryarr("gravity", grav, 0, AMREX_SPACEDIM);
  Vector<Real> gp0(AMREX_SPACEDIM, 0);
  pp.queryarr("gradP0", gp0, 0, AMREX_SPACEDIM);
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    m_background_gp[idim] = gp0[idim];
    m_gravity[idim] = grav[idim];
  }

  // Will automatically add pressure gradient for channel flow to maintain mass
  // flow rate of initial condition
  pp.query("do_periodic_channel", m_do_periodic_channel);
  if (m_do_periodic_channel != 0) {
    pp.get("periodic_channel_dir", m_periodic_channel_dir);
  }

  // Add a "spark", i.e. sphere of heat source to energy equation
  m_n_sparks = pp.countval("sparks");
  if (m_n_sparks > 0) {
    m_spark.resize(m_n_sparks);
    m_spark_time.resize(m_n_sparks);
    m_spark_duration.resize(m_n_sparks);
    m_spark_location.resize(m_n_sparks);
    m_spark_temp.resize(m_n_sparks);
    m_spark_radius.resize(m_n_sparks);
    pp.query("spark_verbose", m_spark_verbose);
    for (int n = 0; n < m_n_sparks; n++) {
      pp.get("sparks", m_spark[n], n);
      std::string spark_prefix = "peleLM." + m_spark[n];
      ParmParse pps(spark_prefix);
      pps.get("time", m_spark_time[n]);
      pps.get("duration", m_spark_duration[n]);
      m_spark_location[n].resize(AMREX_SPACEDIM);
      pps.getarr("location", m_spark_location[n], 0, AMREX_SPACEDIM);
      pps.get("temp", m_spark_temp[n]);
      pps.get("radius", m_spark_radius[n]);
    }
    if (m_spark_verbose > 0) {
      Print() << "Spark list:" << std::endl;
      for (int n = 0; n < m_n_sparks; n++) {
        Print() << "Spark " << n << " name: " << m_spark[n] << std::endl;
        Print() << "Spark " << n << " time: " << m_spark_time[n] << std::endl;
        Print() << "Spark " << n << " duration: " << m_spark_duration[n]
                << std::endl;
        Print() << "Spark " << n << " location: ";
        for (int d = 0; d < AMREX_SPACEDIM; d++) {
          Print() << m_spark_location[n][d] << " ";
        }
        Print() << std::endl;
        Print() << "Spark " << n << " temperature: " << m_spark_temp[n]
                << std::endl;
        Print() << "Spark " << n << " radius: " << m_spark_radius[n]
                << std::endl;
      }
    }
  }

  // -----------------------------------------
  // LES
  // -----------------------------------------
  pp.query("les_model", m_les_model);
  if (m_les_model == "None") {
    m_do_les = false;
  } else {
    if (m_les_model == "Smagorinsky") {
      pp.query("les_cs_smag", m_les_cs_smag);
    } else if (m_les_model == "WALE") {
      pp.query("les_cm_wale", m_les_cm_wale);
    } else if (m_les_model == "Sigma") {
      pp.query("les_cs_sigma", m_les_cs_sigma);
      AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == 3); // Sigma only available in 3D
    } else {
      amrex::Abort(
        "LES model must be None, Smagorinsky, WALE or Sigma. Invalid choie: " +
        m_les_model);
    }
    m_do_les = true;
    m_les_verbose = m_verbose;
    pp.query("plot_les", m_plot_les);
    pp.query("les_v", m_les_verbose);
    for (int lev = 0; lev <= max_level; ++lev) {
      m_turb_visc_time.push_back(-1.0E200);
    }
#ifdef PELE_USE_EFIELD
    amrex::Abort("LES implementation is not yet compatible with efield/ions");
#endif
  }

  // -----------------------------------------
  // diffusion
  ParmParse pptrans("transport");
  pptrans.query("use_soret", m_use_soret);
  pp.query("use_wbar", m_use_wbar);
  pp.query("unity_Le", m_unity_Le);
  pp.query("fixed_Le", m_fixed_Le);
  pp.query("fixed_Pr", m_fixed_Pr);
  if (m_unity_Le != 0) {
    m_fixed_Le = 1;
    amrex::Print() << "WARNING: unity_Le is deprecated and will be removed in "
                      "future version, use fixed_Le instead"
                   << std::endl;
  }
  if (m_do_les) { // For LES, Prandtl and Schmidt numbers are fixed
    m_fixed_Le = 1;
    m_fixed_Pr = 1;
    amrex::Real Schmidt = 0.7;
    pp.query("Schmidt", Schmidt);
    m_Schmidt_inv = 1.0 / Schmidt;
  }
  if (m_fixed_Le != 0 && !m_do_les) { // Only ask for Lewis number when not
                                      // LES, determined by Prandtl and
                                      // Schmidt outside of this
    amrex::Real Lewis = 1.0;
    pp.query("Lewis", Lewis);
    m_Lewis_inv = 1.0 / Lewis;
  }
  if (m_fixed_Pr != 0) {
    amrex::Real Prandtl = 0.7;
    pp.query("Prandtl", Prandtl);
    m_Prandtl_inv = 1.0 / Prandtl;
  }
  if (m_fixed_Le != 0 && m_fixed_Pr != 0 && !m_do_les) { // calculate Schmidt in
                                                         // case of no LES from
                                                         // Lewis and Prandtl
    m_Schmidt_inv = m_Lewis_inv * m_Prandtl_inv;
  }
  if (m_do_les) { // calculate Lewis in case of LES
    m_Lewis_inv = m_Prandtl_inv / m_Schmidt_inv;
  }

  if (
    (m_use_wbar != 0 || m_use_soret != 0) &&
    (m_fixed_Le != 0 || m_fixed_Pr != 0)) {
    m_use_wbar = 0;
    m_use_soret = 0;
    amrex::Print() << "WARNING: use_wbar and use_soret set to false because "
                      "fixed_Pr or fixed_Le is true"
                   << std::endl;
  }

  pp.query("deltaT_verbose", m_deltaT_verbose);
  pp.query("deltaT_iterMax", m_deltaTIterMax);
  pp.query("deltaT_tol", m_deltaT_norm_max);
  pp.query("deltaT_crashIfFailing", m_crashOnDeltaTFail);

  // -----------------------------------------
  // initialization
  pp.query("num_divu_iter", m_numDivuIter);
  pp.query("do_init_proj", m_do_init_proj);
  pp.query("num_init_iter", m_init_iter);
  pp.query("initDataPlt_patch_flow_variables", m_do_patch_flow_variables);
  pp.query("initDataPlt_reset_time", m_do_reset_time);

  // -----------------------------------------
  // advance
  // -----------------------------------------
  pp.query("sdc_iterMax", m_nSDCmax);
  pp.query("floor_species", m_floor_species);
  pp.query("dPdt_factor", m_dpdtFactor);
  pp.query("memory_checks", m_checkMem);
  pp.query("divu_dt_factor", m_divu_dtFactor);
  pp.query("divu_dt_rhoMin", m_divu_rhoMin);
  pp.query("divu_dt_method", m_divu_checkFlag);

  // -----------------------------------------
  // Reaction
  // -----------------------------------------
  pp.query("do_react", m_do_react);
  pp.query("use_typ_vals_chem", m_useTypValChem);
  pp.query("typical_values_reset_int", m_resetTypValInt);
  pp.query("typical_values_Ymin", m_typicalYvalMin);
  if (m_do_react != 0) {
    m_plotChemDiag = 0;
    m_plotHeatRelease = 1;
    pp.query("plot_chemDiagnostics", m_plotChemDiag);
    pp.query("plot_heatRelease", m_plotHeatRelease);
  }
  // Enable the chemistry BA to have smaller grid size
  int mgsc_size = pp.countval("max_grid_size_chem");
  if (mgsc_size > 0) {
    if (mgsc_size == 1) {
      int mgsc;
      pp.query("max_grid_size_chem", mgsc);
      AMREX_D_TERM(m_max_grid_size_chem[0] = mgsc;
                   , m_max_grid_size_chem[1] = mgsc;
                   , m_max_grid_size_chem[2] = mgsc);
    } else if (mgsc_size == AMREX_SPACEDIM) {
      Vector<int> mgsc;
      pp.getarr("max_grid_size_chem", mgsc, 0, AMREX_SPACEDIM);
      AMREX_D_TERM(m_max_grid_size_chem[0] = mgsc[0];
                   , m_max_grid_size_chem[1] = mgsc[1];
                   , m_max_grid_size_chem[2] = mgsc[2]);
    } else {
      Abort("peleLM.max_grid_size_chem should have 1 or AMREX_SPACEDIM values");
    }
  }

  // -----------------------------------------
  // Load Balancing
  // -----------------------------------------
  pp.query("do_load_balancing", m_doLoadBalance);
  parseUserKey(pp, "load_balancing_method", lbmethod, m_loadBalanceMethod);
  parseUserKey(pp, "load_balancing_cost_estimate", lbcost, m_loadBalanceCost);
  pp.query(
    "load_balancing_efficiency_threshold", m_loadBalanceEffRatioThreshold);
  parseUserKey(
    pp, "chem_load_balancing_method", lbmethod, m_loadBalanceMethodChem);
  parseUserKey(
    pp, "chem_load_balancing_cost_estimate", lbcost, m_loadBalanceCostChem);

  // Deactivate load balancing for serial runs
#ifdef AMREX_USE_MPI
  if (ParallelContext::NProcsSub() == 1) {
    m_doLoadBalance = 0;
  }
#else
  m_doLoadBalance = 0;
#endif

  // -----------------------------------------
  // Advection
  // -----------------------------------------
  pp.query("advection_scheme", m_advection_key);
  if (m_advection_key == "Godunov_PLM") {
    m_advection_type = "Godunov";
    m_Godunov_ppm = 0;
    ParmParse ppg("godunov");
    ppg.query("use_forceInTrans", m_Godunov_ForceInTrans);
  } else if (m_advection_key == "Godunov_PPM") {
    m_advection_type = "Godunov";
    m_Godunov_ppm = 1;
    m_Godunov_ppm_limiter = PPM::VanLeer;
    ParmParse ppg("godunov");
    ppg.query("use_forceInTrans", m_Godunov_ForceInTrans);
  } else if (m_advection_key == "Godunov_PPM_WENOZ") {
    m_advection_type = "Godunov";
    m_Godunov_ppm = 1;
    m_Godunov_ppm_limiter = PPM::WENOZ;
    ParmParse ppg("godunov");
    ppg.query("use_forceInTrans", m_Godunov_ForceInTrans);
  } else if (m_advection_key == "Godunov_PPM_NOLIM") {
    m_advection_type = "Godunov";
    m_Godunov_ppm = 1;
    m_Godunov_ppm_limiter = PPM::NoLimiter;
    ParmParse ppg("godunov");
    ppg.query("use_forceInTrans", m_Godunov_ForceInTrans);
  } else if (m_advection_key == "Godunov_BDS") {
    m_advection_type = "BDS";
    m_Godunov_ppm = 0;
  } else {
    Abort("Unknown 'advection_scheme'. Recognized options are: Godunov_PLM, "
          "Godunov_PPM or Godunov_BDS");
  }
  m_predict_advection_type =
    "Godunov"; // Only option at this point. This will disappear when
               // predict_velocity support BDS.

  // -----------------------------------------
  // Linear solvers tols
  // -----------------------------------------
  ParmParse ppnproj("nodal_proj");
  ppnproj.query("mg_max_coarsening_level", m_nodal_mg_max_coarsening_level);
  ppnproj.query("atol", m_nodal_mg_atol);
  ppnproj.query("rtol", m_nodal_mg_rtol);
  ppnproj.query("hypre_namespace", m_hypre_namespace_nodal);

  ParmParse ppmacproj("mac_proj");
  ppmacproj.query("mg_max_coarsening_level", m_mac_mg_max_coarsening_level);
  ppmacproj.query("atol", m_mac_mg_atol);
  ppmacproj.query("rtol", m_mac_mg_rtol);
  ppmacproj.query("hypre_namespace", m_hypre_namespace_mac);

  // -----------------------------------------
  // Temporals
  // -----------------------------------------
  pp.query("do_temporals", m_do_temporals);
  if (m_do_temporals != 0) {
    pp.query("temporal_int", m_temp_int);
    pp.query("do_extremas", m_do_extremas);
    pp.query("do_mass_balance", m_do_massBalance);
    pp.query("do_species_balance", m_do_speciesBalance);
    pp.query("do_patch_mfr", m_do_patch_mfr);
  }

  // -----------------------------------------
  // Time stepping control
  // -----------------------------------------
  ParmParse ppa("amr");
  ppa.query("max_wall_time", m_max_wall_time); // hours
  ppa.query("max_step", m_max_step);
  ppa.query("stop_time", m_stop_time);
  ppa.query("message_int", m_message_int);
  ppa.query("fixed_dt", m_fixed_dt);
  ppa.query("init_dt", m_init_dt);
  ppa.query("cfl", m_cfl);
  ppa.query("dt_shrink", m_dtshrink);
  ppa.query("dt_change_max", m_dtChangeMax);
  ppa.query("max_dt", m_max_dt);
  ppa.query("min_dt", m_min_dt);
  m_nfiles = std::max(1, std::min(ParallelDescriptor::NProcs(), 256));
  ppa.query("n_files", m_nfiles);

  if (max_level > 0 || (m_doLoadBalance != 0)) {
    ppa.query("regrid_int", m_regrid_int);
    ppa.query("regrid_on_restart", m_regrid_on_restart);
  }

#ifdef AMREX_USE_EB
  if (max_level > 0) {
    // Default EB refine type is Static
    pp.query("refine_EB_type", m_EB_refine_type);
    if (m_EB_refine_type != "Static" && m_EB_refine_type != "Adaptive") {
      Abort("refine_EB_type can only be 'Static' or 'Adaptive'");
    }
    // Default EB refinement level is max_level
    m_EB_refine_LevMax = max_level;
    pp.query("refine_EB_max_level", m_EB_refine_LevMax);
    pp.query("refine_EB_buffer", m_derefineEBBuffer);
    if (m_EB_refine_type == "Adaptive") {
      m_EB_refine_LevMin = 0;
      pp.query("refine_EB_min_level", m_EB_refine_LevMin);
      m_EB_refine_LevAdapt = m_EB_refine_LevMin;
    }
    if (m_EB_refine_LevMax < max_level) {
      m_signDistNeeded = 1;
    }
  }
  pp.query("isothermal_EB", m_isothermalEB);
#endif

  // -----------------------------------------
  // Evaluate mode variables
  // -----------------------------------------
  if (runMode() == "evaluate") {
    m_evaluatePlotVarCount = (pp.countval("evaluate_vars"));
    if (m_evaluatePlotVarCount != 0) {
      m_evaluatePlotVars.resize(m_evaluatePlotVarCount);
      for (int ivar = 0; ivar < m_evaluatePlotVarCount; ivar++) {
        pp.get("evaluate_vars", m_evaluatePlotVars[ivar], ivar);
      }
    }
  }

#ifdef PELE_USE_EFIELD
  // -----------------------------------------
  // EFIELD
  // -----------------------------------------
  ppef.query("JFNK_newtonTol", m_ef_newtonTol);
  ppef.query("JFNK_maxNewton", m_ef_maxNewtonIter);
  ppef.query("JFNK_lambda", m_ef_lambda_jfnk);
  ppef.query("JFNK_diffType", m_ef_diffT_jfnk);
  AMREX_ASSERT(m_ef_diffT_jfnk == 1 || m_ef_diffT_jfnk == 2);
  ppef.query("GMRES_rel_tol", m_ef_GMRES_reltol);
  ppef.query("GMRES_abs_tol", m_ef_GMRES_abstol);
  ppef.query("PC_approx", m_ef_PC_approx);
  ppef.query("PC_damping", m_ABecCecOmega);
  ppef.query("advection_scheme_order", m_nEAdvOrder);
  AMREX_ASSERT(m_nEAdvOrder == 1 || m_nEAdvOrder == 2);

  ppef.query("tabulated_Ke", m_electronKappaTab);
  ppef.query("fixed_Ke", m_fixedKappaE);

  ppef.query("restart_nonEF", m_restart_nonEF);
  ppef.query("restart_electroneutral", m_restart_electroneutral);
  ppef.query("restart_resetTime", m_restart_resetTime);
  ppef.query("plot_extras", m_do_extraEFdiags);

  // Getting the ions fluxes on the domain boundaries
  // Species balance data is needed, so override if not activated
  if (m_do_temporals) {
    ppef.query("do_ionsBalance", m_do_ionsBalance);
    if (m_do_ionsBalance) {
      m_do_speciesBalance = 1;
    }
  }
#endif

#ifdef PELE_USE_SPRAY
  SprayReadParameters();
#endif
#ifdef PELE_USE_SOOT
  do_soot_solve = true;
  pp.query("do_soot_solve", do_soot_solve);
  if ((m_verbose != 0) && do_soot_solve) {
    Print() << "Simulation performed with soot modeling \n";
  }
  soot_model->readSootParams();
#endif
#ifdef PELE_USE_RADIATION
  do_rad_solve = false;
  pp.query("do_rad_solve", do_rad_solve);
  if ((m_verbose != 0) && do_rad_solve) {
    Print() << "Simulation performed with radiation modeling \n";
  }
#endif
}

void
PeleLM::readIOParameters()
{
  BL_PROFILE_VAR("PeleLMeX::readIOParameters()", readIOParameters);

  ParmParse pp("amr");

  pp.query("check_file", m_check_file);
  pp.query("check_int", m_check_int);
  pp.query("check_per", m_check_per);
  pp.query("restart", m_restart_chkfile);
  pp.query("initDataPlt", m_restart_pltfile);
  pp.query("initDataPltSource", pltfileSource);
  pp.query("plot_file", m_plot_file);
  pp.query("plot_int", m_plot_int);
  if (pp.contains("plot_per")) {
    int do_exact = 0;
    pp.query("plot_per_exact", do_exact);
    if (do_exact != 0) {
      pp.query("plot_per", m_plot_per_exact);
    } else {
      pp.query("plot_per", m_plot_per_approx);
    }
  }
  pp.query("plot_grad_p", m_plot_grad_p);
  m_derivePlotVarCount = (pp.countval("derive_plot_vars"));
  if (m_derivePlotVarCount != 0) {
    m_derivePlotVars.resize(m_derivePlotVarCount);
    for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++) {
      pp.get("derive_plot_vars", m_derivePlotVars[ivar], ivar);
    }
  }
  pp.query("plot_zeroEBcovered", m_plot_zeroEBcovered);
  pp.query("plot_speciesState", m_plotStateSpec);
  m_initial_grid_file = "";
  m_regrid_file = "";
  pp.query("initial_grid_file", m_initial_grid_file);
  pp.query("regrid_file", m_regrid_file);
  pp.query("file_stepDigits", m_ioDigits);
  pp.query("use_hdf5_plt", m_write_hdf5_pltfile);
  pp.query("regrid_interp_method", m_regrid_interp_method);
  AMREX_ASSERT(m_regrid_interp_method == 0 || m_regrid_interp_method == 1);
}

void
PeleLM::variablesSetup()
{
  BL_PROFILE("PeleLMeX::variablesSetup()");

  //----------------------------------------------------------------
  // Variables ordering is defined through compiler macro in PeleLMeX_Index.H
  // Simply print on screen the state layout and append to the stateComponents
  // list
  Print() << "\n";
  Print() << PrettyLine;
  Print() << " State components \n";
  Print() << PrettyLine;

  Print() << " Velocity X: " << VELX;
  stateComponents.emplace_back(VELX, "x_velocity");
#if AMREX_SPACEDIM > 1
  Print() << ", Velocity Y: " << VELY;
  stateComponents.emplace_back(VELY, "y_velocity");
#if AMREX_SPACEDIM > 2
  Print() << ", Velocity Z: " << VELZ;
  stateComponents.emplace_back(VELZ, "z_velocity");
#endif
#endif
  Print() << " \n";

  if (m_incompressible == 0) {
    Print() << " Density: " << DENSITY << "\n";
    stateComponents.emplace_back(DENSITY, "density");
    Print() << " First species: " << FIRSTSPEC << "\n";
    Vector<std::string> names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
      names);
    for (int n = 0; n < NUM_SPECIES; n++) {
      stateComponents.emplace_back(FIRSTSPEC + n, "rho.Y(" + names[n] + ")");
      reactComponents.emplace_back(n, "I_R(" + names[n] + ")");
    }
    Print() << " Enthalpy: " << RHOH << "\n";
    stateComponents.emplace_back(RHOH, "rhoh");
    Print() << " Temperature: " << TEMP << "\n";
    stateComponents.emplace_back(TEMP, "temp");
    Print() << " thermo. pressure: " << RHORT << "\n";
    stateComponents.emplace_back(RHORT, "RhoRT");
#ifdef PELE_USE_EFIELD
    Print() << " nE: " << NE << "\n";
    stateComponents.emplace_back(NE, "nE");
    Print() << " PhiV: " << PHIV << "\n";
    stateComponents.emplace_back(PHIV, "PhiV");
#endif
#ifdef PELE_USE_SOOT
    for (int mom = 0; mom < NUMSOOTVAR; mom++) {
      std::string sootname = soot_model->sootVariableName(mom);
      Print() << " " << sootname << ": " << FIRSTSOOT + mom << "\n";
      stateComponents.emplace_back(FIRSTSOOT + mom, sootname);
    }
    setSootIndx();
#endif
  }

  if (m_nAux > 0) {
    Print() << " First passive scalar: " << FIRSTAUX << "\n";
    for (int n = 0; n < m_nAux; n++) {
      stateComponents.emplace_back(FIRSTAUX + n, "Aux_" + std::to_string(n));
    }
  }

  if (m_incompressible != 0) {
    Print() << " => Total number of state variables: " << AMREX_SPACEDIM
            << "\n";
  } else {
    Print() << " => Total number of state variables: " << NVAR << "\n";
  }
  Print() << PrettyLine;
  Print() << "\n";

  //----------------------------------------------------------------
  // Set advection/diffusion types
  if (m_incompressible != 0) {
    m_AdvTypeState.resize(AMREX_SPACEDIM);
    m_DiffTypeState.resize(AMREX_SPACEDIM);
  } else {
    m_AdvTypeState.resize(NVAR);
    m_DiffTypeState.resize(NVAR);
  }

  // Velocity - follow incflo
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    m_AdvTypeState[VELX + idim] = 0;  // NonConservative
    m_DiffTypeState[VELX + idim] = 1; // Diffusive
  }

  if (m_incompressible == 0) {
    m_AdvTypeState[DENSITY] = 1;
    m_DiffTypeState[DENSITY] = 0;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      m_AdvTypeState[FIRSTSPEC + n] = 1;
      m_DiffTypeState[FIRSTSPEC + n] = 1;
    }
    m_AdvTypeState[RHOH] = 1;
    m_DiffTypeState[RHOH] = 1;
    m_AdvTypeState[TEMP] = 0;
    m_DiffTypeState[TEMP] = 0;
    m_AdvTypeState[RHORT] = 0;
    m_DiffTypeState[RHORT] = 0;
#ifdef PELE_USE_EFIELD
    m_AdvTypeState[NE] = 0;
    m_DiffTypeState[NE] = 0;
    m_AdvTypeState[PHIV] = 0;
    m_DiffTypeState[PHIV] = 0;
#endif
#ifdef PELE_USE_SOOT
    for (int mom = 0; mom < NUMSOOTVAR; mom++) {
      m_AdvTypeState[FIRSTSOOT + mom] = 0;
      m_DiffTypeState[FIRSTSOOT + mom] = 0;
    }
#endif
  }

  //----------------------------------------------------------------
  // Typical values container
  if (m_incompressible != 0) {
    typical_values.resize(AMREX_SPACEDIM, -1.0);
  } else {
    typical_values.resize(NVAR, -1.0);
  }

  if (m_incompressible == 0) {
    // -----------------------------------------
    // Combustion
    // -----------------------------------------
    ParmParse pp("peleLM");
    std::string fuel_name;
    pp.query("fuel_name", fuel_name);
    fuel_name = "rho.Y(" + fuel_name + ")";
    if (isStateVariable(fuel_name)) {
      fuelID = stateVariableIndex(fuel_name);
      fuelID -= FIRSTSPEC;
    }
  }

  if (max_level > 0 && !m_initial_grid_file.empty()) {
    readGridFile(m_initial_grid_file, m_initial_ba);
    if (verbose > 0) {
      amrex::Print() << "Read initial_ba. Size is " << m_initial_ba.size()
                     << "\n";
    }
  }
  if (max_level > 0 && !m_regrid_file.empty() && m_regrid_int > 0) {
    readGridFile(m_regrid_file, m_regrid_ba);
    if (verbose > 0) {
      amrex::Print() << "Read regrid_ba. Size is " << m_regrid_ba.size()
                     << "\n";
    }
  }
}

void
PeleLM::readGridFile(
  std::string grid_file, amrex::Vector<amrex::BoxArray>& input_ba)
{
#define STRIP                \
  while (is.get() != '\n') { \
  }
  std::ifstream is(grid_file.c_str(), std::ios::in);
  if (!is.good()) {
    amrex::FileOpenFailed(grid_file);
  }
  int in_finest, ngrid;

  is >> in_finest;
  STRIP;
  input_ba.resize(in_finest);
  use_fixed_upto_level = in_finest;
  if (in_finest > max_level) {
    amrex::Error(
      "You have fewer levels in your inputs file then in your grids file!");
  }

  for (int lev = 1; lev <= in_finest; lev++) {
    BoxList bl;
    is >> ngrid;
    STRIP;
    for (int i = 0; i < ngrid; i++) {
      Box bx;
      is >> bx;
      STRIP;
      bx.refine(ref_ratio[lev - 1]);
      bl.push_back(bx);
    }
    input_ba[lev - 1].define(bl);
  }
  is.close();

#undef STRIP
}

void
PeleLM::derivedSetup()
{
  BL_PROFILE("PeleLMeX::derivedSetup()");

  if (m_incompressible == 0) {

    // Get species names
    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
      spec_names);

    // Set species mass fractions
    Vector<std::string> var_names_massfrac(NUM_SPECIES);
    for (int n = 0; n < NUM_SPECIES; n++) {
      var_names_massfrac[n] = "Y(" + spec_names[n] + ")";
    }
    derive_lst.add(
      "mass_fractions", IndexType::TheCellType(), NUM_SPECIES,
      var_names_massfrac, pelelmex_dermassfrac, the_same_box);

    for (int n = 0; n < NUM_SPECIES; n++) {
      var_names_massfrac[n] = "X(" + spec_names[n] + ")";
    }
    derive_lst.add(
      "mole_fractions", IndexType::TheCellType(), NUM_SPECIES,
      var_names_massfrac, pelelmex_dermolefrac, the_same_box);

    // Species diffusion coefficients
    for (int n = 0; n < NUM_SPECIES; n++) {
      var_names_massfrac[n] = "D_" + spec_names[n];
    }
    if (m_use_soret != 0) {
      var_names_massfrac.resize(2 * NUM_SPECIES);
      for (int n = 0; n < NUM_SPECIES; n++) {
        var_names_massfrac[n + NUM_SPECIES] = "theta_" + spec_names[n];
      }
      derive_lst.add(
        "diffcoeff", IndexType::TheCellType(), 2 * NUM_SPECIES,
        var_names_massfrac, pelelmex_derdiffc, the_same_box);
    } else {
      derive_lst.add(
        "diffcoeff", IndexType::TheCellType(), NUM_SPECIES, var_names_massfrac,
        pelelmex_derdiffc, the_same_box);
    }

    // Rho - sum rhoYs
    derive_lst.add(
      "rhominsumrhoY", IndexType::TheCellType(), 1, pelelmex_derrhomrhoy,
      the_same_box);

    // Heat Release
    derive_lst.add(
      "HeatRelease", IndexType::TheCellType(), 1, pelelmex_derheatrelease,
      the_same_box);

    // Thermal diffusivity
    derive_lst.add(
      "lambda", IndexType::TheCellType(), 1, pelelmex_derlambda, the_same_box);

    // Mixture fraction
    derive_lst.add(
      "mixture_fraction", IndexType::TheCellType(), 1, pelelmex_dermixfrac,
      the_same_box);

    // Progress variable
    derive_lst.add(
      "progress_variable", IndexType::TheCellType(), 1, pelelmex_derprogvar,
      the_same_box);
  }

  // Distribution Map
  derive_lst.add(
    "DistributionMap", IndexType::TheCellType(), 1, pelelmex_derdmap,
    the_same_box);

  // Cell average pressure
  derive_lst.add(
    "avg_pressure", IndexType::TheCellType(), 1, pelelmex_deravgpress,
    the_same_box);

  // Viscosity
  derive_lst.add(
    "viscosity", IndexType::TheCellType(), 1, pelelmex_dervisc, the_same_box);

  // Velocity magnitude
  derive_lst.add(
    "mag_vel", IndexType::TheCellType(), 1, pelelmex_dermgvel, the_same_box);

  // Vorticity magnitude
  derive_lst.add(
    "mag_vort", IndexType::TheCellType(), 1, pelelmex_dermgvort,
    grow_box_by_two);

  // Spatial coordinates
  {
    Vector<std::string> var_names({AMREX_D_DECL("x", "y", "z")});
    derive_lst.add(
      "coordinates", IndexType::TheCellType(), AMREX_SPACEDIM, var_names,
      pelelmex_dercoord, the_same_box);
  }

  // Vorticity components
  {
    const int vort_ncomp = 2 * AMREX_SPACEDIM - 3;
#if (AMREX_SPACEDIM == 2)
    Vector<std::string> var_names({"VortZ"});
#elif (AMREX_SPACEDIM == 3)
    Vector<std::string> var_names({"VortX", "VortY", "VortZ"});
#endif
    derive_lst.add(
      "vorticity", IndexType::TheCellType(), vort_ncomp, var_names,
      pelelmex_dervort, grow_box_by_two);
  }

  // UserDefined derived
  {
    Vector<std::string> var_names = pelelmex_setuserderives();
    derive_lst.add(
      "derUserDefined", IndexType::TheCellType(),
      static_cast<int>(var_names.size()), var_names, pelelmex_deruserdef,
      the_same_box);
  }

#if (AMREX_SPACEDIM == 3)
  // Q-criterion
  derive_lst.add(
    "Qcrit", IndexType::TheCellType(), 1, pelelmex_derQcrit, grow_box_by_two);
#endif

  // Kinetic energy
  derive_lst.add(
    "kinetic_energy", IndexType::TheCellType(), 1, pelelmex_derkineticenergy,
    the_same_box);

  // Enstrophy
  derive_lst.add(
    "enstrophy", IndexType::TheCellType(), 1, pelelmex_derenstrophy,
    grow_box_by_two);

#ifdef PELE_USE_EFIELD
  // Charge distribution
  derive_lst.add(
    "chargedistrib", IndexType::TheCellType(), 1, pelelmex_derchargedist,
    the_same_box);

  // Electric field
  derive_lst.add(
    "efieldx", IndexType::TheCellType(), 1, pelelmex_derefx, grow_box_by_one);
#if (AMREX_SPACEDIM > 1)
  derive_lst.add(
    "efieldy", IndexType::TheCellType(), 1, pelelmex_derefy, grow_box_by_one);
#if (AMREX_SPACEDIM > 2)
  derive_lst.add(
    "efieldz", IndexType::TheCellType(), 1, pelelmex_derefz, grow_box_by_one);
#endif
#endif
  // Lorentz forces
  derive_lst.add(
    "LorentzFx", IndexType::TheCellType(), 1, pelelmex_derLorentzx,
    grow_box_by_one);
#if (AMREX_SPACEDIM > 1)
  derive_lst.add(
    "LorentzFy", IndexType::TheCellType(), 1, pelelmex_derLorentzy,
    grow_box_by_one);
#if (AMREX_SPACEDIM > 2)
  derive_lst.add(
    "LorentzFz", IndexType::TheCellType(), 1, pelelmex_derLorentzz,
    grow_box_by_one);
#endif
#endif
#endif
  auto it = m_derivePlotVars.begin();
  while (it != m_derivePlotVars.end()) {
    if (!derive_lst.canDerive(*it)) {
      it = m_derivePlotVars.erase(it);
    } else {
      it++;
    }
  }
  m_derivePlotVarCount = static_cast<int>(m_derivePlotVars.size());
}

void
PeleLM::evaluateSetup()
{
  BL_PROFILE("PeleLMeX::evaluateSetup()");

  // Get species names
  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    spec_names);

  // divU
  evaluate_lst.add("divU", IndexType::TheCellType(), 1, the_same_box);

  // projected velocity field
  {
    Vector<std::string> var_names = {
      AMREX_D_DECL("x_velProj", "y_velProj", "z_velProj")};
    evaluate_lst.add(
      "velProj", IndexType::TheCellType(), AMREX_SPACEDIM, var_names,
      the_same_box);
  }

  // velocity force
  {
    Vector<std::string> var_names = {
      AMREX_D_DECL("x_velForce", "y_velForce", "z_velForce")};
    evaluate_lst.add(
      "velForce", IndexType::TheCellType(), AMREX_SPACEDIM, var_names,
      the_same_box);
  }

  // divTau
  {
    Vector<std::string> var_names = {
      AMREX_D_DECL("x_divTau", "y_divTau", "z_divTau")};
    evaluate_lst.add(
      "divTau", IndexType::TheCellType(), AMREX_SPACEDIM, var_names,
      the_same_box);
  }

  // scalar diffusion term
  {
    Vector<std::string> var_names(NUM_SPECIES + 2);
    for (int n = 0; n < NUM_SPECIES; n++) {
      var_names[n] = "D(" + spec_names[n] + ")";
    }
    var_names[NUM_SPECIES] = "D(RhoH)";
    var_names[NUM_SPECIES + 1] = "D(Temp)";
    evaluate_lst.add(
      "diffTerm", IndexType::TheCellType(), NUM_SPECIES + 2, var_names,
      the_same_box);
  }

  // advection terms
  {
    Vector<std::string> var_names(
      NVAR - 2); // Skip temperature and RhoRT, unused
    AMREX_D_TERM(var_names[VELX] = "A(VELX)";, var_names[VELY] = "A(VELY)";
                 , var_names[VELZ] = "A(VELZ)");
    var_names[DENSITY] = "A(Rho)";
    for (int n = 0; n < NUM_SPECIES; n++) {
      var_names[FIRSTSPEC + n] = "A(" + spec_names[n] + ")";
    }
    var_names[RHOH] = "A(RhoH)";
    evaluate_lst.add(
      "advTerm", IndexType::TheCellType(), NVAR - 2, var_names, the_same_box);
  }

  // Chemical state and external chem. forcing (used in ReactEval)
  {
    Vector<std::string> var_names(2 * (NUM_SPECIES + 1) + 1);
    for (int n = 0; n < NUM_SPECIES; n++) {
      var_names[n] = "rhoY(" + spec_names[n] + ")";
    }
    var_names[NUM_SPECIES] = "rhoH";
    var_names[NUM_SPECIES + 1] = "Temp";
    for (int n = 0; n < NUM_SPECIES; n++) {
      var_names[NUM_SPECIES + 2 + n] = "F_rhoY(" + spec_names[n] + ")";
    }
    var_names[2 * NUM_SPECIES + 2] = "F_rhoH";
    evaluate_lst.add(
      "chemTest", IndexType::TheCellType(), 2 * (NUM_SPECIES + 1) + 1,
      var_names, the_same_box);
  }

  // instantaneous reaction rate
  {
    Vector<std::string> var_names(NUM_SPECIES);
    for (int n = 0; n < NUM_SPECIES; n++) {
      var_names[n] = "I_R(" + spec_names[n] + ")";
    }
    evaluate_lst.add(
      "instRR", IndexType::TheCellType(), NUM_SPECIES, var_names, the_same_box);
  }

  // cell-centered transport coefficients
  {
    Vector<std::string> var_names(NUM_SPECIES + 2);
    for (int n = 0; n < NUM_SPECIES; n++) {
      var_names[n] = "rhoD(" + spec_names[n] + ")";
    }
    var_names[NUM_SPECIES] = "Lamdba";
    var_names[NUM_SPECIES + 1] = "Mu";
    evaluate_lst.add(
      "transportCC", IndexType::TheCellType(), NUM_SPECIES + 2, var_names,
      the_same_box);
  }
}

void
PeleLM::taggingSetup()
{
  BL_PROFILE("PeleLMeX::taggingSetup()");

  std::string amr_prefix = "amr";
  ParmParse ppamr(amr_prefix);

  Vector<std::string> refinement_indicators;
  ppamr.queryarr(
    "refinement_indicators", refinement_indicators, 0,
    ppamr.countval("refinement_indicators"));
  for (const auto& refinement_indicator : refinement_indicators) {
    std::string ref_prefix = amr_prefix + "." + refinement_indicator;
    ParmParse ppr(ref_prefix);

    // Tag a given box
    RealBox realbox;
    if (ppr.countval("in_box_lo") != 0) {
      Vector<Real> box_lo(AMREX_SPACEDIM);
      Vector<Real> box_hi(AMREX_SPACEDIM);
      ppr.getarr("in_box_lo", box_lo, 0, static_cast<int>(box_lo.size()));
      ppr.getarr("in_box_hi", box_hi, 0, static_cast<int>(box_hi.size()));
      realbox = RealBox(box_lo.data(), box_hi.data());
    }

    AMRErrorTagInfo info;

    if (realbox.ok()) {
      info.SetRealBox(realbox);
    }

    if (ppr.countval("start_time") > 0) {
      Real min_time;
      ppr.get("start_time", min_time);
      info.SetMinTime(min_time);
    }

    if (ppr.countval("end_time") > 0) {
      Real max_time;
      ppr.get("end_time", max_time);
      info.SetMaxTime(max_time);
    }

    if (ppr.countval("max_level") > 0) {
      int tag_max_level;
      ppr.get("max_level", tag_max_level);
      info.SetMaxLevel(tag_max_level);
    }

    bool itexists = false;
    if (ppr.countval("value_greater") != 0) {
      Real value;
      ppr.get("value_greater", value);
      std::string field;
      ppr.get("field_name", field);
      errTags.push_back(AMRErrorTag(value, AMRErrorTag::GREATER, field, info));
      itexists = derive_lst.canDerive(field) || isStateVariable(field) ||
                 isReactVariable(field);
    } else if (ppr.countval("value_less") != 0) {
      Real value;
      ppr.get("value_less", value);
      std::string field;
      ppr.get("field_name", field);
      errTags.push_back(AMRErrorTag(value, AMRErrorTag::LESS, field, info));
      itexists = derive_lst.canDerive(field) || isStateVariable(field) ||
                 isReactVariable(field);
    } else if (ppr.countval("vorticity_greater") != 0) {
      Real value;
      ppr.get("vorticity_greater", value);
      const std::string field = "mag_vort";
      errTags.push_back(AMRErrorTag(value, AMRErrorTag::VORT, field, info));
      itexists = derive_lst.canDerive(field) || isStateVariable(field) ||
                 isReactVariable(field);
    } else if (ppr.countval("adjacent_difference_greater") != 0) {
      Real value;
      ppr.get("adjacent_difference_greater", value);
      std::string field;
      ppr.get("field_name", field);
      errTags.push_back(AMRErrorTag(value, AMRErrorTag::GRAD, field, info));
      itexists = derive_lst.canDerive(field) || isStateVariable(field) ||
                 isReactVariable(field);
    } else if (realbox.ok()) {
      errTags.push_back(AMRErrorTag(info));
      itexists = true;
    } else {
      Abort(std::string(
              "Unrecognized refinement indicator for " + refinement_indicator)
              .c_str());
    }

    if (!itexists) {
      amrex::Error(
        "PeleLM::taggingSetup(): unknown variable field for criteria " +
        refinement_indicator);
    }
  }
}

void
PeleLM::resizeArray()
{

  if (m_verbose != 0) {
    Print() << " Initializing data for " << max_level + 1 << " levels \n";
  }

  // State data
  m_leveldata_old.resize(max_level + 1);
  m_leveldata_new.resize(max_level + 1);
  m_leveldatareact.resize(max_level + 1);
  m_halfTimeDensity.resize(max_level + 1);
  if (max_level > 0) {
    m_coveredMask.resize(max_level);
  }
  m_baChem.resize(max_level + 1);
  m_dmapChem.resize(max_level + 1);
  m_baChemFlag.resize(max_level + 1);

#ifdef PELE_USE_EFIELD
  m_leveldatanlsolve.resize(max_level + 1);
  m_ionsFluxes.resize(max_level + 1);
#endif

  // External sources
  m_extSource.resize(max_level + 1);

  // Factory
  m_factory.resize(max_level + 1);

  // Time
  m_t_old.resize(max_level + 1);
  m_t_new.resize(max_level + 1);

#ifdef PELE_USE_SPRAY
  m_spraystate.resize(max_level + 1);
  m_spraysource.resize(max_level + 1);
#endif

  // Load balancing
  m_costs.resize(max_level + 1);
  m_loadBalanceEff.resize(max_level + 1);
}
