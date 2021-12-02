#include <PeleLM.H>
#include <AMReX_ParmParse.H>
#include <PeleLMDeriveFunc.H>
#include "PelePhysics.H"
#include <AMReX_buildInfo.H>
#ifdef PELE_USE_EFIELD
#include "EOS_Extension.H"
#endif

#ifdef AMREX_USE_GPU
#include <AMReX_SUNMemory.H>
#endif

using namespace amrex;

static Box the_same_box (const Box& b)    { return b;                }
static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }
static Box grow_box_by_two (const Box& b) { return amrex::grow(b,2); }

void PeleLM::Setup() {
   BL_PROFILE_VAR("PeleLM::Setup()", Setup);

#ifdef AMREX_USE_GPU
   sundials::MemoryHelper::Initialize();
#endif

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

   // Read PeleLM parameters
   readParameters();

#ifdef AMREX_USE_EB
   makeEBGeometry();
#endif

   // Setup the state variables
   variablesSetup();

   // Derived variables
   derivedSetup();

   // Tagging setup
   taggingSetup();

   // Initialize Level Hierarchy data
   resizeArray();

   // Initialize EOS and others
   if (!m_incompressible) {
      amrex::Print() << " Initialization of Transport ... \n";
      trans_parms.allocate();
      if (m_do_react) {
         int reactor_type = 2;
         int ncells_chem = 1;
         amrex::Print() << " Initialization of chemical reactor ... \n";
         m_chem_integrator = "ReactorNull";
         ParmParse pp("peleLM");
         pp.query("chem_integrator",m_chem_integrator);
         m_reactor = pele::physics::reactions::ReactorBase::create(m_chem_integrator);
         m_reactor->init(reactor_type, ncells_chem);
         // For ReactorNull, we need to also skip instantaneous RR used in divU
         if (m_chem_integrator == "ReactorNull") {
            m_skipInstantRR = 1;
            m_plotChemDiag = 0;
            m_plotHeatRelease = 0;
         }
      }

#ifdef PELE_USE_EFIELD
      pele::physics::eos::charge_mass(zk.arr);
      for (int n = 0; n < NUM_SPECIES; n++) {
         zk[n] *= 1000.0;    // CGS->MKS
      }
#endif
   }

   // Initiliaze BCs
   setBoundaryConditions();

   // Problem parameters
   prob_parm = new ProbParm{};
   prob_parm_d = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm));

   // Problem parameters
   readProbParm();

   // Initialize ambient pressure
   // Will be overwriten on restart.
   m_pOld = prob_parm->P_mean;
   m_pNew = prob_parm->P_mean;

   // Copy problem parameters into device copy
   Gpu::copy(Gpu::hostToDevice, prob_parm, prob_parm+1,prob_parm_d); 
}

void PeleLM::readParameters() {
   BL_PROFILE_VAR("PeleLM::readParameters()", readParameters);

   readIOParameters();

   ParmParse pp("peleLM");

   // -----------------------------------------
   // Misc
   // -----------------------------------------
   pp.query("v", m_verbose);

   // -----------------------------------------
   // Boundary conditions
   // -----------------------------------------
   int isOpenDomain = 0;

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
         isOpenDomain = 1;
      } else if (lo_bc_char[dir] == "Outflow") {
         lo_bc[dir] = 2;
         isOpenDomain = 1;
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
         isOpenDomain = 1;
      } else if (hi_bc_char[dir] == "Outflow") {
         hi_bc[dir] = 2;
         isOpenDomain = 1;
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

   // Store BCs in m_phys_bc.
   for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      m_phys_bc.setLo(idim,lo_bc[idim]);
      m_phys_bc.setHi(idim,hi_bc[idim]);
   }

   // Activate closed chamber if !isOpenDomain
   // enable overwriting
   m_closed_chamber = (isOpenDomain) ? 0 : 1;
   pp.query("closed_chamber", m_closed_chamber);
   if (verbose && m_closed_chamber) {
      Print() << " Simulation performed with the closed chamber algorithm \n";
   }

#ifdef PELE_USE_EFIELD
   ParmParse ppef("ef");

   // Get the phiV bc
   ppef.getarr("phiV_lo_bc",lo_bc_char,0,AMREX_SPACEDIM);
   ppef.getarr("phiV_hi_bc",hi_bc_char,0,AMREX_SPACEDIM);
   for (int idim = 0; idim < AMREX_SPACEDIM; idim++) 
   {
      if (lo_bc_char[idim] == "Interior"){
         m_phiV_bc.setLo(idim,0);
      } else if (lo_bc_char[idim] == "Dirichlet") {
         m_phiV_bc.setLo(idim,1);
      } else if (lo_bc_char[idim] == "Neumann") {
         m_phiV_bc.setLo(idim,2);
      } else {
         amrex::Abort("Wrong PhiV bc. Should be : Interior, Dirichlet or Neumann");
      }    
      if (hi_bc_char[idim] == "Interior"){
         m_phiV_bc.setHi(idim,0);
      } else if (hi_bc_char[idim] == "Dirichlet") {
         m_phiV_bc.setHi(idim,1);
      } else if (hi_bc_char[idim] == "Neumann") {
         m_phiV_bc.setHi(idim,2);
      } else {
         amrex::Abort("Wrong PhiV bc. Should be : Interior, Dirichlet or Neumann");
      }
   }

   // Get the polarity of BCs
   ppef.getarr("phiV_polarity_lo",lo_bc_char,0,AMREX_SPACEDIM);
   ppef.getarr("phiV_polarity_hi",hi_bc_char,0,AMREX_SPACEDIM);
   for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      if (lo_bc_char[idim] == "Neutral"){
         m_phiV_bcpol.setLo(idim,0);
      } else if (lo_bc_char[idim] == "Anode") {     // Pos. elec = 1
         m_phiV_bcpol.setLo(idim,1);
         AMREX_ASSERT(m_phiV_bc.lo(idim)==1);
      } else if (lo_bc_char[idim]  == "Cathode") {   // Neg. elec = 2
         m_phiV_bcpol.setLo(idim,2);
         AMREX_ASSERT(m_phiV_bc.lo(idim)==1);
      } else {
         amrex::Abort("Wrong PhiV polarity. Should be : Neutral, Anode or Cathode");
      }
      if (hi_bc_char[idim] == "Neutral"){
         m_phiV_bcpol.setHi(idim,0);
      } else if (hi_bc_char[idim] == "Anode") {     // Pos. elec = 1
         m_phiV_bcpol.setHi(idim,1);
         AMREX_ASSERT(m_phiV_bc.hi(idim)==1);
      } else if (hi_bc_char[idim] == "Cathode") {   // Neg. elec = 2
         m_phiV_bcpol.setHi(idim,2);
         AMREX_ASSERT(m_phiV_bc.hi(idim)==1);
      } else {
         amrex::Abort("Wrong PhiV polarity. Should be : Neutral, Anode or Cathode");
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
   if (m_incompressible) {
      m_has_divu = 0;
      m_do_react = 0;
   }
   pp.query("rho", m_rho);
   pp.query("mu", m_mu);
   pp.query("gravity", m_gravity);

   // -----------------------------------------
   // diffusion
   pp.query("use_wbar",m_use_wbar);
   pp.query("deltaT_verbose",m_deltaT_verbose);
   pp.query("deltaT_iterMax",m_deltaTIterMax);
   pp.query("deltaT_tol",m_deltaT_norm_max);

   // -----------------------------------------
   // initialization
   pp.query("num_divu_iter",m_numDivuIter);
   pp.query("do_init_proj",m_do_init_proj);
   pp.query("num_init_iter",m_init_iter);

   // -----------------------------------------
   // advance
   // -----------------------------------------
   pp.query("sdc_iterMax",m_nSDCmax);
   pp.query("floor_species",m_floor_species);

   // -----------------------------------------
   // Reaction
   // -----------------------------------------
   pp.query("do_react",m_do_react);
   pp.query("use_typ_vals_chem",m_useTypValChem);
   pp.query("typical_values_reset_int",m_resetTypValInt);
   if (m_do_react) {
      m_plotChemDiag = 0;
      m_plotHeatRelease = 1;
      pp.query("plot_chemDiagnostics",m_plotChemDiag);
      pp.query("plot_heatRelease",m_plotHeatRelease);
   }

   // -----------------------------------------
   // Advection
   // -----------------------------------------
   ParmParse ppg("godunov");
   ppg.query("use_ppm",m_Godunov_ppm);
   ppg.query("use_forceInTrans", m_Godunov_ForceInTrans);

   // -----------------------------------------
   // Linear solvers tols
   // -----------------------------------------
   ParmParse ppnproj("nodal_proj");
   ppnproj.query("mg_max_coarsening_level",m_nodal_mg_max_coarsening_level);
   ppnproj.query("atol",m_nodal_mg_atol);
   ppnproj.query("rtol",m_nodal_mg_rtol);
   ppnproj.query("hypre_namespace",m_hypre_namespace_nodal);

   ParmParse ppmacproj("mac_proj");
   ppmacproj.query("mg_max_coarsening_level",m_mac_mg_max_coarsening_level);
   ppmacproj.query("atol",m_mac_mg_atol);
   ppmacproj.query("rtol",m_mac_mg_rtol);
   ppmacproj.query("hypre_namespace",m_hypre_namespace_mac);

   // -----------------------------------------
   // Temporals
   // -----------------------------------------
   pp.query("do_temporals",m_do_temporals);
   pp.query("temporal_int",m_temp_int);
   pp.query("mass_balance",m_do_massBalance);

   // -----------------------------------------
   // Time stepping control
   // -----------------------------------------
   ParmParse ppa("amr");
   ppa.query("max_step", m_max_step);
   ppa.query("stop_time", m_stop_time);
   ppa.query("fixed_dt", m_fixed_dt);
   ppa.query("cfl", m_cfl);
   ppa.query("dt_shrink", m_dtshrink);
   ppa.query("dt_change_max", m_dtChangeMax);
   ppa.query("max_dt", m_max_dt);

   if ( max_level > 0 ) {
      ppa.query("regrid_int", m_regrid_int);
   }

#ifdef PELE_USE_EFIELD
   // -----------------------------------------
   // EFIELD
   // -----------------------------------------
   ppef.query("JFNK_newtonTol",m_ef_newtonTol);
   ppef.query("JFNK_maxNewton",m_ef_maxNewtonIter);
   ppef.query("JFNK_lambda",m_ef_lambda_jfnk);
   ppef.query("GMRES_rel_tol",m_ef_GMRES_reltol);
   ppef.query("PC_approx",m_ef_PC_approx);

   ppef.query("restart_nonEF",m_restart_nonEF);
   ppef.query("restart_electroneutral",m_restart_electroneutral);
   ppef.query("restart_resetTime",m_restart_resetTime);
#endif

}

void PeleLM::readIOParameters() {
   BL_PROFILE_VAR("PeleLM::readIOParameters()", readIOParameters);

   ParmParse pp("amr");

   pp.query("check_file", m_check_file);
   pp.query("check_int" , m_check_int);
   pp.query("restart" , m_restart_file);
   pp.query("plot_file", m_plot_file);
   pp.query("plot_int" , m_plot_int);
   m_derivePlotVarCount = (pp.countval("derive_plot_vars"));
   if (m_derivePlotVarCount != 0) {
      m_derivePlotVars.resize(m_derivePlotVarCount);
      for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++) {
         pp.get("derive_plot_vars", m_derivePlotVars[ivar],ivar);
      }
   }
   pp.query("plot_speciesState" , m_plotStateSpec);

}

void PeleLM::variablesSetup() {
   BL_PROFILE_VAR("PeleLM::variablesSetup()", variablesSetup);

   std::string PrettyLine = std::string(78, '=') + "\n";

   //----------------------------------------------------------------
   // Variables ordering is defined through compiler macro in PeleLM_Index.H 
   // Simply print on screen the state layout and append to the stateComponents list
   Print() << "\n";
   Print() << PrettyLine;
   Print() << " State components \n";
   Print() << PrettyLine;

   Print() << " Velocity X: " << VELX;
   stateComponents.emplace_back(VELX,"x_velocity");
#if AMREX_SPACEDIM > 1
   Print() << ", Velocity Y: " << VELY;
   stateComponents.emplace_back(VELY,"y_velocity");
#if AMREX_SPACEDIM > 2
   Print() << ", Velocity Z: " << VELZ;
   stateComponents.emplace_back(VELZ,"z_velocity");
#endif
#endif
   Print() << " \n";

   if (! m_incompressible) {
      Print() << " Density: " << DENSITY << "\n";
      stateComponents.emplace_back(DENSITY,"density");
      Print() << " First species: " << FIRSTSPEC << "\n";
      Vector<std::string> names;
      pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(names);
      for (int n = 0; n < NUM_SPECIES; n++ ) {
         stateComponents.emplace_back(FIRSTSPEC+n,"rho.Y("+names[n]+")");
      }
      Print() << " Enthalpy: " << RHOH << "\n";
      stateComponents.emplace_back(RHOH,"rhoh");
      Print() << " Temperature: " << TEMP << "\n";
      stateComponents.emplace_back(TEMP,"temp");
      Print() << " thermo. pressure: " << RHORT << "\n";
      stateComponents.emplace_back(RHORT,"RhoRT");
#ifdef PELE_USE_EFIELD
      Print() << " nE: " << NE << "\n";
      stateComponents.emplace_back(NE,"nE");
      Print() << " PhiV: " << PHIV << "\n";
      stateComponents.emplace_back(PHIV,"PhiV");
#endif
   }

   if (m_nAux > 0) {
      Print() << " First passive scalar: " << FIRSTAUX << "\n";
      for (int n = 0; n < m_nAux; n++ ) {
         stateComponents.emplace_back(FIRSTAUX+n,"Aux_"+std::to_string(n));
      }
   }

   Print() << " => Total number of state variables: " << NVAR << "\n";
   Print() << PrettyLine;
   Print() << "\n";

   //----------------------------------------------------------------
   // Set advection/diffusion types
   m_AdvTypeState.resize(NVAR);
   m_DiffTypeState.resize(NVAR);

   // Velocity - follow incflo
   for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
#ifdef AMREX_USE_EB
      m_AdvTypeState[VELX+idim] = 1;   // Conservative
#else
      m_AdvTypeState[VELX+idim] = 0;   // NonConservative
#endif
      m_DiffTypeState[VELX+idim] = 1;  // Diffusive
   }

   if (! m_incompressible) {
      m_AdvTypeState[DENSITY] = 1;
      m_DiffTypeState[DENSITY] = 0;
      for (int n = 0; n < NUM_SPECIES; ++n) {
         m_AdvTypeState[FIRSTSPEC+n] = 1;
         m_DiffTypeState[FIRSTSPEC+n] = 1;
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
   }

   //----------------------------------------------------------------
   // Typical values container
   typical_values.resize(NVAR,-1.0);
}

void PeleLM::derivedSetup()
{
   BL_PROFILE_VAR("PeleLM::derivedSetup()", derivedSetup);

   if (!m_incompressible) {

      // Get species names
      Vector<std::string> spec_names;
      pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);

      // Set species mass fractions
      Vector<std::string> var_names_massfrac(NUM_SPECIES);
      for (int n = 0 ; n < NUM_SPECIES; n++) {
         var_names_massfrac[n] = "Y("+spec_names[n]+")";
      }
      derive_lst.add("mass_fractions",IndexType::TheCellType(),NUM_SPECIES,
                     var_names_massfrac,pelelm_dermassfrac,the_same_box);

   }

   // Cell average pressure
   derive_lst.add("avg_pressure",IndexType::TheCellType(),1,pelelm_deravgpress,the_same_box);

   // Vorticity magnitude
   derive_lst.add("mag_vort",IndexType::TheCellType(),1,pelelm_dermgvort,grow_box_by_two);

   // Kinetic energy
   derive_lst.add("kinetic_energy",IndexType::TheCellType(),1,pelelm_derkineticenergy,the_same_box);

#ifdef PELE_USE_EFIELD
   // Charge distribution
   derive_lst.add("chargedistrib",IndexType::TheCellType(),1,pelelm_derchargedist,the_same_box);

   // Electric field
   derive_lst.add("efieldx",IndexType::TheCellType(),1,pelelm_derefx,grow_box_by_one);
#if (AMREX_SPACEDIM > 1)
   derive_lst.add("efieldy",IndexType::TheCellType(),1,pelelm_derefy,grow_box_by_one);
#if (AMREX_SPACEDIM > 2)
   derive_lst.add("efieldz",IndexType::TheCellType(),1,pelelm_derefz,grow_box_by_one);
#endif
#endif
   // Lorentz forces
   derive_lst.add("LorentzFx",IndexType::TheCellType(),1,pelelm_derLorentzx,grow_box_by_one);
#if (AMREX_SPACEDIM > 1)
   derive_lst.add("LorentzFy",IndexType::TheCellType(),1,pelelm_derLorentzy,grow_box_by_one);
#if (AMREX_SPACEDIM > 2)
   derive_lst.add("LorentzFz",IndexType::TheCellType(),1,pelelm_derLorentzz,grow_box_by_one);
#endif
#endif
#endif

}

void PeleLM::taggingSetup()
{
   BL_PROFILE_VAR("PeleLM::taggingSetup()", taggingSetup);

   std::string amr_prefix = "amr";
   ParmParse ppamr(amr_prefix);

   Vector<std::string> refinement_indicators;
   ppamr.queryarr("refinement_indicators",refinement_indicators,0,ppamr.countval("refinement_indicators"));
   for (int n = 0; n<refinement_indicators.size(); ++n) {
      std::string ref_prefix = amr_prefix + "." + refinement_indicators[n];
      ParmParse ppr(ref_prefix);

      // Tag a given box
      RealBox realbox;
      if (ppr.countval("in_box_lo")) {
         Vector<Real> box_lo(AMREX_SPACEDIM);
         Vector<Real> box_hi(AMREX_SPACEDIM);
         ppr.getarr("in_box_lo",box_lo,0,box_lo.size());
         ppr.getarr("in_box_hi",box_hi,0,box_hi.size());
         realbox = RealBox(&(box_lo[0]),&(box_hi[0]));
      }

      AMRErrorTagInfo info;

      if (realbox.ok()) {
         info.SetRealBox(realbox);
      }

      if (ppr.countval("start_time") > 0) {
         Real min_time; ppr.get("start_time",min_time);
         info.SetMinTime(min_time);
      }

      if (ppr.countval("end_time") > 0) {
         Real max_time; ppr.get("end_time",max_time);
         info.SetMaxTime(max_time);
      }

      if (ppr.countval("max_level") > 0) {
         int tag_max_level; ppr.get("max_level",tag_max_level);
         info.SetMaxLevel(tag_max_level);
      }

      bool itexists = false;
      if (ppr.countval("value_greater")) {
         Real value; ppr.get("value_greater",value);
         std::string field; ppr.get("field_name",field);
         errTags.push_back(AMRErrorTag(value,AMRErrorTag::GREATER,field,info));
         itexists = derive_lst.canDerive(field) || isStateVariable(field);
      } else if (ppr.countval("value_less")) {
         Real value; ppr.get("value_less",value);
         std::string field; ppr.get("field_name",field);
         errTags.push_back(AMRErrorTag(value,AMRErrorTag::LESS,field,info));
         itexists = derive_lst.canDerive(field) || isStateVariable(field);
      } else if (ppr.countval("vorticity_greater")) {
         Real value; ppr.get("vorticity_greater",value);
         const std::string field="mag_vort";
         errTags.push_back(AMRErrorTag(value,AMRErrorTag::VORT,field,info));
         itexists = derive_lst.canDerive(field) || isStateVariable(field);
      } else if (ppr.countval("adjacent_difference_greater")) {
         Real value; ppr.get("adjacent_difference_greater",value);
         std::string field; ppr.get("field_name",field);
         errTags.push_back(AMRErrorTag(value,AMRErrorTag::GRAD,field,info));
         itexists = derive_lst.canDerive(field) || isStateVariable(field);
      } else if (realbox.ok()) {
        errTags.push_back(AMRErrorTag(info));
      } else {
        Abort(std::string("Unrecognized refinement indicator for " + refinement_indicators[n]).c_str());
      }

      if ( !itexists ) {
         amrex::Error("PeleLM::taggingSetup(): unknown variable field for criteria "+refinement_indicators[n]);
      }
   }
}

void PeleLM::resizeArray() {

   if (m_verbose) {
      Print() << " Initializing data for " << max_level+1 << " levels \n";
   }

   // State data
   m_leveldata_old.resize(max_level+1);
   m_leveldata_new.resize(max_level+1);
   m_leveldatareact.resize(max_level+1);
   m_halfTimeDensity.resize(max_level+1);
   if (max_level > 0) m_coveredMask.resize(max_level);
   m_baChem.resize(max_level+1);
   m_dmapChem.resize(max_level+1);
   m_baChemFlag.resize(max_level+1);

#ifdef PELE_USE_EFIELD
   m_leveldatanlsolve.resize(max_level+1);
#endif

   // Factory
   m_factory.resize(max_level+1);

   // Time
   m_t_old.resize(max_level+1);
   m_t_new.resize(max_level+1);

}
