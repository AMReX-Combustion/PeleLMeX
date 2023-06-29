#include <PeleLM.H>
#include <AMReX_ParmParse.H>


void PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  amrex::Real phi_in         = 0.0;
  amrex::Real vonKarman      = 0.0;
  amrex::Real channel_height = 0.0;
  amrex::Real U_b            = 0.0;
  amrex::Real T_in           = 0.0;
  amrex::Real P_mean         = 0.0;
  amrex::Real moles_NH3      = 0.0;
  amrex::Real moles_H2       = 0.0;
  amrex::Real moles_N2_fuel  = 0.0;
  amrex::Real zmax           = 0.0;
  amrex::Real zmin           = 0.0;
  amrex::Real molefrac[NUM_SPECIES] = {{0.0}};
  amrex::Real massfrac[NUM_SPECIES] = {{0.0}};

  pp.query("P_mean"       , P_mean);
  pp.query("T_in"         , T_in);
  pp.query("T_wall"       , PeleLM::prob_parm->T_wall);
  pp.query("U_b"          , U_b);
  pp.query("phi_in"       , phi_in);
  pp.query("channel_height", channel_height);
  
  // pp.query("U_crossflow"  ,PeleLM::prob_parm->U_crossflow);
  // pp.query("z_height_tanh",PeleLM::prob_parm->z_height_tanh);
  // pp.query("delta"        ,PeleLM::prob_parm->delta);

  pp.query("coflow_start_diam",PeleLM::prob_parm->coflow_start_diam);
  pp.query("U_coflow"         ,PeleLM::prob_parm->U_coflow);

  pp.query("z_height_tanh",PeleLM::prob_parm->z_height_tanh);
  pp.query("delta"        ,PeleLM::prob_parm->delta);
  pp.query("U_crossflow"  ,PeleLM::prob_parm->U_crossflow);

  pp.query("moles_NH3", moles_NH3);
  pp.query("moles_H2" , moles_H2);
  pp.query("moles_N2" , moles_N2_fuel);

  pp.query("vonKarman",vonKarman);
  pp.query("B"        ,PeleLM::prob_parm->B);

  amrex::ParmParse ppEB("EB");
  ppEB.query("in_diam",PeleLM::prob_parm->bluff_inner_diam);

  // PeleLM::prob_parm->phi_in = phi_in;
  // PeleLM::prob_parm->moles_NH3 = moles_NH3;
  // PeleLM::prob_parm->moles_H2 = moles_H2;
  // PeleLM::prob_parm->moles_N2_fuel = moles_N2_fuel;
  
  PeleLM::prob_parm->one_over_vonKarman = 1./vonKarman;
  PeleLM::prob_parm->channel_height = channel_height;
  PeleLM::prob_parm->U_b = U_b;
  PeleLM::prob_parm->T_in = T_in;
  PeleLM::prob_parm->P_mean = P_mean;

  // ------ Initializing NH3-H2-N2 premixed composition -------

  // amrex::Real equivRatio = phi_in;
  // amrex::Real a;
  // amrex::Real stoich_mole_frac = (3.0*moles_NH3+2.0*moles_H2)/4.0;

  // a = stoich_mole_frac/equivRatio;

  // amrex::Real moles_N2  = moles_N2_fuel+ a * 0.79/0.21; //Nitrogen from fuel blend + air
  // amrex::Real moles_O2  = a;
  // amrex::Real sum_oxi  = moles_NH3 + moles_H2 + moles_N2 + moles_O2;

  // molefrac[NH3_ID]  = moles_NH3 / sum_oxi;
  // molefrac[H2_ID]   = moles_H2  / sum_oxi;
  // molefrac[N2_ID]   = moles_N2  / sum_oxi;
  // molefrac[O2_ID]   = moles_O2  / sum_oxi;

  // ------ Initializing H2 premixed composition -------
  amrex::Real a = 0.5;
  molefrac[O2_ID] = 1.0 / ( 1.0 + phi_in / a + 0.79 / 0.21 );
  molefrac[H2_ID] = phi_in * molefrac[O2_ID] / a;
  molefrac[N2_ID] = 1.0 - molefrac[O2_ID] - molefrac[H2_ID];

  // molefrac[N2_ID]   = 0.79;
  // molefrac[O2_ID]   = 0.21;

  auto eos = pele::physics::PhysicsType::eos();

  for (int n = 0; n < NUM_SPECIES; n++)
    massfrac[n] = 0.0;

  eos.X2Y(molefrac,massfrac);

  for (int n = 0; n < NUM_SPECIES; n++){
    (PeleLM::prob_parm->Ys)[n] = massfrac[n];
  }

  amrex::Real rho_cgs;
  eos.PYT2R(P_mean*10., massfrac, T_in, rho_cgs);

  // ------ Initializing viscosity -------
  // Get the transport data pointer
  auto const* ltransparm = trans_parms.device_trans_parm();
  auto trans = pele::physics::PhysicsType::transport();

  const bool wtr_get_xi = true;
  const bool wtr_get_mu = true;
  const bool wtr_get_lam = true;
  const bool wtr_get_Ddiag = true;
  const bool wtr_get_chi = true;
  amrex::Real muloc, xiloc, lamloc;
  amrex::Real Ddiag[NUM_SPECIES] = {0.0};
  amrex::Real chi_mix[NUM_SPECIES] = {0.0};

  // trans.transport(wtr_get_xi, wtr_get_mu, wtr_get_lam, wtr_get_Ddiag, wtr_get_chi, T_in, rho_cgs, massfrac,
  //           Ddiag, chi_mix, muloc, xiloc, lamloc, ltransparm);

  muloc = 1.751067e-04 ;

  amrex::Real mu = muloc/rho_cgs*1.e-4;
  PeleLM::prob_parm->mu = mu;
  // amrex::Print() << "---> Kinematic viscosity [m2/s] = " << mu << std::endl;


  amrex::Real rho_si = rho_cgs*1.e+3;
  amrex::Real half_height = channel_height/2.;

  //Reynolds number based on bulk velocity and half air passage width
  amrex::Real Re = U_b*half_height/mu;
  
  //Friction coefficient taken from Schlichting
  amrex::Real Cf = pow(2.*log10(Re)-0.65,-2.3);

  amrex::Real wall_stress = Cf*0.5*rho_si*pow(U_b,2);
  amrex::Real friction_vel = sqrt(wall_stress/rho_si);

  PeleLM::prob_parm->friction_vel = friction_vel;

}
