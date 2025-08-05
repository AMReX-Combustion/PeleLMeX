#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

namespace m2c = pele::physics::utilities::mks2cgs;
namespace c2m = pele::physics::utilities::cgs2mks;

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");
  auto eos = pele::physics::PhysicsType::eos();

  // Gas phase properties
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("T0_gas", PeleLM::prob_parm->T0_gas);
  pp.query("vel_gas", PeleLM::prob_parm->vel_gas);
  pp.query("N2_gas", PeleLM::prob_parm->Y_N2);
  pp.query("O2_gas", PeleLM::prob_parm->Y_O2);

  // Particle properties
  amrex::Real Re = -1.;
  amrex::Real drop_dia = 0.;
  amrex::Real T_d = 0.;
  pp.query("part_Re", Re);
  pp.query("part_dia", drop_dia);
  pp.query("part_temp", T_d);

  // Check that velocity is well defined
  if (Re > 0. && std::abs(PeleLM::prob_parm->vel_gas) > 0.) {
    amrex::Abort("Cannot specify droplet velocity and Reynolds number");
  }

  // Calculate transport properties from Simple transport model
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[N2_ID] = PeleLM::prob_parm->Y_N2;
  massfrac[O2_ID] = PeleLM::prob_parm->Y_O2;
  amrex::Real T_g = PeleLM::prob_parm->T0_gas;
  amrex::Real T_eff = (2. * T_d + T_g) / 3.;
  amrex::Real p_cgs = m2c::P(PeleLM::prob_parm->P_mean);
  amrex::Real rho_cgs = 0.;
  eos.PYT2R(p_cgs, massfrac, T_eff, rho_cgs);

  const bool get_xi = false, get_mu = true, get_lam = false, get_Ddiag = false,
             get_chi = false;
  amrex::Real dummy_xi, dummy_chi_mix;
  amrex::Real mu_cgs, lambda_cgs = 0.;
  amrex::Real rhoDi_cgs[NUM_SPECIES] = {0.0};
  auto trans = pele::physics::PhysicsType::transport();
  const auto* trans_parm = &(PeleLM::trans_parms.host_parm());
  trans.transport(
    get_xi, get_mu, get_lam, get_Ddiag, get_chi, T_eff, rho_cgs, massfrac,
    rhoDi_cgs, &dummy_chi_mix, mu_cgs, dummy_xi, lambda_cgs, trans_parm);

  amrex::Real mu = c2m::Mu(mu_cgs);
  amrex::Real rho = c2m::Rho(rho_cgs);

  // Set umax from Re
  if (Re > 0.) {
    // Get gas velocity from Re
    PeleLM::prob_parm->vel_gas = mu * Re / (rho * drop_dia);
  } else {
    Re = PeleLM::prob_parm->vel_gas * rho * drop_dia / mu;
  }

  std::ofstream ofs("ic.txt", std::ofstream::out);
  amrex::Print(ofs) << "Re = " << Re << "\n"
                    << "vel_gas = " << prob_parm->vel_gas << "\n"
                    << "mu_r = " << mu << "\n"
                    << "rho_r = " << rho << "\n"
                    << "drop_dia = " << drop_dia << std::endl;
  ofs.close();
}

void
PeleLM::freeProbParm()
{
}
