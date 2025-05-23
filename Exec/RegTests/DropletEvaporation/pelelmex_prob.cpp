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
  pp.query("init_T", PeleLM::prob_parm->T0);
  pp.query("init_vel", PeleLM::prob_parm->vel);
  pp.query("init_N2", PeleLM::prob_parm->Y_N2);
  pp.query("init_O2", PeleLM::prob_parm->Y_O2);

  // Particle properties
  amrex::Real Re = -1.;
  amrex::Real drop_dia = 0.;
  amrex::Real T_d = 0.;
  pp.query("part_Re", Re);
  pp.query("part_dia", drop_dia);
  pp.query("part_temp", T_d);

  if (Re > 0. && std::abs(PeleLM::prob_parm->vel) > 0.) {
    amrex::Abort("Cannot specify droplet velocity and Reynolds number");
  }
  if (Re > 0.) {
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[N2_ID] = PeleLM::prob_parm->Y_N2;
    massfrac[O2_ID] = PeleLM::prob_parm->Y_O2;
    amrex::Real T_g = PeleLM::prob_parm->T0;
    amrex::Real T_eff = (2.*T_d + T_g) / 3.;
    amrex::Real p_cgs = m2c::P(PeleLM::prob_parm->P_mean);
    amrex::Real rho_cgs = 0.;
    eos.PYT2R(p_cgs, massfrac, T_eff, rho_cgs);

    // Calculate transport properties from Simple transport model
    bool FA = false;
    bool TR = true;
    amrex::Real dummy_xi, dummy_chi_mix;
    amrex::Real mu_cgs, lambda_cgs = 0.; 
    amrex::Real rhoDi_cgs[NUM_SPECIES], Di[NUM_SPECIES] = {0.0};
    auto trans = pele::physics::PhysicsType::transport();
    const auto* trans_parm = &(PeleLM::trans_parms.host_parm());
    trans.transport
    (
      FA, TR, TR, TR, FA, T_eff, 
      rho_cgs, massfrac, rhoDi_cgs, &dummy_chi_mix, 
      mu_cgs, dummy_xi, lambda_cgs, trans_parm
    );
    
    // Get gas velocity from Re
    amrex::Real umax = c2m::Mu(mu_cgs) * Re / (c2m::Rho(rho_cgs) * drop_dia);
    PeleLM::prob_parm->vel = umax;
    amrex::Print() << "Re = " << Re << "\n"
                   << "gas_vel = " << prob_parm->vel << "\n"
                   << "mu_cgs = " << mu_cgs << "\n"
                   << "mu = " << c2m::Mu(mu_cgs) << "\n"
                   << "rho_cgs = "<< rho_cgs << "\n"
                   << "rho = "<< c2m::Rho(rho_cgs) << "\n"
                   << "dia = "<< drop_dia << std::endl;
  }

}

void
PeleLM::freeProbParm()
{
}
