
#include "SprayParticles.H"
#include "SprayInjection.H"
#include "pelelm_prob.H"

using namespace amrex;

class CounterFlowJet : public SprayJet
{
public:
  // Use default constructor
  CounterFlowJet(const std::string jet_name, const Geometry& geom)
    : SprayJet(jet_name, geom)
  {
  }

  bool get_new_particle(
    const Real time,
    const Real& phi_radial,
    const Real& cur_radius,
    Real& umag,
    Real& theta_spread,
    Real& phi_swirl,
    Real& dia_part,
    Real& T_part,
    Real* Y_part) override;
};

bool
CounterFlowJet::get_new_particle(
    const Real time,
    const Real& phi_radial,
    const Real& cur_radius,
    Real& umag,
    Real& theta_spread,
    Real& phi_swirl,
    Real& dia_part,
    Real& T_part,
    Real* Y_part)
{
  umag = m_jetVel;
  T_part = m_jetT;
  dia_part = m_dropDist->get_dia();
  phi_swirl = 0.;
  // Random spread angle
  Real tan_vel_comp = (2. * amrex::Random() - 1.) * 0.05;
  theta_spread = std::asin(tan_vel_comp);
  for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp) {
    Y_part[sp] = 0.;
  }
  Y_part[0] = 1.;
  return true;
}

bool
SprayParticleContainer::injectParticles(Real time,
                                        Real dt,
                                        int nstep,
                                        int lev,
                                        int finest_level,
                                        ProbParm const& prob_parm)
{
  if (lev != 0) {
      return false;
  }
  SprayJet* js = m_sprayJets[0].get();
  if (!js->jet_active(time)) {
    return false;
  }
  sprayInjection(time, js, dt, lev);
  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(const bool init_parts,
                                           ProbParm const& prob_parm)
{
  m_sprayJets.resize(1);
  std::string jet_name = "jet1";
  m_sprayJets[0] = std::make_unique<CounterFlowJet>(jet_name, Geom(0));
  m_sprayJets[0]->set_inj_proc(0);
  // This ensures the initial time step size stays reasonable
  m_injectVel = m_sprayJets[0]->jet_vel();
  // Start without any particles
  return;
}
