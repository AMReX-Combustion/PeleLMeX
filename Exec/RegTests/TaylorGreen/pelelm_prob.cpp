#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
    amrex::ParmParse pp("prob");
    
    pp.query("T_mean", prob_parm->T_mean);
    pp.query("P_mean", prob_parm->P_mean);
    pp.query("reynolds", prob_parm->reynolds);
    pp.query("mach", prob_parm->mach);
    pp.query("prandtl", prob_parm->prandtl);
    pp.query("convecting", prob_parm->convecting);
    pp.query("omega_x", prob_parm->omega_x);
    pp.query("omega_y", prob_parm->omega_y);
    pp.query("omega_z", prob_parm->omega_z);

    prob_parm->L = 0.01 / M_PI;

    // Mixture composition
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real P_cgs = prob_parm->P_mean * 10.0;

    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, prob_parm->T_mean, rho_cgs);
    prob_parm->rho0 = rho_cgs * 1.0e3;

    // Velocity (based on Mach number)
    amrex::Real cs;
    eos.RTY2Cs(rho_cgs, prob_parm->T_mean, massfrac, cs);
    prob_parm->v0 = prob_parm->mach * cs * 0.01;

    // Transport
    amrex::Real cpmix_cgs;
    eos.TY2Cp(prob_parm->T_mean, massfrac, cpmix_cgs);
    auto& trans_parm = PeleLM::trans_parms.host_trans_parm();
    trans_parm.const_bulk_viscosity = 0.0;
    trans_parm.const_diffusivity = 0.0;
    trans_parm.const_viscosity = rho_cgs * prob_parm->v0 * 100
                                 * prob_parm->L * 100 / prob_parm->reynolds;
    trans_parm.const_conductivity = trans_parm.const_viscosity * cpmix_cgs / prob_parm->prandtl;
    PeleLM::trans_parms.sync_to_device();

    // Debug statement
    amrex::Print() << " ### Taylor-Green params ############ \n";
    amrex::Print() << " Mach : " << prob_parm->mach << "\n";
    amrex::Print() << " Reynolds : " << prob_parm->reynolds << "\n";
    amrex::Print() << " v0 : " << prob_parm->v0 << "\n";
    amrex::Print() << " L : " << prob_parm->L << "\n";
    amrex::Print() << " t_s : " << prob_parm->L/prob_parm->v0 << "\n";
    amrex::Print() << " mu [CGS] : " << trans_parm.const_viscosity << "\n";
    amrex::Print() << " #################################### \n";
}
