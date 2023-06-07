#ifndef PELELM_PROB_EB_H_
#define PELELM_PROB_EB_H_

#include "pelelm_prob_parm.H"
#include "PelePhysics.H"

struct EBhandle {
    AMREX_GPU_DEVICE
    void bctype_eb(const amrex::Real xEBface[AMREX_SPACEDIM],
                   amrex::Real & EBflagType,
                   amrex::GeometryData const& /*geomdata*/,
                   ProbParm const& /*prob_parm*/) const
    {
        if ( xEBface[0] > 0.0098 ) {
          EBflagType = 1.0;
        } else {
          EBflagType = 0.0;
        }
    }

    AMREX_GPU_DEVICE
    void bcnormal_eb(const amrex::Real xEBface[AMREX_SPACEDIM],
                     amrex::Real s_ext[NVAR],
                     const amrex::Real /*time*/,
                     amrex::GeometryData const& /*geomdata*/,
                     ProbParm const& prob_parm) const
    {
        // Isothermal wall
        s_ext[TEMP] = prob_parm.Twall;

        // Inflow
        if (std::abs(xEBface[0]) <= 0.0015) {
            auto eos = pele::physics::PhysicsType::eos();
            s_ext[VELX] = 0.0;
            s_ext[VELY] = prob_parm.EBinflow_vel;
            amrex::Real massfrac[NUM_SPECIES] = {0.0};
            massfrac[CH4_ID] = prob_parm.EBinflow_Yfuel;
            massfrac[O2_ID] = (1.0 - massfrac[CH4_ID]) * 0.233;
            massfrac[N2_ID] =  1.0 - massfrac[CH4_ID] - massfrac[O2_ID];

            amrex::Real P_cgs = prob_parm.P_mean * 10.0;

            s_ext[TEMP] = prob_parm.EBinflow_T;
            amrex::Real rho_cgs, RhoH_temp;
            eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
            s_ext[DENSITY] = rho_cgs * 1.0e3;

            eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
            s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

            for (int n = 0; n < NUM_SPECIES; n++) {
                s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
            }
        }
    }
};
#endif
