#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMFData.H>
#include <PelePhysics.H>

// -----------------------------------------------------------
// Search for the closest index in an array to a given value
// using the bisection technique.
// INPUTS/OUTPUTS:
// xtable(0:n-1) => array to search in (ascending order)
// n             => array size
// x             => x location
// idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
// -----------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
locate(const amrex::Real* xtable, const int n, const amrex::Real& x, int& idxlo)
{
  // If x is out of bounds, return boundary index
  if (x >= xtable[n - 1]) {
    idxlo = n - 1;
    return;
  }
  if (x <= xtable[0]) {
    idxlo = 0;
    return;
  }

  // Do the bisection
  idxlo = 0;
  int idxhi = n - 1;
  bool notdone = true;
  while (notdone) {
    if (idxhi - idxlo <= 1) {
      notdone = false;
    } else {
      const int idxmid = (idxhi + idxlo) / 2;
      if (x >= xtable[idxmid]) {
        idxlo = idxmid;
      } else {
        idxhi = idxmid;
      }
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void pelelm_initdata(int i, int j, int k,
                     int is_incompressible,
                     amrex::Array4<amrex::Real> const& state,
                     amrex::Array4<amrex::Real> const& aux,
                     amrex::GeometryData const& geomdata,
                     ProbParm const& prob_parm,
                     pele::physics::PMF::PmfData::DataContainer const * pmf_data)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    amrex::Real x[3] = {
      prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0] * 100.0,
      prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1] * 100.0,
      prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2] * 100.0};

    // Fill in the velocities
    amrex::Real uinterp[3] = {0.0};

    // Interpolation factors
    amrex::Real mod[3] = {0.0};
    int idx[3] = {0};
    int idxp1[3] = {0};
    amrex::Real slp[3] = {0.0};
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      mod[idim] = std::fmod(x[idim], prob_parm.Linput);
      locate(prob_parm.d_xarray, prob_parm.input_resolution, mod[idim], idx[idim]);
      idxp1[idim] = (idx[idim] + 1) % prob_parm.input_resolution;
      slp[idim] = (mod[idim] - prob_parm.d_xarray[idx[idim]]) / prob_parm.d_xdiff[idx[idim]];
    }

    int inSize = prob_parm.input_resolution;

    const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
    const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);
    const amrex::Real f2 = (1 - slp[0]) * slp[1] * (1 - slp[2]);
    const amrex::Real f3 = (1 - slp[0]) * (1 - slp[1]) * slp[2];
    const amrex::Real f4 = slp[0] * (1 - slp[1]) * slp[2];
    const amrex::Real f5 = (1 - slp[0]) * slp[1] * slp[2];
    const amrex::Real f6 = slp[0] * slp[1] * (1 - slp[2]);
    const amrex::Real f7 = slp[0] * slp[1] * slp[2];

    uinterp[0] =   prob_parm.d_uinput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0
                 + prob_parm.d_uinput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] * f1
                 + prob_parm.d_uinput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] * f2
                 + prob_parm.d_uinput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] * f3
                 + prob_parm.d_uinput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] * f4
                 + prob_parm.d_uinput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f5
                 + prob_parm.d_uinput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] * f6
                 + prob_parm.d_uinput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f7;

    uinterp[1] =   prob_parm.d_vinput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0
                 + prob_parm.d_vinput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] * f1
                 + prob_parm.d_vinput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] * f2
                 + prob_parm.d_vinput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] * f3
                 + prob_parm.d_vinput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] * f4
                 + prob_parm.d_vinput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f5
                 + prob_parm.d_vinput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] * f6
                 + prob_parm.d_vinput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f7;

    uinterp[2] =   prob_parm.d_winput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0
                 + prob_parm.d_winput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] * f1
                 + prob_parm.d_winput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] * f2
                 + prob_parm.d_winput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] * f3
                 + prob_parm.d_winput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] * f4
                 + prob_parm.d_winput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f5
                 + prob_parm.d_winput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] * f6
                 + prob_parm.d_winput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] * f7;

    //
    // Fill Velocity
    //
    state(i,j,k,VELX) = uinterp[0];
    state(i,j,k,VELY) = uinterp[1];
    state(i,j,k,VELZ) = uinterp[2];

    if ( is_incompressible ) {
       return;
    }

    state(i,j,k,TEMP) = prob_parm.T_mean;

    amrex::Real Yt[NUM_SPECIES];
    Yt[O2_ID] = 0.233;
    Yt[N2_ID] = 0.767;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2R(P_cgs, Yt, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,TEMP), Yt, state(i,j,k,RHOH));
    state(i,j,k,RHOH) *= 1.0e-4 * state(i,j,k,DENSITY);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i,j,k,FIRSTSPEC+n) = Yt[n] * state(i,j,k,DENSITY);
    }
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
}
#endif
