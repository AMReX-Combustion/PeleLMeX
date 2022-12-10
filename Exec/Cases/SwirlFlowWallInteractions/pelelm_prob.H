#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_Y_from_Phi(ProbParm const& prob_parm,
               const amrex::Real phi_in,
               amrex::Real Y[])
{
   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real Xt[NUM_SPECIES] = {0.0};
   amrex::Real a = 0.0;
   if ( prob_parm.fuelID == CH4_ID ) {
      a = 2.0;
   }
   Xt[prob_parm.oxidID] = 1.0 / ( 1.0 + phi_in / a + 0.79 / 0.21 );
   Xt[prob_parm.fuelID] = phi_in * Xt[prob_parm.oxidID] / a;
   Xt[prob_parm.bathID] = 1.0 - Xt[prob_parm.oxidID] - Xt[prob_parm.fuelID];

   eos.X2Y(Xt,Y);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
set_Y_from_Ksi(ProbParm const& prob_parm,
               const amrex::Real ksi_in,
               amrex::Real Y[])
{
   amrex::Real phi = ksi_in / ( 1.0 - amrex::min(ksi_in,0.999999) ) * (1.0 - prob_parm.Zst ) / prob_parm.Zst;

   set_Y_from_Phi(prob_parm,phi,Y);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k,
                 int is_incompressible,
                 amrex::Array4<amrex::Real> const& state,
                 amrex::Array4<amrex::Real> const& aux,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    AMREX_D_TERM(,
                 const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1]);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real pert = 0.0;
    amrex::Real y1, y2, z1, z2;
    amrex::Real massfrac_eq[NUM_SPECIES] = {0.0};

    amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real molefrac_air[NUM_SPECIES] = {0.0};

    amrex::Real Do_swirler = 0.0324;
    amrex::Real Di_swirler = 0.0202;

    amrex::Real Ro_swirler = Do_swirler/2;
    amrex::Real Ri_swirler = Di_swirler/2;
    amrex::Real R = Ri_swirler/Ro_swirler;
    amrex::Real x0 = 0.16;
    amrex::Real y0 = 0.068;
    amrex::Real H = 0.005;
    amrex::Real D_coflow = Do_swirler;
    amrex::Real Ro_coflow = prob_hi[0] ;
    amrex::Real Ri_coflow = Ro_swirler;
    amrex::Real c1 = -0.0067;
    amrex::Real c2 = 3.39e-7;
    amrex::Real c3 = -21.081;
    amrex::Real A = 0.0096;
    amrex::Real B = -3.869e-11;
    amrex::Real F_theta = -7.1884e5;

    massfrac_air[prob_parm.oxidID] = 0.233;
    massfrac_air[prob_parm.bathID] = 0.767;

    molefrac_air[prob_parm.oxidID] = 0.21;
    molefrac_air[prob_parm.bathID] = 0.79;

    amrex::Real y0_1;
    amrex::Real z0_1 = 0.0;
    amrex::Real c1_1 = 50;
    amrex::Real c2_1 = 50;

    amrex::Real y0_2;
    amrex::Real z0_2 = 0.01;
    amrex::Real c1_2 = 150;
    amrex::Real c2_2 = 150;

    y0_1 = 0.112;
    y0_2 = 0.112;

    amrex::Real delta = 0.02;

    amrex::Real cut_off_height = 0.04;

    amrex::Real jet_1 = 0.5*(1-((std::tanh((std::sqrt(std::pow(x-x0,2)+std::pow(y-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x-x0,2)+std::pow(y-y0,2))+(Do_swirler/2))/(H))));

    amrex::Real jet_2 = 0.5*(1-std::tanh((z-cut_off_height)/H));

    amrex::Real jet_total = jet_1 * jet_2;

    z1 = (z - prob_parm.standoff - 0.5*dx[2] + pert)*prob_parm.amplification;
    z2 = (z - prob_parm.standoff + 0.5*dx[2] + pert)*prob_parm.amplification;

    pele::physics::PMF::pmf(pmf_data, z1, z2, pmf_vals);

    state(i,j,k,TEMP) = (pmf_vals[0]*jet_total) + ((1-jet_total)*(310));

    for (int n = 0; n < NUM_SPECIES; ++n) {
        molefrac[n] = (pmf_vals[3 + n]*jet_total)+((1-jet_total)*(molefrac_air[n]));
    }

    eos.X2Y(molefrac, massfrac);

    state(i,j,k,VELX) = 0;
    state(i,j,k,VELY) = 0;
    state(i,j,k,VELZ) = 1.;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,TEMP), massfrac, state(i,j,k,RHOH));
    state(i,j,k,RHOH) *= 1.0e-4 * state(i,j,k,DENSITY);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
        state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
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
  const amrex::Real* prob_hi = geomdata.ProbHi();
  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

  amrex::Real Ro_swirler = 0.0162;
  amrex::Real Ri_swirler = 0.0101;
  amrex::Real R = Ri_swirler/Ro_swirler;
  amrex::Real D_jet = 0.02;
  amrex::Real x0 = 0.16;
  amrex::Real y0 = 0.068;
  amrex::Real H = 0.00001;
  amrex::Real Do_swirler = 2*Ro_swirler;
  amrex::Real Di_swirler = 2*Ri_swirler;
  amrex::Real D_coflow = Do_swirler;
  amrex::Real Ro_coflow = prob_hi[0] ;
  amrex::Real Ri_coflow = Ro_swirler;
  amrex::Real c1 = -0.0067;
  amrex::Real c2 = 3.39e-7;
  amrex::Real c3 = -21.081;
  amrex::Real A = 0.0096;
  amrex::Real B = -3.869e-11;
  amrex::Real F_theta = -7.1884e5;

  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
  massfrac_air[prob_parm.oxidID] = 0.233;
  massfrac_air[prob_parm.bathID] = 0.767;

  amrex::Real dt = 2e-5;

  int refresh_freq = 10;

  int n_const = 5;

  amrex::Real t0 = 0.0;

  amrex::Real A_fluctuation[200] = {-0.101095, 0.329953, 0.321893, -0.508907, 0.267847, 0.239758, 0.243639, 0.173476, 0.837604, 0.524968, 0.780876, 0.428267, -0.637454, -0.769445, -0.946879, 0.244172, 0.209526, 0.589716, 0.687176, 0.525888, -0.290012, -0.879359, -0.155726, 0.213109, 0.289708, -0.287101, 0.195927, -0.451034, -0.749719, -0.214496, 0.0347229, 0.149186, -0.884544, -0.643384, 0.640279, 0.383303, 0.596374, -0.116082, -0.44322, 0.433978, -0.591115, -0.662344, -0.137756, -0.228568, -0.431789, -0.0846348, -0.984396, 0.777738, -0.494919, 0.70278, 0.303626, 0.215069, 0.823421, -0.8521, -0.571822, 0.113128, -0.139201, 0.624105, 0.662094, 0.11108, -0.590391, -0.303183, -0.739735, -0.474935, 0.0534326, 0.900544, 0.908368, -0.350193, -0.215538, -0.534852, -0.916216, 0.193347, -0.197196, -0.053971, 0.964779, 0.371015, 0.861394, 0.980383, 0.148753, -0.633525, 0.683163, -0.547621, 0.581544, 0.506584, -0.399721, -0.990278, -0.380288, 0.461077, 0.633827, -0.718194, -0.427843, -0.956564, -0.0213777, -0.167578, -0.431499, -0.967945, -0.267034, -0.523131, -0.318139, 0.517428, -0.0579822, -0.234354, -0.289225, 0.744822, 0.711675, -0.324446, 0.115837, 0.573069, -0.344063, -0.73541, 0.939544, -0.6609, -0.283031, 0.521089, 0.845684, 0.317248, 0.530811, -0.534604, -0.221675, 0.164638, -0.252798, 0.350482, 0.208074, 0.725824, -0.817095, 0.776575, 0.757879, -0.0841292, -0.746555, -0.56026, -0.566701, 0.195463, 0.205386, 0.144073, -0.0597154, -0.082939, 0.819628, -0.943878, -0.50987, -0.524435, -0.679288, -0.570326, -0.185335, 0.0376813, 0.950763, -0.33965, -0.645071, 0.481574, 0.125746, 0.133255, -0.353788, 0.872947, -0.516263, 0.854287, 0.598771, -0.333358, 0.630862, 0.35665, 0.582513, 0.884307, 0.79639, -0.984189, 0.07977, 0.00177646, 0.159885, -0.979945, 0.918837, -0.0204875, -0.923824, -0.591033, 0.455077, -0.603111, -0.161358, -0.730257, 0.43457, -0.210595, -0.0699077, 0.789499, -0.729021, -0.944162, -0.077246, -0.0828092, 0.928785, 0.406491, -0.228522, 0.527557, -0.926867, -0.59766, -0.115793, 0.655646, -0.713353, -0.319403, 0.671457, 0.366417, 0.682373, -0.168658, 0.386472, 0.601211, 0.810855, 0.462648};

  amrex::Real B_fluctuation[200] = {-0.101095, 0.329953, 0.321893, -0.508907, 0.267847, 0.239758, 0.243639, 0.173476, 0.837604, 0.524968, 0.780876, 0.428267, -0.637454, -0.769445, -0.946879, 0.244172, 0.209526, 0.589716, 0.687176, 0.525888, -0.290012, -0.879359, -0.155726, 0.213109, 0.289708, -0.287101, 0.195927, -0.451034, -0.749719, -0.214496, 0.0347229, 0.149186, -0.884544, -0.643384, 0.640279, 0.383303, 0.596374, -0.116082, -0.44322, 0.433978, -0.591115, -0.662344, -0.137756, -0.228568, -0.431789, -0.0846348, -0.984396, 0.777738, -0.494919, 0.70278, 0.303626, 0.215069, 0.823421, -0.8521, -0.571822, 0.113128, -0.139201, 0.624105, 0.662094, 0.11108, -0.590391, -0.303183, -0.739735, -0.474935, 0.0534326, 0.900544, 0.908368, -0.350193, -0.215538, -0.534852, -0.916216, 0.193347, -0.197196, -0.053971, 0.964779, 0.371015, 0.861394, 0.980383, 0.148753, -0.633525, 0.683163, -0.547621, 0.581544, 0.506584, -0.399721, -0.990278, -0.380288, 0.461077, 0.633827, -0.718194, -0.427843, -0.956564, -0.0213777, -0.167578, -0.431499, -0.967945, -0.267034, -0.523131, -0.318139, 0.517428, -0.0579822, -0.234354, -0.289225, 0.744822, 0.711675, -0.324446, 0.115837, 0.573069, -0.344063, -0.73541, 0.939544, -0.6609, -0.283031, 0.521089, 0.845684, 0.317248, 0.530811, -0.534604, -0.221675, 0.164638, -0.252798, 0.350482, 0.208074, 0.725824, -0.817095, 0.776575, 0.757879, -0.0841292, -0.746555, -0.56026, -0.566701, 0.195463, 0.205386, 0.144073, -0.0597154, -0.082939, 0.819628, -0.943878, -0.50987, -0.524435, -0.679288, -0.570326, -0.185335, 0.0376813, 0.950763, -0.33965, -0.645071, 0.481574, 0.125746, 0.133255, -0.353788, 0.872947, -0.516263, 0.854287, 0.598771, -0.333358, 0.630862, 0.35665, 0.582513, 0.884307, 0.79639, -0.984189, 0.07977, 0.00177646, 0.159885, -0.979945, 0.918837, -0.0204875, -0.923824, -0.591033, 0.455077, -0.603111, -0.161358, -0.730257, 0.43457, -0.210595, -0.0699077, 0.789499, -0.729021, -0.944162, -0.077246, -0.0828092, 0.928785, 0.406491, -0.228522, 0.527557, -0.926867, -0.59766, -0.115793, 0.655646, -0.713353, -0.319403, 0.671457, 0.366417, 0.682373, -0.168658, 0.386472, 0.601211, 0.810855, 0.462648};

  int n = abs(time-t0)/(dt);
  int n_skip_lines = floor(n/refresh_freq);

  amrex::Real C[5];
  amrex::Real D[5];

  for (int i = 0; i < 5; i++) {
      C[i] = abs(A_fluctuation[i + (n_skip_lines * n_const)]);
      D[i] = B_fluctuation[i + (n_skip_lines * n_const)];
  }

  amrex::Real amp = 2;

  amrex::Real fluctuation = amp * ((std::sin(2 * Pi * (std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2)) - D[0] )* 1.5 / Ro_swirler))
                            +(std::sin(2 * Pi * (std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2)) - D[1] )* 1.75 / Ro_swirler))
                            +(std::sin(2 * Pi * (std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2)) - D[2] )* 0.79 / Ro_swirler))
                            +(std::sin(2 * Pi * (std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2)) - D[3] )* 1.1 / Ro_swirler))
                            +(std::sin(2 * Pi * (std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2)) - D[4] )* 0.95 / Ro_swirler)));

  amrex::Real vel_swirler_x = 3 * 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*(((x[1]-y0)/(std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))))*(F_theta)*(((std::pow((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))),2))/3)+(c1*(std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))))+((c2)/((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2)))))));

  amrex::Real vel_fluctuation_swirler_x = fluctuation * 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*(((x[1]-y0)/(std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2)))));

  amrex::Real vel_swirler_y = 3 * 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*((x[0]-x0)/(std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))))*(F_theta)*(((std::pow((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))),2))/3)+(c1*(std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))))+((c2)/((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))))));

  amrex::Real vel_fluctuation_swirler_y = fluctuation * 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*(((x[1]-y0)/(std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2)))));
  amrex::Real vel_jet_x = (fluctuation) *( 1 - (0.5*(1+(std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(D_jet/2))/(H)))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(D_jet/2))/(H))))));

  amrex::Real vel_jet_y = (fluctuation) *( 1 - (0.5*(1+(std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(D_jet/2))/(H)))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(D_jet/2))/(H))))));

  amrex::Real vel_swirler_z = 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*prob_parm.U_b*(c3)*((std::pow(((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2)))/(Ro_swirler)),2))-(((std::pow(R,2)-1)/(std::log(R)))*std::log((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2)))/Ro_swirler))-1);

  amrex::Real vel_fluctuation_swirler_z = fluctuation * 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))));

  amrex::Real vel_coflow_z = prob_parm.U_c * 0.5 * (1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(D_coflow/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(D_coflow/2))/(H))));

  amrex::Real vel_jet_z = (0.8*prob_parm.U_b + fluctuation) *( 1 - (0.5*(1+(std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(D_jet/2))/(H)))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(D_jet/2))/(H))))));

  s_ext[VELX] = -vel_swirler_x + vel_jet_x + vel_fluctuation_swirler_x;
  s_ext[VELY] = vel_swirler_y + vel_jet_y + vel_fluctuation_swirler_y;
  s_ext[VELZ] = (vel_swirler_z + vel_coflow_z + vel_jet_z + vel_fluctuation_swirler_z) * std::tanh(x[1]*200);

  s_ext[TEMP] = prob_parm.T_in;

  // Phi: 0.4 , X[CH4] = 0.3, X[H2] = 0.7

  massfrac_fuel[CH4_ID] = 0.01447;
  massfrac_fuel[H2_ID]  = 0.004244;
  massfrac_fuel[N2_ID]  = 0.752786;
  massfrac_fuel[O2_ID]  = 0.2285;

  amrex::Real Exponential_Spec = 0.5*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))));

  for (int n = 0; n < NUM_SPECIES; ++n) {
      massfrac[n] = (1.0 - Exponential_Spec) * massfrac_air[n] + Exponential_Spec * massfrac_fuel[n];
  }

  amrex::Real rho_cgs, P_cgs, RhoH_temp;
  P_cgs = prob_parm.P_mean * 10.0;

  auto eos = pele::physics::PhysicsType::eos();

  eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
  s_ext[DENSITY] = rho_cgs * 1.0e3;

  eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
  s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int /*i*/, int /*j*/, int /*k*/,
           amrex::Array4<amrex::Real> const& /*beta*/,
           amrex::GeometryData const& /*geomdata*/,
           amrex::Box const& /*dBox*/,
           const int  /*dir*/,
           const int  /*state_comp*/,
           const int  /*nComp*/)
{
}
#endif
