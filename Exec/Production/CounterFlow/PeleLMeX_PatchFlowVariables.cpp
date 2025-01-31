
#include "PeleLMeX.H"
using namespace amrex;

void
patchFlowVariables(
  const amrex::Geometry& geom, ProbParm const& lprobparm, amrex::MultiFab& a_mf)
{

  amrex::Print() << "\nPatching flow variables..";

  const amrex::Real* prob_lo = geom.ProbLo();
  const amrex::Real* prob_hi = geom.ProbHi();
  const amrex::Real* dx = geom.CellSize();

  for (amrex::MFIter mfi(a_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.tilebox();

    auto const& rho_arr = a_mf.array(mfi, DENSITY);
    auto const& rhoY_arr = a_mf.array(mfi, FIRSTSPEC);
    auto const& rhoH_arr = a_mf.array(mfi, RHOH);
    auto const& temp_arr = a_mf.array(mfi, TEMP);
    Real massfrac[NUM_SPECIES] = {0.0};

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      auto eos = pele::physics::PhysicsType::eos();
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      amrex::Real sumYs = 0.0;

      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = rhoY_arr(i, j, k, n);
#ifdef N2_ID
        if (n != N2_ID) {
          sumYs += massfrac[n];
        }
#endif
      }
#ifdef N2_ID
      massfrac[N2_ID] = 1.0 - sumYs;
#endif

      AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                   , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                   , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);

      AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];
                   , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
                   , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

      AMREX_D_TERM(const amrex::Real xc = prob_lo[0] + Lx / 2.0;
                   , const amrex::Real yc = prob_lo[1] + Ly / 2.0;
                   , const amrex::Real zc = prob_lo[2] + Lz / 2.0);

      amrex::Real radiusSq = AMREX_D_TERM(
        (x - xc) * (x - xc), +(y - yc) * (y - yc), +(z - zc) * (z - zc));

      amrex::Real radius = std::sqrt(radiusSq);
      amrex::Real mixingWidth = 0.1 * lprobparm.ignitSphereRad;
      amrex::Real mixingFunction =
        0.5 *
        (1.0 + std::tanh((lprobparm.ignitSphereRad - radius) / mixingWidth));
      temp_arr(i, j, k) = mixingFunction * lprobparm.ignitT +
                          (1.0 - mixingFunction) * lprobparm.T_inert;
    });
  }

  amrex::Print() << "Done\n";
}
