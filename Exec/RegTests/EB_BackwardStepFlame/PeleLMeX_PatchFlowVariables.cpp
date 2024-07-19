
#include "PeleLMeX.H"
using namespace amrex;

void
patchFlowVariables(
  const amrex::Geometry& geom, ProbParm const& lprobparm, amrex::MultiFab& a_mf)
{

  amrex::Print() << "\nPatching flow variables..";
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
    geom.ProbLoArray();
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();

  for (amrex::MFIter mfi(a_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.tilebox();
    auto const& temp_arr = a_mf.array(mfi, TEMP);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      amrex::Real x[3] = {
        prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
        prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
        prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};

      amrex::ignore_unused(x);
      /*User can define how to patch flow variables here.*/
      temp_arr(i, j, k) = lprobparm.T_mean;
    });
  }

  amrex::Print() << "Done\n";
}
