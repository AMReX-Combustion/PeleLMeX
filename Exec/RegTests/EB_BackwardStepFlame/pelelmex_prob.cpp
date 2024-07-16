#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::Patch_Ignition_Source(
  const amrex::Geometry& geom, ProbParm const& lprobparm, amrex::MultiFab& a_mf)
{
  // Patch ignition kernel
  if (lprobparm.ignite_flow) {
    amrex::Print() << "\nCalling Ignition kernel patching..\n";
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
      geom.ProbLoArray();
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
      geom.CellSizeArray();

    for (amrex::MFIter mfi(a_mf, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      auto const& rho_arr = a_mf.array(mfi, DENSITY);
      auto const& rhoY_arr = a_mf.array(mfi, FIRSTSPEC);
      auto const& rhoH_arr = a_mf.array(mfi, RHOH);
      auto const& temp_arr = a_mf.array(mfi, TEMP);
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          auto eos = pele::physics::PhysicsType::eos();
          const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
            prob_lo[0] + (i + 0.5) * dx[0], prob_lo[1] + (j + 0.5) * dx[1],
            prob_lo[2] + (k + 0.5) * dx[2])};
          /*Patch ignition kernel here.*/
        });
    }
    amrex::Print() << "Done\n";
  }
}

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");

  pp.query("T_mean", prob_parm->T_mean);
  pp.query("P_mean", prob_parm->P_mean);
  pp.query("Y_fuel", prob_parm->Y_fuel);
  pp.query("Y_oxid", prob_parm->Y_o2);
  pp.query("T_hot", prob_parm->T_hot);
  pp.query("T_wall", prob_parm->Twall);
  pp.query("meanFlowMag", prob_parm->meanFlowMag);
  pp.query("do_ignite", prob_parm->ignite_flow);

  if (prob_parm->ignite_flow) {
    pp.query("ignition_temperature", prob_parm->ignition_temperature);
  }
}
