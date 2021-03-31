#include <PeleLM.H>
#include <PeleLM_K.H>
#include <DiffusionOp.H>

using namespace amrex;

void PeleLM::poissonSolveEF(const TimeStamp &a_time)
{
   BL_PROFILE_VAR("PeleLM::poissonSolveEF()", poissonSolveEF);

   // Get the phiV BCRec
   auto bcRecPhiV = fetchBCRecArray(PHIV,1);

   // Build Poisson RHS: charge distribution
   int nGhost = 0;
   Vector<std::unique_ptr<MultiFab>> rhsPoisson(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      rhsPoisson[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nGhost, MFInfo(), *m_factory[lev]));
      rhsPoisson[lev]->setVal(0.0); // TODO
   }

   // Solve for PhiV
   getDiffusionOp()->diffuse_scalar(getPhiVVect(a_time),0,
                                    GetVecOfConstPtrs(rhsPoisson),0,
                                    {},0,
                                    {},
                                    {},
                                    {},0,bcRecPhiV,1,1.0);
}
