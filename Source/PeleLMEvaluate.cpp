#include <PeleLM.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void PeleLM::Evaluate() {
   BL_PROFILE("PeleLM::Evaluate()");

   //----------------------------------------------------------------
   // Check that requested evaluate entries exist and determine the size
   // of the container and entries names
   int ncomp = 0;
   Vector<std::string> plt_VarsName;
   for (int ivar = 0; ivar < m_evaluatePlotVarCount; ivar++ ) {
      Print() << m_evaluatePlotVars[ivar] << "\n";
      bool itexists =    derive_lst.canDerive(m_evaluatePlotVars[ivar])
                      || evaluate_lst.canDerive(m_evaluatePlotVars[ivar])
                      || isStateVariable(m_evaluatePlotVars[ivar]);
      if ( !itexists ) {
         amrex::Error("PeleLM::evaluate(): unknown variable: "+m_evaluatePlotVars[ivar]);
      }
      if ( derive_lst.canDerive(m_evaluatePlotVars[ivar]) ) {
         const PeleLMDeriveRec* rec = derive_lst.get(m_evaluatePlotVars[ivar]);
         ncomp += rec->numDerive();
         for (int dvar = 0; dvar < rec->numDerive(); dvar++ ) {
            plt_VarsName.push_back(rec->variableName(dvar));
         }
      } else if ( evaluate_lst.canDerive(m_evaluatePlotVars[ivar]) ) {
         const PeleLMDeriveRec* rec = evaluate_lst.get(m_evaluatePlotVars[ivar]);
         ncomp += rec->numDerive();
         for (int dvar = 0; dvar < rec->numDerive(); dvar++ ) {
            plt_VarsName.push_back(rec->variableName(dvar));
         }
      } else if ( isStateVariable(m_evaluatePlotVars[ivar]) ) {
         ncomp += 1;
         plt_VarsName.push_back(m_evaluatePlotVars[ivar]);
      }
   }

   //----------------------------------------------------------------
   // Define the outgoing container
   Vector<MultiFab> mf_plt(finest_level + 1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      mf_plt[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
   }

   //----------------------------------------------------------------
   // Fill the outgoing container
   int cnt = 0;
   for (int ivar = 0; ivar < m_evaluatePlotVarCount; ivar++ ) {
      int cntIncr = 0;

      // Evaluate function calls actual PeleLM::Evolve pieces and may require
      // the entire multi-level hierarchy
      if ( evaluate_lst.canDerive(m_evaluatePlotVars[ivar]) ) {
         MLevaluate(GetVecOfPtrs(mf_plt),cnt,cntIncr,m_evaluatePlotVars[ivar]);

      // Regular derived functions and State entries are called on a per level basis
      // derive function can handle both derived and state entries
      } else if (    derive_lst.canDerive(m_evaluatePlotVars[ivar])
                  || isStateVariable(m_evaluatePlotVars[ivar]) ) {
         for (int lev = 0; lev <= finest_level; ++lev) {
            std::unique_ptr<MultiFab> mf;
            mf = derive(m_evaluatePlotVars[ivar], m_cur_time, lev, 0);
            MultiFab::Copy(mf_plt[lev], *mf, 0, cnt, mf->nComp(), 0);
            cntIncr = mf->nComp();
         }
      }
      cnt += cntIncr;
   }

   //----------------------------------------------------------------
   // Write the evaluated variables to disc
   Vector<int> istep(finest_level + 1, 0);

   std::string plotfilename = "pltEvaluate";
   amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf_plt),
                                  plt_VarsName, Geom(), m_cur_time, istep, refRatio());
}

void
PeleLM::MLevaluate(const Vector<MultiFab *> &a_MFVec,
                   int a_comp,
                   int &nComp,
                   const std::string &a_var)
{

    // This function manually maps the evaluate variables to the function calls
    // used in PeleLM:::Evolve

    if ( a_var == "divU" ) {
        int is_initialization = 0;             // No, use IRR
        int computeDiffusionTerm = 1;          // Needed here
        int do_avgDown = 1;                    // Always

        // Light version of the diffusion data container
        std::unique_ptr<AdvanceDiffData> diffData;
        diffData.reset(new AdvanceDiffData(finest_level, grids, dmap, m_factory,
                       m_nGrowAdv, m_use_wbar, is_initialization));
        calcDivU(is_initialization,computeDiffusionTerm,do_avgDown,AmrNewTime,diffData);
        for (int lev = 0; lev <= finest_level; ++lev) {
           auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
           MultiFab::Copy(*a_MFVec[lev],ldata_p->divu,0,a_comp,1,0);
        }
        nComp = 1;
    } else if ( a_var == "diffTerm" ) {
        // Use the diffusion data holder, get diffusivity and calc D
        // Finally, copy into a_MFVec
        std::unique_ptr<AdvanceDiffData> diffData;
        diffData.reset(new AdvanceDiffData(finest_level, grids, dmap, m_factory, m_nGrowAdv, m_use_wbar));
        calcDiffusivity(AmrNewTime);
        computeDifferentialDiffusionTerms(AmrNewTime,diffData);
        for (int lev = 0; lev <= finest_level; ++lev) {
           MultiFab::Copy(*a_MFVec[lev],diffData->Dnp1[lev],0,a_comp,NUM_SPECIES+2,0);
        }
        nComp = NUM_SPECIES+2;
    } else if ( a_var == "advTerm" ) {
        // TODO
        Abort("advTerm not available yet, soon hopefully");
    } else if ( a_var == "instRR" ) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            std::unique_ptr<MultiFab> I_RR = std::make_unique<MultiFab> (*a_MFVec[lev],amrex::make_alias,a_comp,NUM_SPECIES);
            computeInstantaneousReactionRate(lev, AmrNewTime, I_RR.get());
        }
        nComp = NUM_SPECIES;
    } else if ( a_var == "transportCC" ) {
        // Cell-centered transport coefficients functions go through the level
        // data container. Simply copy once the later has been filled.
        calcViscosity(AmrNewTime);
        calcDiffusivity(AmrNewTime);
        for (int lev = 0; lev <= finest_level; ++lev) {
           auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
           MultiFab::Copy(*a_MFVec[lev],ldata_p->diff_cc,0,a_comp,NUM_SPECIES+1,0);
           MultiFab::Copy(*a_MFVec[lev],ldata_p->visc_cc,0,a_comp+NUM_SPECIES+1,1,0);
        }
        nComp = NUM_SPECIES+2;
    }
}
