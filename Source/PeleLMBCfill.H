#include <PeleLM_Index.H>
#include <pelelm_prob.H>
#include <pmf_data.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Geometry.H>

using namespace amrex;

struct PeleLMCCFillExtDirVel
{

   ProbParm const* lprobparm;
   PmfData const* lpmfdata;
   const int m_nAux;

   AMREX_GPU_HOST
   constexpr PeleLMCCFillExtDirVel(ProbParm const* a_prob_parm, PmfData const* a_pmf_data, int a_nAux)
                                   : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux) {}

   AMREX_GPU_DEVICE
   void operator()(
      const amrex::IntVect& iv,
      amrex::Array4<amrex::Real> const& vel,
      const int dcomp,
      const int /*numcomp*/,
      amrex::GeometryData const& geom,
      const amrex::Real time,
      const amrex::BCRec* bcr,
      const int /*bcomp*/,
      const int orig_comp) const
   {

      // Get geometry data
      const int* domlo = geom.Domain().loVect();
      const int* domhi = geom.Domain().hiVect();
      const amrex::Real* prob_lo = geom.ProbLo();
      const amrex::Real* dx = geom.CellSize();
      const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(prob_lo[0] + (iv[0] + 0.5) * dx[0],
                                                          prob_lo[1] + (iv[1] + 0.5) * dx[1],
                                                          prob_lo[2] + (iv[2] + 0.5) * dx[2])};

      // Get the first velocity component BC
      const int* bc = bcr->data();

      // bcnormal handles all the state components at once
      amrex::Real s_ext[NVAR] = {0.0};

      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
         // Low
         if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

            // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
            bcnormal(x, m_nAux, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

            for (int n = 0; n < AMREX_SPACEDIM; n++) {
               vel(iv,n) = s_ext[VELX+n];
            }

         // High
         } else if ((bc[idir+AMREX_SPACEDIM] == amrex::BCType::ext_dir) and (iv[idir] > domhi[idir])) {

            // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
            bcnormal(x, m_nAux, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

            for (int n = 0; n < AMREX_SPACEDIM; n++) {
               vel(iv,n) = s_ext[VELX+n];
            }

         }
      }
   }
};

struct PeleLMCCFillExtDirSpec
{

   ProbParm const* lprobparm;
   PmfData const* lpmfdata;
   const int m_nAux;

   AMREX_GPU_HOST
   constexpr PeleLMCCFillExtDirSpec(ProbParm const* a_prob_parm, PmfData const* a_pmf_data, int a_nAux)
                                   : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux) {}

   AMREX_GPU_DEVICE
   void operator()(
      const amrex::IntVect& iv,
      amrex::Array4<amrex::Real> const& rhoY,
      const int dcomp,
      const int /*numcomp*/,
      amrex::GeometryData const& geom,
      const amrex::Real time,
      const amrex::BCRec* bcr,
      const int /*bcomp*/,
      const int orig_comp) const
   {

      // Get geometry data
      const int* domlo = geom.Domain().loVect();
      const int* domhi = geom.Domain().hiVect();
      const amrex::Real* prob_lo = geom.ProbLo();
      const amrex::Real* dx = geom.CellSize();
      const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(prob_lo[0] + (iv[0] + 0.5) * dx[0],
                                                          prob_lo[1] + (iv[1] + 0.5) * dx[1],
                                                          prob_lo[2] + (iv[2] + 0.5) * dx[2])};

      // Get the first species component BC
      const int* bc = bcr->data();

      // bcnormal handles all the state components at once
      amrex::Real s_ext[NVAR] = {0.0};

      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
         // Low
         if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

            // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
            bcnormal(x, m_nAux, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

            for (int n = 0; n < NUM_SPECIES; n++) {
               rhoY(iv,n) = s_ext[FIRSTSPEC+n];
            }

         // High
         } else if ((bc[idir+AMREX_SPACEDIM] == amrex::BCType::ext_dir) and (iv[idir] > domhi[idir])) {

            // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
            bcnormal(x, m_nAux, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

            for (int n = 0; n < NUM_SPECIES; n++) {
               rhoY(iv,n) = s_ext[FIRSTSPEC+n];
            }

         }
      }
   }
};

struct PeleLMCCFillExtDirDens
{

   ProbParm const* lprobparm;
   PmfData const* lpmfdata;
   const int m_nAux;

   AMREX_GPU_HOST
   constexpr PeleLMCCFillExtDirDens(ProbParm const* a_prob_parm, PmfData const* a_pmf_data, int a_nAux)
                                   : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux) {}

   AMREX_GPU_DEVICE
   void operator()(
      const amrex::IntVect& iv,
      amrex::Array4<amrex::Real> const& rho,
      const int dcomp,
      const int /*numcomp*/,
      amrex::GeometryData const& geom,
      const amrex::Real time,
      const amrex::BCRec* bcr,
      const int /*bcomp*/,
      const int orig_comp) const
   {

      // Get geometry data
      const int* domlo = geom.Domain().loVect();
      const int* domhi = geom.Domain().hiVect();
      const amrex::Real* prob_lo = geom.ProbLo();
      const amrex::Real* dx = geom.CellSize();
      const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(prob_lo[0] + (iv[0] + 0.5) * dx[0],
                                                          prob_lo[1] + (iv[1] + 0.5) * dx[1],
                                                          prob_lo[2] + (iv[2] + 0.5) * dx[2])};

      // Get the first species component BC
      const int* bc = bcr->data();

      // bcnormal handles all the state components at once
      amrex::Real s_ext[NVAR] = {0.0};

      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
         // Low
         if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

            // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
            bcnormal(x, m_nAux, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

            rho(iv) = s_ext[DENSITY];

         // High
         } else if ((bc[idir+AMREX_SPACEDIM] == amrex::BCType::ext_dir) and (iv[idir] > domhi[idir])) {

            // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
            bcnormal(x, m_nAux, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

            rho(iv) = s_ext[DENSITY];

         }
      }
   }
};

struct PeleLMCCFillExtDirRhoH
{

   ProbParm const* lprobparm;
   PmfData const* lpmfdata;
   const int m_nAux;

   AMREX_GPU_HOST
   constexpr PeleLMCCFillExtDirRhoH(ProbParm const* a_prob_parm, PmfData const* a_pmf_data, int a_nAux)
                                   : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux) {}

   AMREX_GPU_DEVICE
   void operator()(
      const amrex::IntVect& iv,
      amrex::Array4<amrex::Real> const& rhoH,
      const int dcomp,
      const int /*numcomp*/,
      amrex::GeometryData const& geom,
      const amrex::Real time,
      const amrex::BCRec* bcr,
      const int /*bcomp*/,
      const int orig_comp) const
   {

      // Get geometry data
      const int* domlo = geom.Domain().loVect();
      const int* domhi = geom.Domain().hiVect();
      const amrex::Real* prob_lo = geom.ProbLo();
      const amrex::Real* dx = geom.CellSize();
      const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(prob_lo[0] + (iv[0] + 0.5) * dx[0],
                                                          prob_lo[1] + (iv[1] + 0.5) * dx[1],
                                                          prob_lo[2] + (iv[2] + 0.5) * dx[2])};

      // Get the first species component BC
      const int* bc = bcr->data();

      // bcnormal handles all the state components at once
      amrex::Real s_ext[NVAR] = {0.0};

      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
         // Low
         if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

            // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
            bcnormal(x, m_nAux, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

            rhoH(iv) = s_ext[RHOH];

         // High
         } else if ((bc[idir+AMREX_SPACEDIM] == amrex::BCType::ext_dir) and (iv[idir] > domhi[idir])) {

            // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
            bcnormal(x, m_nAux, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

            rhoH(iv) = s_ext[RHOH];

         }
      }
   }
};

struct PeleLMCCFillExtDirTemp
{

   ProbParm const* lprobparm;
   PmfData const* lpmfdata;
   const int m_nAux;

   AMREX_GPU_HOST
   constexpr PeleLMCCFillExtDirTemp(ProbParm const* a_prob_parm, PmfData const* a_pmf_data, int a_nAux)
                                   : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux) {}

   AMREX_GPU_DEVICE
   void operator()(
      const amrex::IntVect& iv,
      amrex::Array4<amrex::Real> const& temp,
      const int dcomp,
      const int /*numcomp*/,
      amrex::GeometryData const& geom,
      const amrex::Real time,
      const amrex::BCRec* bcr,
      const int /*bcomp*/,
      const int orig_comp) const
   {

      // Get geometry data
      const int* domlo = geom.Domain().loVect();
      const int* domhi = geom.Domain().hiVect();
      const amrex::Real* prob_lo = geom.ProbLo();
      const amrex::Real* dx = geom.CellSize();
      const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(prob_lo[0] + (iv[0] + 0.5) * dx[0],
                                                          prob_lo[1] + (iv[1] + 0.5) * dx[1],
                                                          prob_lo[2] + (iv[2] + 0.5) * dx[2])};

      // Get the first species component BC
      const int* bc = bcr->data();

      // bcnormal handles all the state components at once
      amrex::Real s_ext[NVAR] = {0.0};

      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
         // Low
         if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

            // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
            bcnormal(x, m_nAux, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

            temp(iv) = s_ext[TEMP];

         // High
         } else if ((bc[idir+AMREX_SPACEDIM] == amrex::BCType::ext_dir) and (iv[idir] > domhi[idir])) {

            // bcnormal() is defined in pelelm_prob.H in problem directory in /Exec
            bcnormal(x, m_nAux, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

            temp(iv) = s_ext[TEMP];

         }
      }
   }
};

struct PeleLMCCFillExtDirDummy
{

   ProbParm const* lprobparm;
   const int m_nAux;

   AMREX_GPU_HOST
   constexpr PeleLMCCFillExtDirDummy(ProbParm const* a_prob_parm, int a_nAux)
                                     : lprobparm(a_prob_parm), m_nAux(a_nAux) {}

   AMREX_GPU_DEVICE
   void operator()(
      const amrex::IntVect& /*iv*/,
      amrex::Array4<amrex::Real> const& /*dummy*/,
      const int /*dcomp*/,
      const int /*numcomp*/,
      amrex::GeometryData const& /*geom*/,
      const amrex::Real /*time*/,
      const amrex::BCRec* /*bcr*/,
      const int /*bcomp*/,
      const int /*orig_comp*/) const
   { }
};
