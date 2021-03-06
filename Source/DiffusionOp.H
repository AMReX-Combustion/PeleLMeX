#ifndef DIFFUSION_OP_H_
#define DIFFUSION_OP_H_

#include <AMReX_MLMG.H>
#include <AMReX_BC_TYPES.H>

#ifdef AMREX_USE_EB
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLEBTensorOp.H>
#endif
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLTensorOp.H>

class PeleLM;

class DiffusionOp
{
   public:

   // Methods

   DiffusionOp (PeleLM* a_pelelm);

   void diffuse_scalar (amrex::Vector<amrex::MultiFab*> const& a_phi,
                        int phi_comp,
                        amrex::Vector<amrex::MultiFab const*> const& a_rhs,
                        int rhs_comp,
                        amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> const& a_fluxes,
                        int flux_comp,
                        amrex::Vector<amrex::MultiFab const*> const& a_acoeff,
                        amrex::Vector<amrex::MultiFab const*> const& a_density,
                        amrex::Vector<amrex::MultiFab const*> const& a_bcoeff,
                        int bcoeff_comp,
                        amrex::Vector<amrex::BCRec> a_bcrec,
                        int ncomp,
                        int isPoissonSolve,
                        amrex::Real dt);

   void computeDiffLap (amrex::Vector<amrex::MultiFab*> const& laps, int lap_comp,
                        amrex::Vector<amrex::MultiFab const*> const& phi, int phi_comp,
                        amrex::Vector<amrex::MultiFab const*> const& density,
                        amrex::Vector<amrex::MultiFab const*> const& bcoeff, int bcoeff_comp,
                        amrex::Vector<amrex::BCRec> a_bcrec,
                        int ncomp);

   void computeDiffFluxes (amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> const& a_fluxes,
                           int flux_comp,
                           amrex::Vector<amrex::MultiFab const*> const& a_phi,
                           int phi_comp,
                           amrex::Vector<amrex::MultiFab const*> const& a_density,
                           amrex::Vector<amrex::MultiFab const*> const& a_bcoeff,
                           int bcoeff_comp,
                           amrex::Vector<amrex::BCRec> a_bcrec,
                           int ncomp,
                           amrex::Real scale,
                           int do_avgDown);

   void computeGradient(const amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> &a_grad,
                        const amrex::Vector<amrex::MultiFab *> &a_laps,
                        const amrex::Vector<amrex::MultiFab const*> &a_phi,
                        const amrex::BCRec &a_bcrec,
                        int do_avgDown);

   void scaleExtensiveFluxes(amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> const& a_fluxes,
                             int flux_comp,
                             int ncomp,
                             amrex::Real fac);

   void avgDownFluxes(amrex::Vector<amrex::Array<amrex::MultiFab*,AMREX_SPACEDIM>> const& a_fluxes,
                      int flux_comp,
                      int ncomp);

   void readParameters ();

   // Data and parameters
   PeleLM* m_pelelm;

   // Scalar diffusion apply operator
   std::unique_ptr<amrex::MLABecLaplacian> m_scal_apply_op;

   // Scalar diffusion solve operator
   std::unique_ptr<amrex::MLABecLaplacian> m_scal_solve_op;

   // Poisson-like operator for gradients
   std::unique_ptr<amrex::MLABecLaplacian> m_gradient_op;

   int m_verbose = 0;

   // Options to control MLMG behavior
   int m_mg_verbose = 0;
   int m_mg_bottom_verbose = 0;
   int m_mg_max_iter = 100;
   int m_mg_bottom_maxiter = 100;
   int m_mg_max_fmg_iter = 0;
   int m_mg_max_coarsening_level = 100;
   int m_mg_maxorder = 4;

   int m_num_pre_smooth = 2;
   int m_num_post_smooth = 2;

   amrex::Real m_mg_rtol = 1.0e-11;
   amrex::Real m_mg_atol = 1.0e-14;

   std::string m_mg_bottom_solver = "bicgstab";

};

class DiffusionTensorOp
{
public:

   DiffusionTensorOp (PeleLM* a_pelelm);

   void compute_divtau (amrex::Vector<amrex::MultiFab*> const& a_divtau,
                        amrex::Vector<amrex::MultiFab const*> const& a_vel,
                        amrex::Vector<amrex::MultiFab const*> const& a_density,
                        amrex::Vector<amrex::MultiFab const*> const& a_beta,
                        const amrex::BCRec &a_bcrec,
                        amrex::Real scale);

   void diffuse_velocity(amrex::Vector<amrex::MultiFab*> const& a_vel,
                         amrex::Vector<amrex::MultiFab const*> const& a_density,
                         amrex::Vector<amrex::MultiFab const*> const& a_beta,
                         const amrex::BCRec &a_bcrec,
                         amrex::Real dt);

private:

   void readParameters ();

   // Data and parameters
   PeleLM* m_pelelm;

   // Apply operator
   std::unique_ptr<amrex::MLTensorOp> m_apply_op;
   std::unique_ptr<amrex::MLTensorOp> m_solve_op;

   int m_verbose = 0;

   // Options to control MLMG behavior
   int m_mg_verbose = 0;
   int m_mg_bottom_verbose = 0;
   int m_mg_max_iter = 100;
   int m_mg_bottom_maxiter = 100;
   int m_mg_max_fmg_iter = 0;
   int m_mg_max_coarsening_level = 100;
   int m_mg_maxorder = 2;

   int m_num_pre_smooth = 2;
   int m_num_post_smooth = 2;

   amrex::Real m_mg_rtol = 1.0e-11;
   amrex::Real m_mg_atol = 1.0e-14;

   std::string m_mg_bottom_solver = "bicgstab";
};

#endif
