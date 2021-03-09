//
//  Godunov.H
//

#ifndef Godunov_H
#define Godunov_H

#include <AMReX_MultiFabUtil.H>
#include <AMReX_BCRec.H>

namespace godunov {

    constexpr amrex::Real small_vel = 1.e-8;

    void predict_godunov (amrex::Real time, 
                          AMREX_D_DECL(amrex::MultiFab& u_mac,
                                       amrex::MultiFab& v_mac,
                                       amrex::MultiFab& w_mac),
                          amrex::MultiFab const& vel, 
                          amrex::MultiFab const& vel_forces,
                          amrex::Vector<amrex::BCRec> const& h_bcrec,
                                        amrex::BCRec  const* d_bcrec,
                          amrex::Geometry& geom,
                          amrex::Real dt, bool use_ppm, bool use_forces_in_trans);

    void predict_godunov (amrex::Real time, 
                          AMREX_D_DECL(amrex::MultiFab& u_mac,
                                       amrex::MultiFab& v_mac,
                                       amrex::MultiFab& w_mac),
                          amrex::MultiFab const& vel, 
                          amrex::MultiFab const& vel_forces,
                          amrex::Vector<amrex::BCRec> const& h_bcrec,
                                        amrex::BCRec  const* d_bcrec,
                          amrex::Geometry& geom,
                          amrex::Real dt, bool use_ppm, bool use_forces_in_trans,
                          AMREX_D_DECL(amrex::MultiFab const& gmacphi_x,
                                       amrex::MultiFab const& gmacphi_y,
                                       amrex::MultiFab const& gmacphi_z),
                          bool use_mac_phi_in_godunov);

    void predict_plm_x (amrex::Box const& bx,
                        amrex::Array4<amrex::Real> const& Imx, amrex::Array4<amrex::Real> const& Ipx,
                        amrex::Array4<amrex::Real const> const& q,
                        amrex::Array4<amrex::Real  const> const& vcc,
                        amrex::Geometry& geom,
                        amrex::Real dt,
                        amrex::Vector<amrex::BCRec> const& h_bcrec,
                        amrex::BCRec const* d_bcrec);

    void predict_plm_y (amrex::Box const& bx,
                        amrex::Array4<amrex::Real> const& Imy, amrex::Array4<amrex::Real> const& Ipy,
                        amrex::Array4<amrex::Real const> const& q,
                        amrex::Array4<amrex::Real  const> const& vcc,
                        amrex::Geometry& geom,
                        amrex::Real dt,
                        amrex::Vector<amrex::BCRec> const& h_bcrec,
                        amrex::BCRec const* d_bcrec);
#if (AMREX_SPACEDIM == 3)
    void predict_plm_z (amrex::Box const& bx,
                        amrex::Array4<amrex::Real> const& Imz, amrex::Array4<amrex::Real> const& Ipz,
                        amrex::Array4<amrex::Real const> const& q,
                        amrex::Array4<amrex::Real  const> const& vcc,
                        amrex::Geometry& geom,
                        amrex::Real dt,
                        amrex::Vector<amrex::BCRec> const& h_bcrec,
                        amrex::BCRec const* d_bcrec);
#endif

    void predict_ppm (amrex::Box const& bx,
                      AMREX_D_DECL(amrex::Array4<amrex::Real> const& Imx,
                                   amrex::Array4<amrex::Real> const& Imy,
                                   amrex::Array4<amrex::Real> const& Imz),
                      AMREX_D_DECL(amrex::Array4<amrex::Real> const& Ipx,
                                   amrex::Array4<amrex::Real> const& Ipy,
                                   amrex::Array4<amrex::Real> const& Ipz),
                      amrex::Array4<amrex::Real const> const& q,
                      amrex::Array4<amrex::Real const> const& vel,
                      amrex::Geometry& geom,
                      amrex::Real dt,
                      amrex::BCRec const* d_bcrec);

    void make_trans_velocities (AMREX_D_DECL(amrex::Box const& xbx,
                                             amrex::Box const& ybx, 
                                             amrex::Box const& zbx),
                                AMREX_D_DECL(amrex::Array4<amrex::Real> const& u_ad,
                                             amrex::Array4<amrex::Real> const& v_ad,
                                             amrex::Array4<amrex::Real> const& w_ad),
                                AMREX_D_DECL(amrex::Array4<amrex::Real const> const& Imx,
                                             amrex::Array4<amrex::Real const> const& Imy,
                                             amrex::Array4<amrex::Real const> const& Imz),
                                AMREX_D_DECL(amrex::Array4<amrex::Real const> const& Ipx,
                                             amrex::Array4<amrex::Real const> const& Ipy,
                                             amrex::Array4<amrex::Real const> const& Ipz),
                                amrex::Array4<amrex::Real const> const& vel,
                                amrex::Array4<amrex::Real const> const& f,
                                const amrex::Box& domain,
                                amrex::Real dt,
                                amrex::BCRec const* d_bcrec,
                                bool use_forces_in_trans);

    void predict_godunov_on_box (amrex::Box const& bx, int ncomp,
                                 AMREX_D_DECL(amrex::Box const& xbx, 
                                              amrex::Box const& ybx, 
                                              amrex::Box const& zbx),
                                 AMREX_D_DECL(amrex::Array4<amrex::Real> const& qx,
                                              amrex::Array4<amrex::Real> const& qy,
                                              amrex::Array4<amrex::Real> const& qz),
                                 amrex::Array4<amrex::Real const> const& q,
                                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& u_ad,
                                              amrex::Array4<amrex::Real const> const& v_ad,
                                              amrex::Array4<amrex::Real const> const& w_ad),
                                 AMREX_D_DECL(amrex::Array4<amrex::Real> const& Imx,
                                              amrex::Array4<amrex::Real> const& Imy,
                                              amrex::Array4<amrex::Real> const& Imz),
                                 AMREX_D_DECL(amrex::Array4<amrex::Real> const& Ipx,
                                              amrex::Array4<amrex::Real> const& Ipy,
                                              amrex::Array4<amrex::Real> const& Ipz),
                                 amrex::Array4<amrex::Real const> const& f,
                                 const amrex::Box& domain,
                                 const amrex::Real* dx,
                                 amrex::Real dt,
                                 amrex::BCRec const* d_bcrec,
                                 bool use_forces_in_trans,
                                 amrex::Real* p);

    void predict_godunov_on_box (amrex::Box const& bx, int ncomp,
                                 AMREX_D_DECL(amrex::Box const& xbx, 
                                              amrex::Box const& ybx, 
                                              amrex::Box const& zbx),
                                 AMREX_D_DECL(amrex::Array4<amrex::Real> const& qx,
                                              amrex::Array4<amrex::Real> const& qy,
                                              amrex::Array4<amrex::Real> const& qz),
                                 amrex::Array4<amrex::Real const> const& q,
                                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& u_ad,
                                              amrex::Array4<amrex::Real const> const& v_ad,
                                              amrex::Array4<amrex::Real const> const& w_ad),
                                 AMREX_D_DECL(amrex::Array4<amrex::Real> const& Imx,
                                              amrex::Array4<amrex::Real> const& Imy,
                                              amrex::Array4<amrex::Real> const& Imz),
                                 AMREX_D_DECL(amrex::Array4<amrex::Real> const& Ipx,
                                              amrex::Array4<amrex::Real> const& Ipy,
                                              amrex::Array4<amrex::Real> const& Ipz),
                                 amrex::Array4<amrex::Real const> const& f,
                                 const amrex::Box& domain,
                                 const amrex::Real* dx,
                                 amrex::Real dt,
                                 amrex::BCRec const* d_bcrec,
                                 bool use_forces_in_trans,
                                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& gmacphi_x,
                                              amrex::Array4<amrex::Real const> const& gmacphi_y,
                                              amrex::Array4<amrex::Real const> const& gmacphi_z),
                                 bool use_mac_phi_in_godunov,
                                 amrex::Real* p);

    void compute_godunov_advection (amrex::Box const& bx, int ncomp,
                                    amrex::Array4<amrex::Real> const& dqdt,
                                    amrex::Array4<amrex::Real const> const& q,
                                    AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                                 amrex::Array4<amrex::Real const> const& vmac,
                                                 amrex::Array4<amrex::Real const> const& wmac),
                                    amrex::Array4<amrex::Real const> const& fq,
                                    amrex::Array4<amrex::Real const> const& divu,
                                    amrex::Real dt,
                                    amrex::BCRec const* d_bcrec, 
                                    int const* iconserv,
                                    amrex::Real* p, 
                                    bool use_ppm, 
                                    bool use_forces_in_trans,
                                    amrex::Geometry& geom,
                                    bool is_velocity = false);

    void compute_godunov_advection (amrex::Box const& bx, int ncomp,
                                    amrex::Array4<amrex::Real> const& dqdt,
                                    amrex::Array4<amrex::Real const> const& q,
                                    AMREX_D_DECL(amrex::Array4<amrex::Real const> const& umac,
                                                 amrex::Array4<amrex::Real const> const& vmac,
                                                 amrex::Array4<amrex::Real const> const& wmac),
                                    amrex::Array4<amrex::Real const> const& fq,
                                    amrex::Array4<amrex::Real const> const& divu,
                                    AMREX_D_DECL(amrex::Array4<amrex::Real> const& qx,
                                                 amrex::Array4<amrex::Real> const& qy,
                                                 amrex::Array4<amrex::Real> const& qz),
                                    bool knownEdgeState,
                                    amrex::Real dt,
                                    amrex::BCRec const* d_bcrec, 
                                    int const* iconserv,
                                    amrex::Real* p, 
                                    bool use_ppm, 
                                    bool use_forces_in_trans,
                                    amrex::Geometry& geom,
                                    bool is_velocity = false);

} // namespace godunov

#endif /* Godunov_H */
