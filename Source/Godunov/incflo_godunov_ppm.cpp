#include <Godunov.H>
#include <incflo_godunov_ppm.H>

using namespace amrex;

void godunov::predict_ppm (Box const& bx,
                           AMREX_D_DECL(amrex::Array4<amrex::Real> const& Imx,
                                        amrex::Array4<amrex::Real> const& Imy,
                                        amrex::Array4<amrex::Real> const& Imz),
                           AMREX_D_DECL(amrex::Array4<amrex::Real> const& Ipx,
                                        amrex::Array4<amrex::Real> const& Ipy,
                                        amrex::Array4<amrex::Real> const& Ipz),
                           Array4<Real const> const& q,
                           Array4<Real const> const& vel,
                           Geometry& geom,
                           Real dt,
                           BCRec const* pbc)
{
    const Box& domain = geom.Domain();
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    const auto dx = geom.CellSizeArray();
    AMREX_D_TERM(Real l_dtdx = dt / dx[0];,
                 Real l_dtdy = dt / dx[1];,
                 Real l_dtdz = dt / dx[2];);

    amrex::ParallelFor(bx, AMREX_SPACEDIM, 
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        AMREX_D_TERM(Godunov_ppm_pred_x(i,j,k,n,l_dtdx,vel(i,j,k,0),q,Imx,Ipx,pbc[n],dlo.x,dhi.x);,
                     Godunov_ppm_pred_y(i,j,k,n,l_dtdy,vel(i,j,k,1),q,Imy,Ipy,pbc[n],dlo.y,dhi.y);,
                     Godunov_ppm_pred_z(i,j,k,n,l_dtdz,vel(i,j,k,2),q,Imz,Ipz,pbc[n],dlo.z,dhi.z););
    });
}
