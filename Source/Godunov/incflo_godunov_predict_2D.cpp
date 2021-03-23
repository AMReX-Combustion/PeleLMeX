
#include <incflo_godunov_trans_bc.H>
#include <Godunov.H>

using namespace amrex;

void godunov::predict_godunov (Real /*time*/, 
                               MultiFab& u_mac, MultiFab& v_mac,
                               MultiFab const& vel, 
                               MultiFab const& vel_forces,
                               Vector<BCRec> const& h_bcrec,
                                      BCRec  const* d_bcrec,
                               Geometry& geom, Real l_dt, 
#ifdef AMREX_USE_EB
                               bool /*use_ppm*/,
#else
                               bool use_ppm,
#endif
                               bool use_forces_in_trans,
                               MultiFab const& gmacphi_x, MultiFab const& gmacphi_y,
                               bool use_mac_phi_in_godunov)
{
    Box const& domain = geom.Domain();
    const Real* dx    = geom.CellSize();

    const int ncomp = AMREX_SPACEDIM;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox scratch;
        for (MFIter mfi(vel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Box const& bxg1 = amrex::grow(bx,1);
            Box const& xbx = mfi.nodaltilebox(0);
            Box const& ybx = mfi.nodaltilebox(1);

            Array4<Real> const& a_umac = u_mac.array(mfi);
            Array4<Real> const& a_vmac = v_mac.array(mfi);

            Array4<Real const> const& gmacphi_x_arr = gmacphi_x.const_array(mfi);
            Array4<Real const> const& gmacphi_y_arr = gmacphi_y.const_array(mfi);

            Array4<Real const> const& a_vel = vel.const_array(mfi);
            Array4<Real const> const& a_f = vel_forces.const_array(mfi);

            scratch.resize(bxg1, ncomp*(4*AMREX_SPACEDIM)+AMREX_SPACEDIM);
            Real* p = scratch.dataPtr();

            Array4<Real> Imx = makeArray4(p,bxg1,ncomp);
            p +=         Imx.size();
            Array4<Real> Ipx = makeArray4(p,bxg1,ncomp);
            p +=         Ipx.size();
            Array4<Real> Imy = makeArray4(p,bxg1,ncomp);
            p +=         Imy.size();
            Array4<Real> Ipy = makeArray4(p,bxg1,ncomp);
            p +=         Ipy.size();
            Array4<Real> u_ad = makeArray4(p,Box(bx).grow(1,1).surroundingNodes(0),1);
            p +=         u_ad.size();
            Array4<Real> v_ad = makeArray4(p,Box(bx).grow(0,1).surroundingNodes(1),1);
            p +=         v_ad.size();

#ifndef AMREX_USE_EB
            if (use_ppm){
                godunov::predict_ppm (bxg1, Imx, Imy, Ipx, Ipy, a_vel, a_vel,
                                      geom, l_dt, d_bcrec);
            }
            else
#endif
            {
                godunov::predict_plm_x (Box(u_ad), Imx, Ipx, a_vel, a_vel,
                                        geom, l_dt, h_bcrec, d_bcrec);
                godunov::predict_plm_y (Box(v_ad), Imy, Ipy, a_vel, a_vel,
                                        geom, l_dt, h_bcrec, d_bcrec);
            }

            make_trans_velocities(Box(u_ad), Box(v_ad),
                                  u_ad, v_ad,
                                  Imx, Imy, Ipx, Ipy, a_vel, a_f, 
                                  domain, l_dt, d_bcrec, use_forces_in_trans);

            predict_godunov_on_box(bx, ncomp, xbx, ybx, 
                                   a_umac, a_vmac,
                                   a_vel, u_ad, v_ad,
                                   Imx, Imy, Ipx, Ipy, a_f, 
                                   domain, dx, l_dt, d_bcrec, 
                                   use_forces_in_trans, 
                                   gmacphi_x_arr, gmacphi_y_arr, 
                                   use_mac_phi_in_godunov, p);

            Gpu::streamSynchronize();  // otherwise we might be using too much memory
        }
    }
}

void godunov::predict_godunov (Real /*time*/, 
                               MultiFab& u_mac, MultiFab& v_mac,
                               MultiFab const& vel, 
                               MultiFab const& vel_forces,
                               Vector<BCRec> const& h_bcrec,
                                      BCRec  const* d_bcrec,
                               Geometry& geom, Real l_dt, 
#ifdef AMREX_USE_EB
                               bool /*use_ppm*/,
#else
                               bool use_ppm,
#endif
                               bool use_forces_in_trans)
{
    Box const& domain = geom.Domain();
    const Real* dx    = geom.CellSize();

    // Switch off gmacphi
    bool use_mac_phi_in_godunov = false;

    const int ncomp = AMREX_SPACEDIM;
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox scratch;
        for (MFIter mfi(vel,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Box const& bxg1 = amrex::grow(bx,1);
            Box const& xbx = mfi.nodaltilebox(0);
            Box const& ybx = mfi.nodaltilebox(1);

            Array4<Real> const& a_umac = u_mac.array(mfi);
            Array4<Real> const& a_vmac = v_mac.array(mfi);

            // Setup dummy gmacphi
            Array4<Real const> const& gmacphi_x_arr = u_mac.const_array(mfi);
            Array4<Real const> const& gmacphi_y_arr = v_mac.const_array(mfi);

            Array4<Real const> const& a_vel = vel.const_array(mfi);
            Array4<Real const> const& a_f = vel_forces.const_array(mfi);

            scratch.resize(bxg1, ncomp*(4*AMREX_SPACEDIM)+AMREX_SPACEDIM);
            Real* p = scratch.dataPtr();

            Array4<Real> Imx = makeArray4(p,bxg1,ncomp);
            p +=         Imx.size();
            Array4<Real> Ipx = makeArray4(p,bxg1,ncomp);
            p +=         Ipx.size();
            Array4<Real> Imy = makeArray4(p,bxg1,ncomp);
            p +=         Imy.size();
            Array4<Real> Ipy = makeArray4(p,bxg1,ncomp);
            p +=         Ipy.size();
            Array4<Real> u_ad = makeArray4(p,Box(bx).grow(1,1).surroundingNodes(0),1);
            p +=         u_ad.size();
            Array4<Real> v_ad = makeArray4(p,Box(bx).grow(0,1).surroundingNodes(1),1);
            p +=         v_ad.size();

#ifndef AMREX_USE_EB
            if (use_ppm){
                godunov::predict_ppm (bxg1, Imx, Imy, Ipx, Ipy, a_vel, a_vel,
                                      geom, l_dt, d_bcrec);
            }
            else
#endif
            {
                godunov::predict_plm_x (Box(u_ad), Imx, Ipx, a_vel, a_vel,
                                        geom, l_dt, h_bcrec, d_bcrec);
                godunov::predict_plm_y (Box(v_ad), Imy, Ipy, a_vel, a_vel,
                                        geom, l_dt, h_bcrec, d_bcrec);
            }

            make_trans_velocities(Box(u_ad), Box(v_ad),
                                  u_ad, v_ad,
                                  Imx, Imy, Ipx, Ipy, a_vel, a_f, 
                                  domain, l_dt, d_bcrec, use_forces_in_trans);

            predict_godunov_on_box(bx, ncomp, xbx, ybx, 
                                   a_umac, a_vmac,
                                   a_vel, u_ad, v_ad,
                                   Imx, Imy, Ipx, Ipy, a_f, 
                                   domain, dx, l_dt, d_bcrec, 
                                   use_forces_in_trans, 
                                   gmacphi_x_arr, gmacphi_y_arr, 
                                   use_mac_phi_in_godunov, p);

            Gpu::streamSynchronize();  // otherwise we might be using too much memory
        }
    }
}

void godunov::make_trans_velocities (Box const& xbx, Box const& ybx, 
                                     Array4<Real> const& u_ad,
                                     Array4<Real> const& v_ad,
                                     Array4<Real const> const& Imx,
                                     Array4<Real const> const& Imy,
                                     Array4<Real const> const& Ipx,
                                     Array4<Real const> const& Ipy,
                                     Array4<Real const> const& vel,
                                     Array4<Real const> const& f,
                                     const Box& domain,
                                     Real l_dt, 
                                     BCRec  const* pbc,
                                     bool l_use_forces_in_trans)
{
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);

    // BCRec const* pbc = get_velocity_bcrec_device_ptr();

    amrex::ParallelFor(xbx, ybx, //zbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about x-velocity on x-faces here
        constexpr int n = 0;

        Real lo = Ipx(i-1,j,k,n);
        Real hi = Imx(i  ,j,k,n);

        if (l_use_forces_in_trans) {
            lo += 0.5*l_dt*f(i-1,j,k,n);
            hi += 0.5*l_dt*f(i  ,j,k,n);
        }

        auto bc = pbc[n];
        Godunov_trans_xbc(i, j, k, n, vel, lo, hi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
        u_ad(i,j,k) = ltm ? 0. : st;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // We only care about y-velocity on y-faces here
        constexpr int n = 1;

        Real lo = Ipy(i,j-1,k,n);
        Real hi = Imy(i,j  ,k,n);

        if (l_use_forces_in_trans) {
            lo += 0.5*l_dt*f(i,j-1,k,n);
            hi += 0.5*l_dt*f(i,j  ,k,n);
        }

        auto bc = pbc[n];
        Godunov_trans_ybc(i, j, k, n, vel, lo, hi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);

        Real st = ( (lo+hi) >= 0.) ? lo : hi;
        bool ltm = ( (lo <= 0. && hi >= 0.) || (amrex::Math::abs(lo+hi) < small_vel) );
        v_ad(i,j,k) = ltm ? 0. : st;
    });
}

void godunov::predict_godunov_on_box (Box const& bx, int ncomp,
                                      Box const& xbx, Box const& ybx, 
                                      Array4<Real> const& qx,
                                      Array4<Real> const& qy,
                                      Array4<Real const> const& q,
                                      Array4<Real const> const& u_ad,
                                      Array4<Real const> const& v_ad,
                                      Array4<Real> const& Imx,
                                      Array4<Real> const& Imy,
                                      Array4<Real> const& Ipx,
                                      Array4<Real> const& Ipy,
                                      Array4<Real const> const& f,
                                      const Box& domain,
                                      const Real* dx_arr,
                                      Real l_dt,
                                      BCRec  const* pbc,
                                      bool l_use_forces_in_trans,
                                      Array4<Real const> const& gmacphi_x,
                                      Array4<Real const> const& gmacphi_y,
                                      bool l_use_mac_phi_in_godunov,
                                      Real* p)
{
    const Dim3 dlo = amrex::lbound(domain);
    const Dim3 dhi = amrex::ubound(domain);
    Real dx = dx_arr[0];
    Real dy = dx_arr[1];

    Box xebox = Box(bx).grow(1,1).surroundingNodes(0);
    Box yebox = Box(bx).grow(0,1).surroundingNodes(1);
    Array4<Real> xlo = makeArray4(p, xebox, ncomp);
    p += xlo.size();
    Array4<Real> xhi = makeArray4(p, xebox, ncomp);
    p += xhi.size();
    Array4<Real> ylo = makeArray4(p, yebox, ncomp);
    p += ylo.size();
    Array4<Real> yhi = makeArray4(p, yebox, ncomp);
    p += yhi.size();

    amrex::ParallelFor(
        xebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipx(i-1,j,k,n);
            Real hi = Imx(i  ,j,k,n);

            if (l_use_forces_in_trans) {
                lo += 0.5*l_dt*f(i-1,j,k,n);
                hi += 0.5*l_dt*f(i  ,j,k,n);
            }

            Real uad = u_ad(i,j,k);
            auto bc = pbc[n];

            Godunov_trans_xbc(i, j, k, n, q, lo, hi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);

            xlo(i,j,k,n) = lo;
            xhi(i,j,k,n) = hi;

            Real st = (uad >= 0.) ? lo : hi;
            Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
            Imx(i, j, k, n) = fu*st + (1.0 - fu) *0.5 * (hi + lo); // store xedge
        },
        yebox, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real lo = Ipy(i,j-1,k,n);
            Real hi = Imy(i,j  ,k,n);

            if (l_use_forces_in_trans) {
                lo += 0.5*l_dt*f(i,j-1,k,n);
                hi += 0.5*l_dt*f(i,j  ,k,n);
            }

            Real vad = v_ad(i,j,k);
            auto bc = pbc[n];

            Godunov_trans_ybc(i, j, k, n, q, lo, hi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);

            ylo(i,j,k,n) = lo;
            yhi(i,j,k,n) = hi;


            Real st = (vad >= 0.) ? lo : hi;
            Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
            Imy(i, j, k, n) = fu*st + (1.0 - fu)*0.5*(hi + lo); // store yedge
        });

    Array4<Real> divu = makeArray4(Ipx.dataPtr(), grow(bx,1), 1);
    amrex::ParallelFor(Box(divu), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        divu(i,j,k) = 0.0;
    });

    // We can reuse the space in Ipy and Ipz.

    //
    // X-Flux
    //
    Box const xbxtmp = Box(xbx).enclosedCells().grow(0,1);
    Array4<Real> yzlo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(xbxtmp,1), 1);
    // Add d/dy term to z-faces
    // Start with {zlo,zhi} --> {zylo, zyhi} and upwind using w_ad to {zylo}
    // Add d/dz to y-faces
    // Start with {ylo,yhi} --> {yzlo, yzhi} and upwind using v_ad to {yzlo}
    amrex::ParallelFor(Box(yzlo),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 0;
        const auto bc = pbc[n];
        Real l_yzlo, l_yzhi;

        l_yzlo = ylo(i,j,k,n);
        l_yzhi = yhi(i,j,k,n);
        Real vad = v_ad(i,j,k);
        Godunov_trans_ybc(i, j, k, n, q, l_yzlo, l_yzhi, bc.lo(1), bc.hi(1), dlo.y, dhi.y, true);

        Real st = (vad >= 0.) ? l_yzlo : l_yzhi;
        Real fu = (amrex::Math::abs(vad) < small_vel) ? 0.0 : 1.0;
        yzlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_yzhi + l_yzlo);
    });
    //
    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 0;
        auto bc = pbc[n];
        Real stl = xlo(i,j,k,n) - (0.25*l_dt/dy)*(v_ad(i-1,j+1,k  )+v_ad(i-1,j,k))*
                                                 (yzlo(i-1,j+1,k  )-yzlo(i-1,j,k));
        Real sth = xhi(i,j,k,n) - (0.25*l_dt/dy)*(v_ad(i  ,j+1,k  )+v_ad(i  ,j,k))*
                                                 (yzlo(i  ,j+1,k  )-yzlo(i  ,j,k));
        if (!l_use_forces_in_trans) {
            stl += 0.5 * l_dt * f(i-1,j,k,n);
            sth += 0.5 * l_dt * f(i  ,j,k,n);
        }

        Real gphi_x;

        if (l_use_mac_phi_in_godunov)
        {
            // Note that the getFluxes call returns (-1/rho G^Mac phi) so we use the negative here
            gphi_x = -gmacphi_x(i,j,k);
            stl -= 0.5 * l_dt * gphi_x;
            sth -= 0.5 * l_dt * gphi_x;
        }

        Godunov_cc_xbc_lo(i, j, k, n, q, stl, sth, bc.lo(0), dlo.x, true);
        Godunov_cc_xbc_hi(i, j, k, n, q, stl, sth, bc.hi(0), dhi.x, true);

        // Prevent backflow
        if ( (i==dlo.x) and (bc.lo(0) == BCType::foextrap || bc.lo(0) == BCType::hoextrap) )
        {
            sth = amrex::min(sth,0.);
            stl = sth;
        }
        if ( (i==dhi.x+1) and (bc.hi(0) == BCType::foextrap || bc.hi(0) == BCType::hoextrap) )
        {
             stl = amrex::max(stl,0.);
             sth = stl;
        }
        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qx(i,j,k) = ltm ? 0. : st;

        if (l_use_mac_phi_in_godunov)
            qx(i,j,k) += 0.5 * l_dt * gphi_x;
    });
    //
    // Y-Flux
    //
    Box const ybxtmp = Box(ybx).enclosedCells().grow(1,1);
    Array4<Real> xzlo = makeArray4(Ipy.dataPtr(), amrex::surroundingNodes(ybxtmp,0), 1);

    // Add d/dz to x-faces
    // Start with {xlo,xhi} --> {xzlo, xzhi} and upwind using u_ad to {xzlo}
    // Add d/dx term to z-faces
    // Start with {zlo,zhi} --> {zxlo, zxhi} and upwind using w_ad to {zxlo}
    amrex::ParallelFor(Box(xzlo),
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 1;
        const auto bc = pbc[n];
        Real l_xzlo, l_xzhi;

        l_xzlo = xlo(i,j,k,n);
        l_xzhi = xhi(i,j,k,n);

        Real uad = u_ad(i,j,k);
        Godunov_trans_xbc(i, j, k, n, q, l_xzlo, l_xzhi, bc.lo(0), bc.hi(0), dlo.x, dhi.x, true);

        Real st = (uad >= 0.) ? l_xzlo : l_xzhi;
        Real fu = (amrex::Math::abs(uad) < small_vel) ? 0.0 : 1.0;
        xzlo(i,j,k) = fu*st + (1.0 - fu) * 0.5 * (l_xzhi + l_xzlo);
    });
    //
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        constexpr int n = 1;
        auto bc = pbc[n];
        Real stl = ylo(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j-1,k  )+u_ad(i,j-1,k))*
                                                 (xzlo(i+1,j-1,k  )-xzlo(i,j-1,k));
        Real sth = yhi(i,j,k,n) - (0.25*l_dt/dx)*(u_ad(i+1,j  ,k  )+u_ad(i,j  ,k))*
                                                 (xzlo(i+1,j  ,k  )-xzlo(i,j  ,k));
        if (!l_use_forces_in_trans) {
           stl += 0.5 * l_dt * f(i,j-1,k,n);
           sth += 0.5 * l_dt * f(i,j  ,k,n);
        }

        Real gphi_y;

        if (l_use_mac_phi_in_godunov)
        {
            // Note that the getFluxes call returns (-1/rho G^Mac phi) so we use the negative here
            gphi_y = -gmacphi_y(i,j,k);
            stl -= 0.5 * l_dt * gphi_y;
            sth -= 0.5 * l_dt * gphi_y;
        }

        Godunov_cc_ybc_lo(i, j, k, n, q, stl, sth, bc.lo(1), dlo.y, true);
        Godunov_cc_ybc_hi(i, j, k, n, q, stl, sth, bc.hi(1), dhi.y, true);

        // Prevent backflow
        if ( (j==dlo.y) and (bc.lo(1) == BCType::foextrap || bc.lo(1) == BCType::hoextrap) )
        {
            sth = amrex::min(sth,0.);
            stl = sth;
        }
        if ( (j==dhi.y+1) and (bc.hi(1) == BCType::foextrap || bc.hi(1) == BCType::hoextrap) )
        {
            stl = amrex::max(stl,0.);
            sth = stl;
        }

        Real st = ( (stl+sth) >= 0.) ? stl : sth;
        bool ltm = ( (stl <= 0. && sth >= 0.) || (amrex::Math::abs(stl+sth) < small_vel) );
        qy(i,j,k) = ltm ? 0. : st;

        if (l_use_mac_phi_in_godunov)
            qy(i,j,k) += 0.5 * l_dt * gphi_y;
    });

}
