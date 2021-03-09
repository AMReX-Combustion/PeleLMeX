#include <AMReX_Slopes_K.H>
#include <Godunov.H>

using namespace amrex;

namespace {
    std::pair<bool,bool> has_extdir_or_ho (BCRec const* bcrec, int ncomp, int dir)
    {
        std::pair<bool,bool> r{false,false};
        for (int n = 0; n < ncomp; ++n) {
            r.first = r.first
                 or (bcrec[n].lo(dir) == BCType::ext_dir)
                 or (bcrec[n].lo(dir) == BCType::hoextrap);
            r.second = r.second
                 or (bcrec[n].hi(dir) == BCType::ext_dir)
                 or (bcrec[n].hi(dir) == BCType::hoextrap);
        }
        return r;
    }
}

void godunov::predict_plm_x (Box const& xebox,
                             Array4<Real> const& Imx, Array4<Real> const& Ipx,
                             Array4<Real const> const& q,
                             Array4<Real const> const& vcc,
                             Geometry& geom,
                             Real dt,
                             Vector<BCRec> const& h_bcrec,
                             BCRec const* pbc)
{
    const Real dx = geom.CellSize(0);
    const Real dtdx = dt/dx;

    int ncomp = AMREX_SPACEDIM;

    const Box& domain_box = geom.Domain();
    const int domain_ilo = domain_box.smallEnd(0);
    const int domain_ihi = domain_box.bigEnd(0);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::x));
    bool has_extdir_or_ho_lo = extdir_lohi.first;
    bool has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_ilo >= xebox.smallEnd(0)-1) or
        (has_extdir_or_ho_hi and domain_ihi <= xebox.bigEnd(0)))
    {
        amrex::ParallelFor(xebox, ncomp, [q,vcc,domain_ilo,domain_ihi,Imx,Ipx,dtdx,pbc]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const auto& bc = pbc[n];
            bool extdir_or_ho_ilo = (bc.lo(0) == BCType::ext_dir) or
                                    (bc.lo(0) == BCType::hoextrap);
            bool extdir_or_ho_ihi = (bc.hi(0) == BCType::ext_dir) or
                                    (bc.hi(0) == BCType::hoextrap);

            int order = 4;

            Real upls = q(i  ,j,k,n) + 0.5 * (-1.0 - vcc(i  ,j,k,0) * dtdx) *
                amrex_calc_xslope_extdir(i  ,j,k,n,order,q,extdir_or_ho_ilo,extdir_or_ho_ihi,domain_ilo,domain_ihi);
            Real umns = q(i-1,j,k,n) + 0.5 * ( 1.0 - vcc(i-1,j,k,0) * dtdx) *
                amrex_calc_xslope_extdir(i-1,j,k,n,order,q,extdir_or_ho_ilo,extdir_or_ho_ihi,domain_ilo,domain_ihi);

            Ipx(i-1,j,k,n) = umns;
            Imx(i  ,j,k,n) = upls;
        });
    }
    else
    {
        amrex::ParallelFor(xebox, ncomp, [q,vcc,Ipx,Imx,dtdx]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            int order = 4;

            Real upls = q(i  ,j,k,n) + 0.5 * (-1.0 - vcc(i  ,j,k,0) * dtdx) *
                amrex_calc_xslope(i  ,j,k,n,order,q);
            Real umns = q(i-1,j,k,n) + 0.5 * ( 1.0 - vcc(i-1,j,k,0) * dtdx) *
                amrex_calc_xslope(i-1,j,k,n,order,q);

            Ipx(i-1,j,k,n) = umns;
            Imx(i  ,j,k,n) = upls;
        });
    }
}

void godunov::predict_plm_y (Box const& yebox,
                            Array4<Real> const& Imy, Array4<Real> const& Ipy,
                            Array4<Real const> const& q,
                            Array4<Real const> const& vcc,
                            Geometry& geom,
                            Real dt,
                            Vector<BCRec> const& h_bcrec,
                            BCRec const* pbc)
{
    int ncomp = AMREX_SPACEDIM;

    const Real dy = geom.CellSize(1);
    const Real dtdy = dt/dy;

    const Box& domain_box = geom.Domain();
    const int domain_jlo = domain_box.smallEnd(1);
    const int domain_jhi = domain_box.bigEnd(1);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp,  static_cast<int>(Direction::y));
    bool has_extdir_or_ho_lo = extdir_lohi.first;
    bool has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_jlo >= yebox.smallEnd(1)-1) or
        (has_extdir_or_ho_hi and domain_jhi <= yebox.bigEnd(1)))
    {
        amrex::ParallelFor(yebox, ncomp, [q,vcc,domain_jlo,domain_jhi,Imy,Ipy,dtdy,pbc]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const auto& bc = pbc[n];
            bool extdir_or_ho_jlo = (bc.lo(1) == BCType::ext_dir) or
                                    (bc.lo(1) == BCType::hoextrap);
            bool extdir_or_ho_jhi = (bc.hi(1) == BCType::ext_dir) or
                                    (bc.hi(1) == BCType::hoextrap);

            int order = 4;

            Real vpls = q(i,j  ,k,n) + 0.5 * (-1.0 - vcc(i,j  ,k,1) * dtdy) *
                amrex_calc_yslope_extdir(i,j  ,k,n,order,q,extdir_or_ho_jlo,extdir_or_ho_jhi,domain_jlo,domain_jhi);
            Real vmns = q(i,j-1,k,n) + 0.5 * ( 1.0 - vcc(i,j-1,k,1) * dtdy) *
                amrex_calc_yslope_extdir(i,j-1,k,n,order,q,extdir_or_ho_jlo,extdir_or_ho_jhi,domain_jlo,domain_jhi);

            Ipy(i,j-1,k,n) = vmns;
            Imy(i,j  ,k,n) = vpls;
        });
    }
    else
    {
        amrex::ParallelFor(yebox, ncomp, [q,vcc,Ipy,Imy,dtdy]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            int order = 4;

            Real vpls = q(i,j  ,k,n) + 0.5 * (-1.0 - vcc(i,j  ,k,1) * dtdy) *
                amrex_calc_yslope(i,j  ,k,n,order,q);
            Real vmns = q(i,j-1,k,n) + 0.5 * ( 1.0 - vcc(i,j-1,k,1) * dtdy) *
                amrex_calc_yslope(i,j-1,k,n,order,q);

            Ipy(i,j-1,k,n) = vmns;
            Imy(i,j  ,k,n) = vpls;
        });
    }
}

#if (AMREX_SPACEDIM == 3)
void godunov::predict_plm_z (Box const& zebox,
                             Array4<Real> const& Imz, Array4<Real> const& Ipz,
                             Array4<Real const> const& q,
                             Array4<Real const> const& vcc,
                             Geometry& geom,
                             Real dt,
                             Vector<BCRec> const& h_bcrec,
                             BCRec const* pbc)
 {
    int ncomp = AMREX_SPACEDIM;

    const Real dz = geom.CellSize(2);
    const Real dtdz = dt/dz;

    const Box& domain_box = geom.Domain();
    const int domain_klo = domain_box.smallEnd(2);
    const int domain_khi = domain_box.bigEnd(2);

    // At an ext_dir boundary, the boundary value is on the face, not cell center.
    auto extdir_lohi = has_extdir_or_ho(h_bcrec.data(), ncomp, static_cast<int>(Direction::z));
    bool has_extdir_or_ho_lo = extdir_lohi.first;
    bool has_extdir_or_ho_hi = extdir_lohi.second;

    if ((has_extdir_or_ho_lo and domain_klo >= zebox.smallEnd(2)-1) or
        (has_extdir_or_ho_hi and domain_khi <= zebox.bigEnd(2)))
    {
        amrex::ParallelFor(zebox, ncomp, [q,vcc,domain_klo,domain_khi,Ipz,Imz,dtdz,pbc]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const auto& bc = pbc[n];
            bool extdir_or_ho_klo = (bc.lo(2) == BCType::ext_dir) or
                                    (bc.lo(2) == BCType::hoextrap);
            bool extdir_or_ho_khi = (bc.hi(2) == BCType::ext_dir) or
                                    (bc.hi(2) == BCType::hoextrap);

            int order = 4;

            Real wpls = q(i,j,k  ,n) + 0.5 * (-1.0 - vcc(i,j,k  ,2) * dtdz) *
                amrex_calc_zslope_extdir(i,j,k,n,order,q,extdir_or_ho_klo,extdir_or_ho_khi,domain_klo,domain_khi);
            Real wmns = q(i,j,k-1,n) + 0.5 * ( 1.0 - vcc(i,j,k-1,2) * dtdz) *
                amrex_calc_zslope_extdir(i,j,k-1,n,order,q,extdir_or_ho_klo,extdir_or_ho_khi,domain_klo,domain_khi);


            Ipz(i,j,k-1,n) = wmns;
            Imz(i,j,k  ,n) = wpls;
        });
    }
    else
    {
        amrex::ParallelFor(zebox, ncomp, [q,vcc,Ipz,Imz,dtdz]
        AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            int order = 4;

            Real wpls = q(i,j,k  ,n) + 0.5 * (-1.0 - vcc(i,j,k  ,2) * dtdz) *
                amrex_calc_zslope(i,j,k  ,n,order,q);
            Real wmns = q(i,j,k-1,n) + 0.5 * ( 1.0 - vcc(i,j,k-1,2) * dtdz) *
                amrex_calc_zslope(i,j,k-1,n,order,q);

            Ipz(i,j,k-1,n) = wmns;
            Imz(i,j,k  ,n) = wpls;
        });
    }
}

#endif
