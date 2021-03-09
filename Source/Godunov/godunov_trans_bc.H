#ifndef GODUNOV_TRANS_BC_H
#define GODUNOV_TRANS_BC_H

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_BCRec.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_Array.H>

/* This header file contains the inlined __host__ __device__ functions required for
   the scalar advection routines for 3D Godunov. It also contains function declarations
   for controlling host functions.  */

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_xbc(const int i, const int j, const int k, const int n,
                              const amrex::Array4<const amrex::Real> &s,
                                    amrex::Real &lo,
                                    amrex::Real &hi,
                              const int bclo, const int bchi,
                              const int domlo, const int domhi, bool is_velocity)
{
    using namespace amrex;


    // Low X
    if (i <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            lo = s(domlo-1,j,k,n);
            if( n == 0 and is_velocity ) hi=lo;
        }
        else if(bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
            lo = hi;

        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }

    // High X
    if (i > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            hi = s(domhi+1,j,k,n) ;
            if( n == 0 and is_velocity ) lo=hi;
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
            hi = lo;

        else if(bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_ybc(const int i, const int j, const int k, const int n,
                              const amrex::Array4<const amrex::Real> &s,
                                    amrex::Real &lo,
                                    amrex::Real &hi,
                              const int bclo, const int bchi,
                              const int domlo, const int domhi, bool is_velocity )
{
    using namespace amrex;


    // Low Y
    if (j <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            lo = s(i,domlo-1,k,n);
            if( n == 1 and is_velocity ) hi = lo;
        }

        else if (bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
            lo = hi;

        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }

    // High Y
    if (j > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            hi = s(i,domhi+1,k,n);
            if( n == 1 and is_velocity ) lo = hi ;
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
            hi = lo;

        else if (bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_trans_zbc(const int i, const int j, const int k, const int n,
                              const amrex::Array4<const amrex::Real> &s,
                                    amrex::Real &lo,
                                    amrex::Real &hi,
                              const int bclo, const int bchi,
                              const int domlo, const int domhi, bool is_velocity)
{
    using namespace amrex;


    // Low Z
    if (k <= domlo)
    {
        if (bclo==BCType::ext_dir)
        {
            lo =s(i,j,domlo-1,n);
            if( n == 2 and is_velocity ) hi = lo;
        }
        else if (bclo == BCType::foextrap || bclo == BCType::hoextrap || bclo == BCType::reflect_even)
            lo = hi;

        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }

    // High Z
    if (k > domhi)
    {
        if (bchi==BCType::ext_dir)
        {
            hi = s(i,j,domhi+1,n);
            if( n == 2 and is_velocity ) lo = hi ;
        }
        else if (bchi == BCType::foextrap || bchi == BCType::hoextrap || bchi == BCType::reflect_even)
            hi = lo;

        else if (bchi == BCType::reflect_odd)
        {
            hi = 0.;
            lo = 0.;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_xbc_lo(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s,
                       amrex::Real &lo, amrex::Real &hi,
                       const int bclo, const int domlo, bool is_velocity)
{
    using namespace amrex;

    if (i == domlo)
    {
        if (bclo == BCType::ext_dir )
        {
            lo = s(domlo-1, j, k, n);
            if ( n==0 and is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even )
            lo = hi;

        else if (bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_xbc_hi(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s,
                       amrex::Real &lo, amrex::Real &hi,
                       const int bchi, const int domhi, bool is_velocity )
{
    using namespace amrex;

    if (i == domhi+1)
    {
        if (bchi == BCType::ext_dir)
        {
            hi = s(domhi+1, j, k, n);
            if ( n==0 and is_velocity ) lo = hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even)
            hi = lo;

        else if (bchi == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_ybc_lo(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s,
                       amrex::Real &lo, amrex::Real &hi,
                       const int bclo, const int domlo, bool is_velocity  )
{
    using namespace amrex;

    if (j == domlo)
    {
        if (bclo == BCType::ext_dir)
        {
            lo = s(i, domlo-1, k, n);
            if ( n==1 and is_velocity ) hi=lo;
        }
        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even)
            lo = hi;

        else if (bclo == BCType::reflect_odd){
            hi = 0.;
            lo = hi;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_ybc_hi(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s,
                       amrex::Real &lo, amrex::Real &hi,
                       const int bchi, const int domhi, bool is_velocity  )
{
    using namespace amrex;

    if (j == domhi+1)
    {
        if (bchi == BCType::ext_dir)
        {
            hi = s(i, domhi+1, k, n);
            if ( n==1 and is_velocity ) lo=hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even)
            hi = lo;

        else if(bchi == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_zbc_lo(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s,
                       amrex::Real &lo, amrex::Real &hi,
                       const int bclo, const int domlo, bool is_velocity )
{
    using namespace amrex;

    if (k == domlo)
    {
        if (bclo == BCType::ext_dir)
        {
            lo = s(i, j, domlo-1, n);
            if ( n==2 and is_velocity ) hi=lo;
        }

        else if ( bclo == BCType::foextrap || bclo == BCType::hoextrap ||
                  bclo == BCType::reflect_even)
            lo = hi;

        else if(bclo == BCType::reflect_odd)
        {
            hi = 0.;
            lo = hi;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_cc_zbc_hi(const int i, const int j, const int k, const int n,
                       const amrex::Array4<const amrex::Real> &s,
                       amrex::Real &lo, amrex::Real &hi,
                       const int bchi, const int domhi, bool is_velocity  )
{
    using namespace amrex;

    if (k == domhi+1)
    {
        if (bchi == BCType::ext_dir)
        {
            hi = s(i,j,domhi+1, n);
            if ( n==2 and is_velocity ) lo=hi;
        }
        else if ( bchi == BCType::foextrap || bchi == BCType::hoextrap ||
                  bchi == BCType::reflect_even)
            hi = lo;

        else if(bchi == BCType::reflect_odd)
        {
            lo = 0.;
            hi = lo;
        }
    }
}

}
#endif
