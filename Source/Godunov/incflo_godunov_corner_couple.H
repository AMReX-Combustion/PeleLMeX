#ifndef GODUNOV_CORNER_COUPLE_H
#define GODUNOV_CORNER_COUPLE_H

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
void Godunov_corner_couple_yx (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dx)*(state(i+1,j-1,k,n)*mac(i+1,j-1,k)
                               - state(i  ,j-1,k,n)*mac(i  ,j-1,k))
            - dt/(3.e0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k)
            + dt/(3.e0*dx)*s(i,j-1,k,n)*(mac(i+1,j-1,k)-mac(i,j-1,k));

        hi1  = hi - dt/(3.e0*dx)*(state(i+1,j,k,n)*mac(i+1,j,k)
                                - state(i  ,j,k,n)*mac(i  ,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i+1,j,k) - mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dx)*(mac  (i+1,j-1,k)   +   mac(i,j-1,k))
                              * (state(i+1,j-1,k,n) - state(i,j-1,k,n));
        hi1 = hi - dt/(6.e0*dx)*(mac  (i+1,j  ,k)   +   mac(i,j  ,k))
                              * (state(i+1,j  ,k,n) - state(i,j  ,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_zx (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dx,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dx)*(state(i+1,j,k-1,n)*mac(i+1,j,k-1)
                               - state(i  ,j,k-1,n)*mac(i  ,j,k-1))
            - dt/(3.e0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1)
            + dt/(3.e0*dx)*s(i,j,k-1,n)*(mac(i+1,j,k-1)-mac(i,j,k-1));

        hi1 = hi - dt/(3.e0*dx)*(state(i+1,j,k,n)*mac(i+1,j,k)
                               - state(i  ,j,k,n)*mac(i  ,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dx)*s(i,j,k,n)*(mac(i+1,j,k) - mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dx)*(  mac(i+1,j,k-1)   +   mac(i,j,k-1))
                              * (state(i+1,j,k-1,n) - state(i,j,k-1,n));
        hi1 = hi - dt/(6.e0*dx)*(  mac(i+1,j,k  )   +   mac(i,j,k))
                              * (state(i+1,j,k  ,n) - state(i,j,k  ,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_xy (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dy,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dy)*(state(i-1,j+1,k,n)*mac(i-1,j+1,k)
                               - state(i-1,j  ,k,n)*mac(i-1,j  ,k))
            - dt/3.e0     *s(i-1,j,k,n)*divu_cc(i-1,j,k)
            + dt/(3.e0*dy)*s(i-1,j,k,n)*(mac(i-1,j+1,k)-mac(i-1,j,k));

        hi1 = hi - dt/(3.e0*dy)*(state(i,j+1,k,n)*mac(i,j+1,k)
                               - state(i,j  ,k,n)*mac(i,j  ,k))
            - dt/3.e0     *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dy)*s(i,j,k,n)*(mac(i,j+1,k)-mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dy)*(  mac(i-1,j+1,k)   +   mac(i-1,j,k))
                              * (state(i-1,j+1,k,n) - state(i-1,j,k,n));
        hi1 = hi - dt/(6.e0*dy)*(mac  (i  ,j+1,k)   +   mac(i  ,j,k))
                              * (state(i  ,j+1,k,n) - state(i  ,j,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_zy (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dy,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dy)*(state(i,j+1,k-1,n)*mac(i,j+1,k-1)
                               - state(i,j  ,k-1,n)*mac(i,j  ,k-1))
            - dt/(3.e0)   *s(i,j,k-1,n)*divu_cc(i,j,k-1)
            + dt/(3.e0*dy)*s(i,j,k-1,n)*(mac(i,j+1,k-1)-mac(i,j,k-1));

        hi1 = hi - dt/(3.e0*dy)*(state(i,j+1,k,n)*mac(i,j+1,k)
                               - state(i,j  ,k,n)*mac(i,j  ,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dy)*s(i,j,k,n)*(mac(i,j+1,k) - mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dy)*(  mac(i,j+1,k-1)   +   mac(i,j,k-1))
                              * (state(i,j+1,k-1,n) - state(i,j,k-1,n));
        hi1 = hi - dt/(6.e0*dy)*(  mac(i,j+1,k)     +   mac(i,j,k))
                              * (state(i,j+1,k  ,n) - state(i,j,k  ,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_xz (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dz,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dz)*(state(i-1,j,k+1,n)*mac(i-1,j,k+1)
                               - state(i-1,j,k  ,n)*mac(i-1,j,k))
            - dt/(3.e0)   *s(i-1,j,k,n)*divu_cc(i-1,j,k)
            + dt/(3.e0*dz)*s(i-1,j,k,n)*(mac(i-1,j,k+1) - mac(i-1,j,k));

        hi1 = hi - dt/(3.e0*dz)*(state(i,j,k+1,n)*mac(i,j,k+1)
                               - state(i,j,k  ,n)*mac(i,j,k))
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k)
            + dt/(3.e0*dz)*s(i,j,k,n)*(mac(i,j,k+1) - mac(i,j,k));
    } else {
        lo1 = lo - dt/(6.e0*dz)*(  mac(i-1,j,k+1)   +   mac(i-1,j,k))
                              * (state(i-1,j,k+1,n) - state(i-1,j,k,n));
        hi1 = hi - dt/(6.e0*dz)*(  mac(i  ,j,k+1)   +   mac(i  ,j,k))
                              * (state(i  ,j,k+1,n) - state(i  ,j,k,n));
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void Godunov_corner_couple_yz (amrex::Real& lo1, amrex::Real& hi1,
                               int i, int j, int k, int n, amrex::Real dt, amrex::Real dz,
                               bool iconserv, amrex::Real lo, amrex::Real hi,
                               amrex::Array4<amrex::Real const> const& s,
                               amrex::Array4<amrex::Real const> const& divu_cc,
                               amrex::Array4<amrex::Real const> const& mac,
                               amrex::Array4<amrex::Real const> const& state)
{
    if (iconserv) {
        lo1 = lo - dt/(3.e0*dz)*(state(i,j-1,k+1,n)*mac(i,j-1,k+1)
                               - state(i,j-1,k  ,n)*mac(i,j-1,k))
            - dt/(3.e0)   *s(i,j-1,k,n)*divu_cc(i,j-1,k)
            + dt/(3.e0*dz)*s(i,j-1,k,n)*(mac(i,j-1,k+1)-mac(i,j-1,k));

        hi1 = hi - dt/(3.e0*dz)*(state(i,j,k+1,n)*mac(i,j,k+1)
                               - state(i,j,k  ,n)*mac(i,j,k)) 
            - dt/(3.e0)   *s(i,j,k,n)*divu_cc(i,j,k) 
            + dt/(3.e0*dz)*s(i,j,k,n)*(mac(i,j,k+1) - mac(i,j,k));

    } else {
        lo1 = lo - dt/(6.e0*dz)*(  mac(i,j-1,k+1)   +   mac(i,j-1,k))
                              * (state(i,j-1,k+1,n) - state(i,j-1,k,n));
        hi1 = hi - dt/(6.e0*dz)*(  mac(i,j  ,k+1)   +   mac(i,j  ,k))
                              * (state(i,j  ,k+1,n) - state(i,j  ,k,n));
    }
}

}
#endif
