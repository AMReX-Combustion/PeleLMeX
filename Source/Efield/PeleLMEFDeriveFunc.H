#ifndef PELELMEF_DERIVEFUNC_H_
#define PELELMEF_DERIVEFUNC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>

void pelelm_derchargedist(PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                          const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                          const amrex::Geometry& geomdata,
                          amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derefx (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                    const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                    const amrex::Geometry& geomdata,
                    amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derefy (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                    const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                    const amrex::Geometry& geomdata,
                    amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derefz (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                    const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                    const amrex::Geometry& geomdata,
                    amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);
 
void pelelm_derLorentzx (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derLorentzy (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derLorentzz (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);
#endif
