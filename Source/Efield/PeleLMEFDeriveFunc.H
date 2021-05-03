#ifndef PELELMEF_DERIVEFUNC_H_
#define PELELMEF_DERIVEFUNC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>

void pelelm_derchargedist(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                          const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                          const amrex::Geometry& geomdata,
                          amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derefx (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                    const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                   const amrex::Geometry& geomdata,
                   amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derefy (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                    const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                   const amrex::Geometry& geomdata,
                   amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derefz (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                    const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                   const amrex::Geometry& geomdata,
                   amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derLorentzx (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derLorentzy (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derLorentzz (const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);
#endif
