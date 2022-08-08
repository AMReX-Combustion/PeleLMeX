#ifndef PELELM_DERIVEFUNC_H_
#define PELELM_DERIVEFUNC_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>

// Forward declaration of PeleLM
class PeleLM;

void pelelm_dermassfrac (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_dermolefrac (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_dertemp (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                     const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                     const amrex::Geometry& geomdata,
                     amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derheatrelease (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                            const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                            const amrex::Geometry& geomdata,
                            amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derkineticenergy (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                              const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                              const amrex::Geometry& geomdata,
                              amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derenstrophy (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                          const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                          const amrex::Geometry& geomdata,
                          amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_deravgpress (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                         const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                         const amrex::Geometry& geomdata,
                         amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_dermgvort (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                       const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                       const amrex::Geometry& geomdata,
                       amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_dermixfrac (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                        const amrex::Geometry& geomdata,
                        amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derprogvar (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                        const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                        const amrex::Geometry& geomdata,
                        amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_dervisc (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                     const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                     const amrex::Geometry& geomdata,
                     amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derdiffc (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                      const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                      const amrex::Geometry& geomdata,
                      amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

void pelelm_derlambda (PeleLM* a_pelelm, const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int ncomp,
                       const amrex::FArrayBox& statefab, const amrex::FArrayBox& reactfab, const amrex::FArrayBox& pressfab,
                       const amrex::Geometry& geomdata,
                       amrex::Real time, const amrex::Vector<amrex::BCRec> &bcrec, int level);

#ifdef PELE_USE_EFIELD
#include <PeleLMEFDeriveFunc.H>
#endif

#endif
