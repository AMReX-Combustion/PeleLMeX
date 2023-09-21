#ifndef PELELM_DERIVEFUNC_H
#define PELELM_DERIVEFUNC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>

// Forward declaration of PeleLM
class PeleLM;

void pelelmex_dermassfrac(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_dermolefrac(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_dertemp(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derheatrelease(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derkineticenergy(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derenstrophy(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_deravgpress(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derrhomrhoy(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_dermgvel(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_dermgvort(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_dervort(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_dercoord(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derQcrit(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_dermixfrac(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derprogvar(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_dervisc(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derdiffc(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derlambda(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derdmap(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

amrex::Vector<std::string> pelelmex_setuserderives();

void pelelmex_deruserdef(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geom,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

#ifdef PELE_USE_EFIELD
#include <PeleLMEFDeriveFunc.H>
#endif

#endif
