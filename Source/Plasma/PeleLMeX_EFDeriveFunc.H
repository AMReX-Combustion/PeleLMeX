#ifndef PELELMEF_DERIVEFUNC_H
#define PELELMEF_DERIVEFUNC_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>

void pelelmex_derchargedist(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derefx(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derefy(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derefz(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derLorentzx(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derLorentzy(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);

void pelelmex_derLorentzz(
  PeleLM* a_pelelm,
  const amrex::Box& bx,
  amrex::FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const amrex::FArrayBox& statefab,
  const amrex::FArrayBox& reactfab,
  const amrex::FArrayBox& pressfab,
  const amrex::Geometry& geomdata,
  amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcrec,
  int level);
#endif
