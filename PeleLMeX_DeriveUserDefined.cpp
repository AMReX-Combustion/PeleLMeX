#include <PeleLMeX_DeriveFunc.H>
#include <PeleLMeX_Index.H>
#include <PelePhysics.H>
#include <mechanism.H>
#include <PeleLMeX.H>
#include <PeleLMeX_K.H>

//
// User-defined derived variables list
//
amrex::Vector<std::string>
pelelmex_setuserderives()
{
  return {"derUserDefine_null"}; // var_names;
}

//
// User-defined derived definition
//
void
pelelmex_deruserdef(
  PeleLM* /*a_pelelm*/,
  const amrex::Box& /*bx*/,
  amrex::FArrayBox& /*derfab*/,
  int /*dcomp*/,
  int /*ncomp*/,
  const amrex::FArrayBox& /*statefab*/,
  const amrex::FArrayBox& /*reactfab*/,
  const amrex::FArrayBox& /*pressfab*/,
  const amrex::Geometry& /*geom*/,
  amrex::Real /*time*/,
  const amrex::Vector<amrex::BCRec>& /*bcrec*/,
  int /*level*/)
{
  amrex::Abort(
    "Using derUserDefine derived requires providing a definition in local "
    "DeriveUserDefined.cpp");
}
