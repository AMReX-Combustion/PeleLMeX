#include <PeleLMeX_DeriveFunc.H>
#include <PeleLMeX_Index.H>
#include <PelePhysics.H>
#include <mechanism.H>
#include <PeleLMeX.H>
#include <PeleLMeX_K.H>

using namespace amrex;

//
// User-defined derived variables list
//
Vector<std::string>
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
  const Box& /*bx*/,
  FArrayBox& /*derfab*/,
  int /*dcomp*/,
  int /*ncomp*/,
  const FArrayBox& /*statefab*/,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)
{
  Abort("Using derUserDefine derived requires providing a definition in local "
        "DeriveUserDefined.cpp");
}
