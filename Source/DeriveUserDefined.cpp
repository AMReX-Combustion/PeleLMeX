#include <PeleLMDeriveFunc.H>
#include <PeleLM_Index.H>
#include <PelePhysics.H>
#include <mechanism.H>
#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;

//
// User-defined derived variables list
//
Vector<std::string> pelelm_setuserderives()
{
  //Vector<std::string> var_names({"derUserDefine_null"});
  return {"derUserDefine_null"}; //var_names;
}

//
// User-defined derived definition
//
void pelelm_deruserdef (PeleLM* /*a_pelelm*/, const Box& /*bx*/, FArrayBox& /*derfab*/, int /*dcomp*/, int /*ncomp*/,
                        const FArrayBox& /*statefab*/, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                        const Geometry& /*geom*/, Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)
{
    Abort("Using derUserDefine derived requires providing a definition in local DeriveUserDefined.cpp");
}
