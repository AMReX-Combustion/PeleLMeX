#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("Vin", PeleLM::prob_parm->Vin);
  pp.query("Thigh", PeleLM::prob_parm->Thigh);
  pp.query("Tlow", PeleLM::prob_parm->Tlow);

  PeleLM::prob_parm->bathID = N2_ID;
  PeleLM::prob_parm->oxidID = O2_ID;
  // get the fuel name from the peleLM.fuel_name declaration
  amrex::ParmParse pp_pele("peleLM");
  std::string fuelName = "";
  pp_pele.get("fuel_name", fuelName);
#if defined(H2_ID)
  if (fuelName == "H2") {
    PeleLM::prob_parm->fuelID = H2_ID;
  } else
#endif
#if defined(CH4_ID)
    if (fuelName == "CH4") {
    PeleLM::prob_parm->fuelID = CH4_ID;
  } else
#endif
  {
    amrex::Abort("fuel_name not recognised!");
  }
}
