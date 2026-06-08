#include <AMReX_ParmParse.H>
#include "PeleLMeX_Ascent.H"

namespace pele {
PeleLMeXAscent::PeleLMeXAscent()
{
  {
    amrex::ParmParse pp("ascent");
    pp.query("plot_int", plot_int);
  }
}
} // namespace pele
