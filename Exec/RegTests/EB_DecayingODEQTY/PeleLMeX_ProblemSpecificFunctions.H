#ifndef PROBLEMSPECIFICFUNCTIONS_H
#define PROBLEMSPECIFICFUNCTIONS_H

using namespace amrex;

#include <PeleLMeX_Index.H>
#include <PeleLMeX.H>

void set_ode_names(amrex::Vector<std::string>& a_ode_names)
{
#if NUM_ODE > 0
    a_ode_names.resize(NUM_ODE);
    for (int n = 0; n < NUM_ODE; n++) {
      a_ode_names[n] = "MY_ODE_" + std::to_string(n);
    }
#endif
}

#endif