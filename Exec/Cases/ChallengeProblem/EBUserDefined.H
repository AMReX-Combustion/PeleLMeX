#ifndef _EBUSERDEFINED_H_
#define _EBUSERDEFINED_H_

using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
void
EBUserDefined(const Geometry& geom,
              const int required_coarsening_level,
              const int max_coarsening_level)
{
    // ParmParse your geometry parameters
    // --> fixed geometry of the challenge problem

    // Build geometry pieces using EB2::* methods
    EB2::SplineIF Piston;
    std::vector<amrex::RealVect> lnpts;
    amrex::RealVect p;

    Real scaleFact  = 0.25;      // Geometry scaled by a factor 1/4
    scaleFact *= 0.01;           // MKS conversion

    p = amrex::RealVect(D_DECL(49.0*0.1*scaleFact, 7.8583*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(D_DECL(36.193*0.1*scaleFact, 7.8583*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);
    lnpts.clear();

    p = amrex::RealVect(D_DECL(36.193*0.1*scaleFact, 7.8583*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(D_DECL(24.035*0.1*scaleFact, -7.8586*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);
    lnpts.clear();

    p = amrex::RealVect(D_DECL(24.035*0.1*scaleFact, -7.8586*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(D_DECL(20.0*0.1*scaleFact, -7.8586*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);
    lnpts.clear();

    p = amrex::RealVect(D_DECL(20.0*0.1*scaleFact, -7.8586*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(D_DECL(1.9934*0.1*scaleFact, 3.464*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);
    lnpts.clear();

    p = amrex::RealVect(D_DECL(1.9934*0.1*scaleFact, 3.464*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    p = amrex::RealVect(D_DECL(0.09061*0.1*scaleFact, 3.464*0.1*scaleFact, 0.0));
    lnpts.push_back(p);
    Piston.addLineElement(lnpts);

    EB2::CylinderIF cylinder(48.0*0.1*scaleFact, 75.0*0.1*scaleFact, 2, {0.0, 0.0, -20.0*0.1*scaleFact}, true);

    auto revolvePiston  = EB2::lathe(Piston);
    auto PistonCylinder = EB2::makeUnion(revolvePiston, cylinder);

    // Build your geometry shop using EB2::makeShop
    auto gshop = EB2::makeShop(PistonCylinder);

    // Build geom using EB2::Build
    EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
}
#endif
#endif
