#ifndef _EBUSERDEFINED_H_
#define _EBUSERDEFINED_H_

using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
void
EBUserDefined(const Geometry& geom,
              const int required_coarsening_level,
              const int max_coarsening_level)
{
    // ParmParse your geometry parameters
    // --> fixed geometry of the challenge problem
    amrex::ParmParse pp("eb2");

    amrex::Real radius = 48.0;
    pp.query("cylinder_radius", radius);
    amrex::Real zoffset = 0.0;
    pp.query("head_offset", zoffset);
    int useFlatHead = 0;
    pp.query("head_isFlat", useFlatHead);

    // Build geometry pieces using EB2::* methods
    EB2::SplineIF Piston;
    std::vector<amrex::RealVect> lnpts;
    amrex::RealVect p;

    Real scaleFact  = 0.25;      // Geometry scaled by a factor 1/4
    scaleFact *= 0.01;           // MKS conversion

    if (!useFlatHead) {
       p = amrex::RealVect(D_DECL(radius*1.1*0.1*scaleFact, (7.8583+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       p = amrex::RealVect(D_DECL(36.193*0.1*scaleFact, (7.8583+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       Piston.addLineElement(lnpts);
       lnpts.clear();

       p = amrex::RealVect(D_DECL(36.193*0.1*scaleFact, (7.8583+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       p = amrex::RealVect(D_DECL(24.035*0.1*scaleFact, (-7.8586+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       Piston.addLineElement(lnpts);
       lnpts.clear();

       p = amrex::RealVect(D_DECL(24.035*0.1*scaleFact, (-7.8586+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       p = amrex::RealVect(D_DECL(20.0*0.1*scaleFact, (-7.8586+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       Piston.addLineElement(lnpts);
       lnpts.clear();

       p = amrex::RealVect(D_DECL(20.0*0.1*scaleFact, (-7.8586+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       p = amrex::RealVect(D_DECL(1.9934*0.1*scaleFact, (3.464+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       Piston.addLineElement(lnpts);
       lnpts.clear();

       p = amrex::RealVect(D_DECL(1.9934*0.1*scaleFact, (3.464+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       p = amrex::RealVect(D_DECL(0.09061*0.1*scaleFact, (3.464+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       Piston.addLineElement(lnpts);
    } else {
       p = amrex::RealVect(D_DECL(radius*1.1*0.1*scaleFact, (-7.8586+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       p = amrex::RealVect(D_DECL(0.0*0.1*scaleFact, (-7.8586+zoffset)*0.1*scaleFact, 0.0));
       lnpts.push_back(p);
       Piston.addLineElement(lnpts);
    }

    EB2::CylinderIF cylinder(radius*0.1*scaleFact, 100.0*0.1*scaleFact, 2, {0.0, 0.0, -20.0*0.1*scaleFact}, true);

    auto revolvePiston  = EB2::lathe(Piston);
    auto PistonCylinder = EB2::makeUnion(revolvePiston, cylinder);

    // Build your geometry shop using EB2::makeShop
    auto gshop = EB2::makeShop(PistonCylinder);

    // Build geom using EB2::Build
    EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
}
#endif
#endif
