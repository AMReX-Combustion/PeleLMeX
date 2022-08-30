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

    const amrex::Real *problo = geom.ProbLo();
    const amrex::Real *probhi = geom.ProbHi();
    amrex::Real dy = geom.CellSize()[1]; //* pow(2.0,max_coarsening_level);
    amrex::Real dx = geom.CellSize()[0]; //* pow(2.0,max_coarsening_level);

    // The cell_id step is done to have the vertical and horizontal planes 
    // exactly at the cell face and avoid numerical problems with small volume cells
    int cell_id = 0.15875e-2 / dx;
    amrex::EB2::PlaneIF cathode_vert_left({AMREX_D_DECL(-cell_id*dx,0.,0.)},
                                  {AMREX_D_DECL(1.,0.,0.)});

    amrex::EB2::PlaneIF cathode_vert_right({AMREX_D_DECL(cell_id*dx,0.,0.)},
                                   {AMREX_D_DECL(-1.,0.,0.)});

    cell_id = 2.121e-2 / dy;
    amrex::EB2::PlaneIF insulator_tip({AMREX_D_DECL(0.,2.1e-2,0.)},
                                  {AMREX_D_DECL(0.,1.,0.)});

    amrex::EB2::PlaneIF insulator_vert_left({AMREX_D_DECL(-0.225e-2,0.,0.)},
                                  {AMREX_D_DECL(1.,0.,0)});

    amrex::EB2::PlaneIF insulator_vert_right({AMREX_D_DECL(0.225e-2,0.,0.)},
                                   {AMREX_D_DECL(-1.,0.,0.)});
    
    cell_id = 0.077e-2 / dx;
    amrex::EB2::PlaneIF anode_vert_left({AMREX_D_DECL(-cell_id*dx,0.,0.)},
                                  {AMREX_D_DECL(1.,0.,0.)});

    amrex::EB2::PlaneIF anode_vert_right({AMREX_D_DECL(cell_id*dx,0.,0.)},
                                   {AMREX_D_DECL(-1.,0.,0.)});


    // Real radius = 0.001;
    // // pp.get("sphere_radius", radius);

    // bool has_fluid_inside = true;
    // // pp.get("sphere_has_fluid_inside", has_fluid_inside);

    // RealArray center;
    // center[0] = cell_id*dx+radius;
    // cell_id = 2.121e-2 / dy;
    // center[1] = 2.121e-2-radius-1.e-5;
    // center[2] = 0.0;
    // EB2::SphereIF sf_right(radius, center, has_fluid_inside);

    // center[0] = -cell_id*dx-radius;
    // cell_id = 2.121e-2 / dy;
    // center[1] = 2.121e-2-radius+1.e-5;
    // center[2] = 0.0;
    // EB2::SphereIF sf_left(radius, center, has_fluid_inside);

    /* ---------- Cathode cone ---------- */
    amrex::Array<amrex::Real,AMREX_SPACEDIM> point0_cat;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> point1_cat;
    
    point1_cat[0] = -0.15875e-2; //x-coordinate
    point1_cat[1] =  1.8122e-2; //y-coordinate

    point0_cat[0] =  0.00000e-2; //x-coordinate
    point0_cat[1] =  1.87000e-2; //y-coordinate

    amrex::Array<amrex::Real,AMREX_SPACEDIM> norm0;

    norm0[0] =  (point0_cat[0]-point1_cat[0]);
    norm0[1] =  (point0_cat[1]-point1_cat[1]); 
    norm0[2] = 0.0;

    amrex::Real norm = sqrt(norm0[0]*norm0[0]+norm0[1]*norm0[1]);
    norm0[0] = norm0[0]/norm;
    norm0[1] = norm0[1]/norm;
    
    amrex::EB2::PlaneIF cone_cathode_angle_left({AMREX_D_DECL(point1_cat[0],1.46147e-2,0)},
                                   {AMREX_D_DECL(norm0[0],-norm0[1],0)});

    amrex::EB2::PlaneIF cone_cathode_angle_right({AMREX_D_DECL(-point1_cat[0],1.46147e-2,0)},
                                   {AMREX_D_DECL(-norm0[0],-norm0[1],0)});

    // Anode tip is cut based on grid spacing in y to avoid numerical problems
    //get the cathode tip plane based on grid spacing to avoid cutting through a cell
    cell_id = (point0_cat[1]-1.80e-2-0.001e-2) / dy; 

    amrex::EB2::PlaneIF cathode_tip({AMREX_D_DECL(0.,(1.8e-2 + cell_id*dy),0)},
                            {AMREX_D_DECL(0.,-1.,0)});

    // This is adding a shallower cone angle close to the cathode tip
    amrex::Array<amrex::Real,AMREX_SPACEDIM> p0_cat;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> p1_cat;
    
    // p1_cat[0] = -4.1955e-03; //x-coordinate 40 degrees
    // p1_cat[0] = -1.8198e-03; //x-coordinate 20 degrees
    p1_cat[0] = -2.88675e-05; //x-coordinate 30 degrees
    p1_cat[1] =  1.87000e-2-0.005e-2; //y-coordinate

    p0_cat[0] =  0.00000; //x-coordinate
    p0_cat[1] =  1.87000e-2; //y-coordinate

    amrex::Array<amrex::Real,AMREX_SPACEDIM> n0;

    n0[0] =  (p0_cat[0]-p1_cat[0]);
    n0[1] =  (p0_cat[1]-p1_cat[1]); 
    n0[2] = 0.0;

    amrex::Real n = sqrt(n0[0]*n0[0]+n0[1]*n0[1]);
    n0[0] = n0[0]/n;
    n0[1] = n0[1]/n;
    
    amrex::EB2::PlaneIF cathode_chamfer_left({AMREX_D_DECL(p1_cat[0],1.872e-2,0)},
                                   {AMREX_D_DECL(n0[0],-n0[1],0)});

    amrex::EB2::PlaneIF cathode_chamfer_right({AMREX_D_DECL(-p1_cat[0],1.872e-2,0)},
                                   {AMREX_D_DECL(-n0[0],-n0[1],0)});

    // Making intersections
    //auto intersection05  = amrex::EB2::makeIntersection(cathode_vert_left ,cone_cathode_angle_left);
    //auto intersection105 = amrex::EB2::makeIntersection(intersection05    ,cathode_chamfer_left);
    //auto intersection1   = amrex::EB2::makeIntersection(intersection105   ,cathode_tip);
    //auto intersection2   = amrex::EB2::makeIntersection(cathode_vert_right,cone_cathode_angle_right);
    //auto intersection205 = amrex::EB2::makeIntersection(intersection2     ,cathode_chamfer_right);

    auto intersection05 = amrex::EB2::makeIntersection(cathode_vert_left ,cone_cathode_angle_left);
    auto intersection1 = amrex::EB2::makeIntersection(intersection05 ,cathode_tip);
    auto intersection2 = amrex::EB2::makeIntersection(cathode_vert_right,cone_cathode_angle_right);

    auto cathode = amrex::EB2::makeIntersection(intersection1,intersection2);


    auto intersection4 = amrex::EB2::makeIntersection(insulator_vert_left ,insulator_tip);
    auto intersection5 = amrex::EB2::makeIntersection(insulator_vert_right,insulator_tip);
    auto insulator = amrex::EB2::makeIntersection(intersection4,intersection5);

    // amrex::EB2::PlaneIF insulator_tip2({AMREX_D_DECL(0.,2.121e-2-radius,0.)},
    //                               {AMREX_D_DECL(0.,1.,0.)});

    // auto intersection4 = amrex::EB2::makeIntersection(insulator_vert_left ,insulator_tip);
    // auto insul1 = amrex::EB2::makeIntersection(intersection4,sf_left);
    // auto insul2 = amrex::EB2::makeIntersection(insul1,insulator_tip2);
    // auto insul = amrex::EB2::makeIntersection(insul2,sf_right);
    // auto insulator = amrex::EB2::makeIntersection(insul,insulator_tip);



    /* ---------- Anode cone ---------- */
    amrex::Array<amrex::Real,AMREX_SPACEDIM> point0_anode;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> point1_anode;

    
    cell_id = 0.07700e-2 / dx; 
    point1_anode[0] = -cell_id*dx; //x-coordinate
    cell_id = 2.05960e-2 / dy; 
    point1_anode[1] =  cell_id*dy; //y-coordinate

    cell_id = 0.00000 / dx; 
    point0_anode[0] =  cell_id*dx; //x-coordinate
    cell_id = 1.97500e-2 / dy; 
    point0_anode[1] =  cell_id*dy; //y-coordinate

    amrex::Array<amrex::Real,AMREX_SPACEDIM> norm1;

    norm1[0] =  (point0_anode[0]-point1_anode[0]);
    norm1[1] =  (point0_anode[1]-point1_anode[1]); 
    norm1[2] = 0.0;

    norm = sqrt(norm1[0]*norm1[0]+norm1[1]*norm1[1]);
    norm1[0] = norm1[0]/norm;
    norm1[1] = norm1[1]/norm;

    // amrex::Print() << " Anode norm " << norm1[0] << norm1[1];
    // amrex::Abort();

    cell_id = 2.03793e-2 / dy; 
    amrex::EB2::PlaneIF cone_anode_angle_left({AMREX_D_DECL(point1_anode[0],cell_id*dy,0)},
                                   {AMREX_D_DECL(norm1[0],-norm1[1],0)});

    cell_id = 2.03793e-2 / dy; 
    amrex::EB2::PlaneIF cone_anode_angle_right({AMREX_D_DECL(-point1_anode[0],cell_id*dy,0)},
                                   {AMREX_D_DECL(-norm1[0],-norm1[1],0)});

    cell_id = (1.98e-2 - 5.0e-05) / dy;
    amrex::EB2::PlaneIF anode_tip({AMREX_D_DECL(0.,cell_id*dy,0)},
                            {AMREX_D_DECL(0.,1.,0)});

     //cell_id = 1.98 / dy;
     //amrex::EB2::PlaneIF anode_tip({AMREX_D_DECL(0.,cell_id*dy,0)},
     //                        {AMREX_D_DECL(0.,1.,0)});

    auto intersection6_5 = amrex::EB2::makeIntersection(anode_vert_left ,cone_anode_angle_left);
    // auto intersection6 = amrex::EB2::makeIntersection(intersection6_5,anode_tip);
    auto intersection7 = amrex::EB2::makeIntersection(anode_vert_right,cone_anode_angle_right);
    auto anode = amrex::EB2::makeIntersection(intersection6_5,intersection7);

    //  ----> Squared anode tip
    dy = geom.CellSize()[1];
    dx = geom.CellSize()[0];
    cell_id = 1;
    amrex::EB2::PlaneIF cone_anode_vert_left({AMREX_D_DECL(-cell_id*dx,0.,0)},
                            {AMREX_D_DECL(1.,0.,0)});

    amrex::EB2::PlaneIF cone_anode_vert_right({AMREX_D_DECL(cell_id*dx,0.,0)},
                            {AMREX_D_DECL(-1.,0.,0)});

    cell_id = (1.98e-2 - 1*dy) / dy;
    amrex::EB2::PlaneIF tip({AMREX_D_DECL(0.,cell_id*dy,0)},
                            {AMREX_D_DECL(0.,1.,0)});

    auto intersect1 = amrex::EB2::makeIntersection(cone_anode_vert_left,tip);
    auto intersect2 = amrex::EB2::makeIntersection(cone_anode_vert_right,intersect1);
    //  ----> End of squared anode tip

    auto polys = amrex::EB2::makeUnion(cathode,anode,insulator);
    // auto polys = amrex::EB2::makeUnion(cathode,insulator);
    // auto polys = amrex::EB2::makeUnion(cathode);

    int dir = 0;
    // pp.get("rot_dir",dir);
    amrex::Real angle = 270; 
    // pp.get("rot_angle",angle);

    auto pr     = amrex::EB2::lathe(polys);
    auto pr_rot = amrex::EB2::rotate(amrex::EB2::lathe(polys), angle*3.1415/180., dir);
    auto shop   = amrex::EB2::makeShop(pr_rot);

	EB2::Build(shop, geom, required_coarsening_level,
                   max_coarsening_level);
}
#endif
#endif
