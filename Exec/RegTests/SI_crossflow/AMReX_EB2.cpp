
#include <AMReX_EB2_IF_AllRegular.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Torus.H>
#include <AMReX_EB2_IF_Spline.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2.H>
#include <AMReX_ParmParse.H>
#include <AMReX.H>
#include <algorithm>

namespace amrex { namespace EB2 {

AMREX_EXPORT Vector<std::unique_ptr<IndexSpace> > IndexSpace::m_instance;

AMREX_EXPORT int max_grid_size = 64;
AMREX_EXPORT bool extend_domain_face = true;

void Initialize ()
{
    ParmParse pp("eb2");
    pp.queryAdd("max_grid_size", max_grid_size);
    pp.queryAdd("extend_domain_face", extend_domain_face);

    amrex::ExecOnFinalize(Finalize);
}

void Finalize ()
{
    IndexSpace::clear();
}

bool ExtendDomainFace ()
{
    return extend_domain_face;
}

void
IndexSpace::push (IndexSpace* ispace)
{
    auto r = std::find_if(m_instance.begin(), m_instance.end(),
                          [=] (const std::unique_ptr<IndexSpace>& x) -> bool
                          { return x.get() == ispace; });
    if (r == m_instance.end()) {
        m_instance.emplace_back(ispace);
    } else if (r+1 != m_instance.end()) {
        std::rotate(r, r+1, m_instance.end());
    }
}

void
IndexSpace::erase (IndexSpace* ispace)
{
    auto r = std::find_if(m_instance.begin(), m_instance.end(),
                          [=] (const std::unique_ptr<IndexSpace>& x) -> bool
                          { return x.get() == ispace; });
    if (r != m_instance.end()) {
        m_instance.erase(r);
    }
}

const IndexSpace* TopIndexSpaceIfPresent() noexcept {
    if (IndexSpace::size() > 0) {
        return &IndexSpace::top();
    }
    return nullptr;
}

void
Build (const Geometry& geom, int required_coarsening_level,
       int max_coarsening_level, int ngrow, bool build_coarse_level_by_coarsening)
{
    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);

    if (geom_type == "all_regular")
    {
        EB2::AllRegularIF rif;
        EB2::GeometryShop<EB2::AllRegularIF> gshop(rif);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening);
    }
    else if (geom_type == "box")
    {
        RealArray lo;
        pp.get("box_lo", lo);

        RealArray hi;
        pp.get("box_hi", hi);

        bool has_fluid_inside;
        pp.get("box_has_fluid_inside", has_fluid_inside);

        EB2::BoxIF bf(lo, hi, has_fluid_inside);

        EB2::GeometryShop<EB2::BoxIF> gshop(bf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening);
    }
    else if (geom_type == "cylinder")
    {
        RealArray center;
        pp.get("cylinder_center", center);

        Real radius;
        pp.get("cylinder_radius", radius);

        Real height = -1.0;
        pp.queryAdd("cylinder_height", height);

        int direction;
        pp.get("cylinder_direction", direction);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(direction >=0 && direction < 3,
                                         "eb2.cylinder_direction is invalid");

        bool has_fluid_inside;
        pp.get("cylinder_has_fluid_inside", has_fluid_inside);

        EB2::CylinderIF cf(radius, height, direction, center, has_fluid_inside);

        EB2::GeometryShop<EB2::CylinderIF> gshop(cf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening);
    }
    else if (geom_type == "plane")
    {
        RealArray point;
        pp.get("plane_point", point);

        RealArray normal;
        pp.get("plane_normal", normal);

        EB2::PlaneIF pf(point, normal);

        EB2::GeometryShop<EB2::PlaneIF> gshop(pf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening);
    }
    else if (geom_type == "sphere")
    {
        RealArray center;
        pp.get("sphere_center", center);

        Real radius;
        pp.get("sphere_radius", radius);

        bool has_fluid_inside;
        pp.get("sphere_has_fluid_inside", has_fluid_inside);

        EB2::SphereIF sf(radius, center, has_fluid_inside);

        EB2::GeometryShop<EB2::SphereIF> gshop(sf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening);
    }
    else if (geom_type == "torus")
    {
        RealArray center;
        pp.get("torus_center", center);

        Real small_radius;
        pp.get("torus_small_radius", small_radius);

        Real large_radius;
        pp.get("torus_large_radius", large_radius);

        bool has_fluid_inside = true;

        EB2::TorusIF sf(large_radius, small_radius, center, has_fluid_inside);

        EB2::GeometryShop<EB2::TorusIF> gshop(sf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening);
    }
    else if (geom_type == "3D-ignition_vessel") {
        amrex::Print() << "** Initialising 3D ignition_vessel geometry... \n";
        // amrex::ParmParse pp2("ignition_vessel");

        //Real fwl; 
        // pp.get("far_wall_loc",fwl);

        const amrex::Real *problo,*probhi;
        problo = geom.ProbLo();
        probhi = geom.ProbHi();
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
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening);
    }

    else if (geom_type == "3D-ignition_vessel_modified") {
        amrex::Print() << "** Initialising 3D ignition_vessel MODIFIED geometry... \n";
        amrex::ParmParse pp("ignition_vessel");

        // Real fwl; 
        // pp.get("far_wall_loc",fwl);

        const amrex::Real *problo,*probhi;
        problo = geom.ProbLo();
        probhi = geom.ProbHi();
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
        amrex::EB2::PlaneIF insulator_tip({AMREX_D_DECL(0.,cell_id*dy,0.)},
                                      {AMREX_D_DECL(0.,1.,0.)});

        cell_id = 0.225e-2 / dx;
        amrex::EB2::PlaneIF insulator_vert_left({AMREX_D_DECL(-0.225e-2,0.,0.)},
                                      {AMREX_D_DECL(1.,0.,0)});

        amrex::EB2::PlaneIF insulator_vert_right({AMREX_D_DECL(0.225e-2,0.,0.)},
                                       {AMREX_D_DECL(-1.,0.,0.)});
        
        cell_id = 0.077e-2 / dx;
        amrex::EB2::PlaneIF anode_vert_left({AMREX_D_DECL(-cell_id*dx,0.,0.)},
                                      {AMREX_D_DECL(1.,0.,0.)});

        amrex::EB2::PlaneIF anode_vert_right({AMREX_D_DECL(cell_id*dx,0.,0.)},
                                       {AMREX_D_DECL(-1.,0.,0.)});


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

        // auto polys = amrex::EB2::makeUnion(cathode,anode,insulator);
        // auto polys = amrex::EB2::makeUnion(cathode,insulator);
        auto polys = amrex::EB2::makeUnion(cathode,anode,insulator);

        int dir = 0;
        // pp.get("rot_dir",dir);
        amrex::Real angle = 270; 
        // pp.get("rot_angle",angle);

        auto pr     = amrex::EB2::lathe(polys);
        auto pr_rot = amrex::EB2::rotate(amrex::EB2::lathe(polys), angle*3.1415/180., dir);
        auto shop   = amrex::EB2::makeShop(pr_rot);

        EB2::Build(shop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening);
    }   
    else
    {
        amrex::Abort("geom_type "+geom_type+ " not supported");
    }
}

namespace {
static int comp_max_crse_level (Box cdomain, const Box& domain)
{
    int ilev;
    for (ilev = 0; ilev < 30; ++ilev) {
        if (cdomain.contains(domain)) break;
        cdomain.refine(2);
    }
    if (cdomain != domain) ilev = -1;
    return ilev;
}
}

int
maxCoarseningLevel (const Geometry& geom)
{
    const Box& domain = amrex::enclosedCells(geom.Domain());
    const Box& cdomain = IndexSpace::top().coarsestDomain();
    return comp_max_crse_level(cdomain, domain);
}

int
maxCoarseningLevel (IndexSpace const* ebis, const Geometry& geom)
{
    const Box& domain = amrex::enclosedCells(geom.Domain());
    const Box& cdomain = ebis->coarsestDomain();
    return comp_max_crse_level(cdomain,domain);
}

}}
