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

    amrex::Print() << "** Initialising bluff-body geometry... \n";
    
    amrex::ParmParse pp("EB");
    amrex::Real in_diam = 0.0;
    pp.query("in_diam",in_diam);
    amrex::Real out_diam = 0.0;
    pp.query("out_diam",out_diam);
    amrex::Real plane = 0.0;
    pp.query("horizontal_plane" ,plane);
    amrex::Real bluff_diam = 0.0;
    pp.query("bluff_diam"       ,bluff_diam);
    amrex::Real bluff_height = 0.0;
    pp.query("bluff_height"     ,bluff_height);
    // The cell_id step is done to have the vertical and horizontal planes 
    // exactly at the cell face and avoid numerical problems with small volume cells

    EB2::PlaneIF out_vert_right({AMREX_D_DECL(out_diam/2.0,0.,0.)},
                                   {AMREX_D_DECL(1.,0.,0.)});

    EB2::PlaneIF in_vert_right({AMREX_D_DECL(in_diam/2.0,0.,0.)},
                                   {AMREX_D_DECL(-1.,0.,0.)});
    
    int cell_id = (plane-problo[1])/dy;
    EB2::PlaneIF horizontal_plane({AMREX_D_DECL(0., plane,0.)},
                                  {AMREX_D_DECL(0.,-1.,0.)});


    amrex::Array<amrex::Real,AMREX_SPACEDIM> point0_cat;
    amrex::Array<amrex::Real,AMREX_SPACEDIM> point1_cat;
    
    point1_cat[0] =   in_diam/2.0; //x-coordinate
    point1_cat[1] =   bluff_height; //y-coordinate

    point0_cat[0] =  bluff_diam/2.; //x-coordinate
    point0_cat[1] =   0.0; //y-coordinate

    amrex::Array<amrex::Real,AMREX_SPACEDIM> norm0;

    norm0[0] =  (point0_cat[0]-point1_cat[0]);
    norm0[1] =  (point0_cat[1]-point1_cat[1]); 
    norm0[2] = 0.0;

    amrex::Real norm = sqrt(norm0[0]*norm0[0]+norm0[1]*norm0[1]);
    norm0[0] = norm0[0]/norm;
    norm0[1] = norm0[1]/norm;
    
    amrex::Real loc_x;
    pp.query("loc_x", loc_x );
    amrex::Real loc_y;
    pp.query("loc_y", loc_y );    
    EB2::PlaneIF bluff_body({AMREX_D_DECL(loc_x,loc_y,0)},{AMREX_D_DECL(-norm0[0],-norm0[1],0)});

    auto outer_geom = EB2::makeIntersection(out_vert_right,horizontal_plane);
    auto inner_geom = EB2::makeIntersection(in_vert_right,horizontal_plane);
    auto bluff      = EB2::makeIntersection(horizontal_plane,bluff_body);
    // auto inner_geom2 = EB2::makeIntersection(inner_geom,in_vert_right);

    auto polys = EB2::makeUnion(inner_geom,bluff,outer_geom);
    // auto polys = EB2::makeUnion(inner_geom2);

    auto pr     = EB2::lathe(polys);
    auto shop   = EB2::makeShop(pr);

    // EB2::Build(shop, geom, max_coarsening_level, max_coarsening_level);
    EB2::Build(shop, geom, required_coarsening_level,max_coarsening_level);

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBState(const amrex::Real xEBface[AMREX_SPACEDIM],
           amrex::Real s_ext[NVAR],
           const amrex::Real /*time*/,
           amrex::GeometryData const& geomdata,
           ProbParm const& prob_parm)
{
  // Provide the state on the EB centroid provided
  // by xEBface
  // Not the entire has to be specified. For instance, if using isothermal
  // boundary, only s_ext[TEMP] will be used.
  // Aborting here (message will not show on GPUs)
  // Abort("Default setEBState called. Copy and edit this file for your needs");
  
  amrex::Real r = sqrt(pow((0.0-xEBface[0]),2)+pow((0.0-xEBface[1]),2));

  // printf("setEBState z = %e\n",xEBface[2]);
  // Set higher temperature only for the EB face in contact with the recirculation zone
  if (r < 19e-3/2.0 and xEBface[2] > -1.0e-4) {
     s_ext[TEMP] = prob_parm.T_wall;
  }
  else{
     s_ext[TEMP] = prob_parm.T_in;
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBType(const amrex::Real * /*xEBface[AMREX_SPACEDIM]*/,
          amrex::Real & EBflagType,
          amrex::GeometryData const& /*geomdata*/,
          ProbParm const& /*prob_parm*/)
{
  // Provide a factor applied on the EB face diffusion
  // coefficient: 0.0 will turn off the diffusion (homogeneous Neumann BC)
  // while 1.0 will let it unchanged. Intermediate values can be specified
  // for partially Dirichlet BCs
  EBflagType = 1.0;
  // Abort("Default setEBType called. Copy and edit this file for your needs");
}

#endif
#endif
