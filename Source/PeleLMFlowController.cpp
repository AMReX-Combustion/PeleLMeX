#include "PeleLM.H"
#include "pelelm_prob.H"

#include <sys/stat.h>

using namespace amrex;

void
PeleLM::initActiveControl()
{
    ParmParse pp("active_control");

    pp.query("on",m_ctrl_active);

    if (!m_ctrl_active) {
        return;
    }

    Print() << " Initialization of active control \n";

    pp.query("use_temp",m_ctrl_useTemp);
    pp.query("temperature",m_ctrl_temperature);
    pp.query("tau",m_ctrl_tauControl);
    pp.query("v",m_ctrl_verbose);
    pp.query("height",m_ctrl_h);
    pp.query("changeMax",m_ctrl_changeMax);
    pp.query("velMax",m_ctrl_velMax);
    pp.query("pseudo_gravity",m_ctrl_pseudoGravity);
    pp.query("flow_dir",m_ctrl_flameDir);
    pp.query("AC_history",m_ctrl_AChistory);
    pp.query("npoints_average",m_ctrl_NavgPts);
    pp.query("method",m_ctrl_method);

    // Active control checks
    if ( m_ctrl_useTemp && (m_ctrl_temperature <= 0.0) ) {
        amrex::Error("active_control.temperature MUST be set with active_control.use_temp = 1");
    }

    if ( m_ctrl_active && (m_ctrl_tauControl <= 0.0) ) {
        amrex::Error("active_control.tau MUST be set when using active_control");
    }

    if ( m_ctrl_h <= 0.0 ) {
        amrex::Error("active_control.height MUST be set when using active_control");
    }

    if ( m_ctrl_flameDir > 2 ) {
        amrex::Error("active_control.flame_direction MUST be 0, 1 or 2 for X, Y and Z resp.");
    }

    // Initialize flow controller
    if (hasFlowControllerData<ProbParm>::value == false) {
        Abort("ProbParm doesn't have a FCData FlowControllerData member variable");
    }

    // Get FlowControllerData from ProbParm if it exists
    FlowControllerData* fcdata_host{nullptr};
    fcdata_host = getFCDataPtr(*prob_parm, hasFlowControllerData<ProbParm>{});

    // Resize vector for temporal average
    m_ctrl_time_pts.resize(m_ctrl_NavgPts+1,-1.0);
    m_ctrl_velo_pts.resize(m_ctrl_NavgPts+1,-1.0);
    m_ctrl_cntl_pts.resize(m_ctrl_NavgPts+1,-1.0);

    // Extract data from BC: assumes flow comes in from lo side of ctrl_flameDir
    ProbParm const* lprobparm = prob_parm_d;
    pele::physics::PMF::PmfData::DataContainer const* lpmfdata = pmf_data.getDeviceData();

    Gpu::DeviceVector<Real> s_ext_v(NVAR);
    Real* s_ext_d = s_ext_v.data();
    Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(Geom(0).ProbLo(0),
                                           Geom(0).ProbLo(1),
                                           Geom(0).ProbLo(2))};
    x[m_ctrl_flameDir] -= 1.0;
    const int ctrl_flameDir_l = m_ctrl_flameDir;
    const amrex::Real time_l = -1.0;
    const auto geomdata = Geom(0).data();

    Box dumbx({AMREX_D_DECL(0,0,0)},{AMREX_D_DECL(0,0,0)});
    amrex::ParallelFor(dumbx, [x, nAux=m_nAux, s_ext_d, ctrl_flameDir_l, time_l, geomdata, lprobparm, lpmfdata]
    AMREX_GPU_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept
    {
       bcnormal(x, nAux, s_ext_d, ctrl_flameDir_l, 1, time_l, geomdata, *lprobparm, lpmfdata);
    });
    Vector<Real> s_ext(NVAR);
    Gpu::copy(Gpu::deviceToHost, s_ext_v.begin(), s_ext_v.end(), s_ext.begin());

    if ( !m_ctrl_useTemp ) {
        // Get the fuel rhoY
        if (fuelID < 0) {
            Abort("Using activeControl based on fuel mass requires peleLM.fuelName !");
        }

        // Compute some active control parameters
        Real area_tot = 1.0;
        for (int idim{0}; idim < AMREX_SPACEDIM; idim++) {
            if (idim != m_ctrl_flameDir) area_tot *= (Geom(0).ProbHi(idim) - Geom(0).ProbLo(idim));
        }

        m_ctrl_scale = area_tot * s_ext[FIRSTSPEC+fuelID];
        m_ctrl_cfix  = m_ctrl_h * m_ctrl_scale;
    } else {
        // If using temp, scale is 1.0
        m_ctrl_scale = 1.0;
        m_ctrl_cfix = m_ctrl_h;
    }

    // Extract initial ctrl_V_in from bc
    m_ctrl_V_in  = s_ext[m_ctrl_flameDir];
    m_ctrl_V_in_old = m_ctrl_V_in;

    // Fill host ProbParm FCData
    fcdata_host->ctrl_active = 1;
    fcdata_host->ctrl_V_in = m_ctrl_V_in;

    if ( m_ctrl_verbose ) {
        if ( m_ctrl_useTemp ) {
            Print() << " Active control based on temperature iso-level activated."
                    << " Maintaining the flame at " << m_ctrl_h << " in " << m_ctrl_flameDir << " direction. \n";
        } else {
            Print() << " Active control based on fuel mass activated."
                    << " Maintaining the flame at " << m_ctrl_h << " in " << m_ctrl_flameDir << " direction. \n";
        }
    }
}

void
PeleLM::activeControl(int is_restart)
{
    if (!m_ctrl_active) {
        return;
    }

    // -------------------------------------------
    // Get the current target state (either amount of fuel or T-iso position)
    // -------------------------------------------
    Real coft = 0.0;
    if ( !m_ctrl_useTemp ) {
        // Compute the integral of the fuel mass in the domain
        coft = MFSum(GetVecOfConstPtrs(getSpeciesVect(AmrNewTime)),fuelID);;

    } else {
        // Get the low T position
        coft = 1.e37;
        getActiveControlLowT(coft);

    }

    // If first time, set old target state
    if ( m_ctrl_coftOld < 0.0 ) m_ctrl_coftOld = coft;

    if ( m_nstep == 0 ) m_ctrl_nfilled = m_ctrl_NavgPts+1;

    // Manage AC history I/O
    if ( is_restart ) {
        loadActiveControlHistory();
    }

    m_ctrl_V_in_old = m_ctrl_V_in;
    if ( m_dt > 0.0 ) {
       m_ctrl_V_in += m_dt * m_ctrl_dV;
    }

    Real slocal = 0.5 * (m_ctrl_V_in_old + m_ctrl_V_in) - (coft - m_ctrl_coftOld) / ( m_dt * m_ctrl_scale );

    // -------------------------------------------
    // Get s_est averaged from previous N steps
    // -------------------------------------------
    // Update m_ctrl_* Vectors if not restarting
    if ( !is_restart ) {
       m_ctrl_time_pts.insert(m_ctrl_time_pts.begin(),m_cur_time);
       m_ctrl_velo_pts.insert(m_ctrl_velo_pts.begin(),m_ctrl_V_in);
       m_ctrl_cntl_pts.insert(m_ctrl_cntl_pts.begin(),coft);
       if ( m_ctrl_time_pts.size() > m_ctrl_NavgPts ) { // Pop_back only if it's filled
          m_ctrl_time_pts.pop_back();
          m_ctrl_velo_pts.pop_back();
          m_ctrl_cntl_pts.pop_back();
       }
    }

    if ( m_ctrl_nfilled <= 0 ) {
       Real velIntegral = 0.0;
       for (int n = 1; n <= m_ctrl_NavgPts; n++ ) { // Piecewise constant velocity over NavgPts last steps
          velIntegral += 0.5 * ( m_ctrl_velo_pts[n-1] + m_ctrl_velo_pts[n] )
                             * ( m_ctrl_time_pts[n-1] - m_ctrl_time_pts[n] );
       }
       m_ctrl_sest = ( velIntegral - ( m_ctrl_cntl_pts[0] - m_ctrl_cntl_pts[m_ctrl_NavgPts]) / m_ctrl_scale )
                   / ( m_ctrl_time_pts[0] - m_ctrl_time_pts[m_ctrl_NavgPts] );

    } else {
       m_ctrl_nfilled -= 1;
       if ( m_nstep != 0 ) {
          m_ctrl_sest = (1.0 - m_ctrl_corr ) * m_ctrl_sest + m_ctrl_corr * slocal;
       }
    }

    // -------------------------------------------
    // Compute Vnew
    // -------------------------------------------
    Real Vnew = 0.0;

    if (m_ctrl_method == 1) {         // linear
       Real vslope = 2.0 * ( (m_ctrl_cfix - coft)/(m_ctrl_scale * m_ctrl_tauControl) + m_ctrl_sest - m_ctrl_V_in ) / m_ctrl_tauControl;
       Vnew = m_ctrl_V_in + m_dt * vslope;
    } else if (m_ctrl_method == 2) { // Quadratic 1
       Real vslope = 3.0 * ( (m_ctrl_cfix - coft)/(m_ctrl_scale * m_ctrl_tauControl) + m_ctrl_sest - m_ctrl_V_in ) / m_ctrl_tauControl;
       Vnew = m_ctrl_V_in + ( m_dt - 0.5 * m_dt*m_dt / m_ctrl_tauControl ) * vslope;
    } else if (m_ctrl_method == 3) { // Quadratic 2
       Real rhs2 = 2.0 * ( (m_ctrl_cfix - coft)/(m_ctrl_scale * m_ctrl_tauControl) + m_ctrl_sest - m_ctrl_V_in ) / m_ctrl_tauControl;
       Real rhs1 = ( m_ctrl_sest - m_ctrl_V_in ) / m_ctrl_tauControl;
       Real vt_tay = 3.0 * rhs2 - 2.0 * rhs1;
       Real vtt_tay = 6.0 * ( rhs1 - rhs2 ) / m_ctrl_tauControl;
       Vnew = m_ctrl_V_in + m_dt * vt_tay + 0.5 * m_dt * m_dt * vtt_tay;
    }

    // Limit Vnew
    Real dVmax = m_ctrl_changeMax * 1.0;
    Real dVmin = m_ctrl_changeMax * std::max(1.0,m_ctrl_V_in);
    Vnew = std::max(Vnew,0.0);
    Vnew = std::min(std::max(Vnew,m_ctrl_V_in-dVmin),m_ctrl_V_in+dVmax);
    Vnew = std::min(Vnew,m_ctrl_velMax);

    if (!is_restart && m_nstep > 0) {
       m_ctrl_tBase = m_cur_time;

       // Compute
       m_ctrl_dV = ( Vnew - m_ctrl_V_in ) / m_dt;
    }

    m_ctrl_coftOld = coft;

    // Get FlowControllerData from ProbParm if it exists
    FlowControllerData* fcdata_h{nullptr};
    FlowControllerData* fcdata_d{nullptr};
    fcdata_h = getFCDataPtr(*prob_parm, hasFlowControllerData<ProbParm>{});
    fcdata_d = getFCDataPtr(*prob_parm_d, hasFlowControllerData<ProbParm>{});

    // Pass dV and ctrl_V_in to FCData
    fcdata_h->ctrl_V_in = m_ctrl_V_in;
    fcdata_h->ctrl_dV = m_ctrl_dV;
    fcdata_h->ctrl_tBase = m_ctrl_tBase;

    // Update Device version
    Gpu::copy(Gpu::hostToDevice, fcdata_h, fcdata_h+1, fcdata_d);

    if ( m_ctrl_verbose && !is_restart) {
       Print() << "\n------------------------AC CONTROL------------------------ \n";
       Print() << " |     Time: " << m_cur_time << " -     dt: " << m_dt << "\n";
       Print() << " | Position: " << coft/m_ctrl_scale << " - Target: " << m_ctrl_cfix/m_ctrl_scale << "\n";
       Print() << " |    V_new: " << Vnew << " -   V_in: " << m_ctrl_V_in << " - dV: " << m_ctrl_dV << "\n";
       Print() << "---------------------------------------------------------- \n";
    }

    // Append to (or create) AC history file
    if (!is_restart) {
       std::ofstream ACfile(m_ctrl_AChistory.c_str(), std::ofstream::out | std::ofstream::app );
       Print(ACfile).SetPrecision(15) << m_nstep << "  "
                                      << m_cur_time << "  "
                                      << m_ctrl_V_in << "  "
                                      << slocal << "  "
                                      << m_ctrl_dV << "  "
                                      << m_ctrl_sest << "  "
                                      << m_ctrl_coftOld << "\n";
    }
}

void
PeleLM::getActiveControlLowT(Real &a_coft)
{
    for (int lev = 0; lev <= finest_level; lev++) {

        // Get t^{n+1} data pointer
        auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
        const auto geomdata = Geom(lev).data();

        // local FC data
        int  AC_FlameDir   = m_ctrl_flameDir;
        Real AC_Tcross     = m_ctrl_temperature;

        Real lowT = 1.e37;

        if ( lev != finest_level ) {
            lowT = amrex::ReduceMin(ldata_p->state, *m_coveredMask[lev], 0, [=]
                   AMREX_GPU_HOST_DEVICE(Box const& bx, Array4<Real const> const& T_arr,
                                                        Array4<int  const> const& covered_arr) -> Real
                   {
                       using namespace amrex::literals;
                       const auto lo = amrex::lbound(bx);
                       const auto hi = amrex::ubound(bx);
                       Real tmp_pos = 1.e37_rt;
                       const Real* prob_lo = geomdata.ProbLo();
                       const Real* dx = geomdata.CellSize();
                       for         (int k = lo.z; k <= hi.z; ++k) {
                           for     (int j = lo.y; j <= hi.y; ++j) {
                               for (int i = lo.x; i <= hi.x; ++i) {
                                   Real lcl_pos = 1.e37_rt;
                                   if ( T_arr(i,j,k,TEMP) > AC_Tcross && covered_arr(i,j,k) > 0) {
                                       int idx[3] = {i,j,k};
                                       idx[AC_FlameDir] -= 1;
                                       if ( T_arr(idx[0],idx[1],idx[2],TEMP) < AC_Tcross ) {
                                           Real coor[3] = {0.0};
                                           AMREX_D_TERM(coor[0] = prob_lo[0] + (i+0.5)*dx[0];,
                                                        coor[1] = prob_lo[1] + (j+0.5)*dx[1];,
                                                        coor[2] = prob_lo[2] + (k+0.5)*dx[2];);
                                           Real slope = ((T_arr(i,j,k,TEMP) ) - T_arr(idx[0],idx[1],idx[2],TEMP))/dx[AC_FlameDir];
                                           lcl_pos = coor[AC_FlameDir] - dx[AC_FlameDir] + ( AC_Tcross - T_arr(idx[0],idx[1],idx[2],TEMP) ) / slope;
                                       }
                                   }
                                   tmp_pos = amrex::min(tmp_pos,lcl_pos);
                               }
                           }
                       }
                       return tmp_pos;
                   });
        } else {
            lowT = amrex::ReduceMin(ldata_p->state, 0, [=]
                   AMREX_GPU_HOST_DEVICE(Box const& bx, Array4<Real const> const& T_arr) -> Real
                   {
                       using namespace amrex::literals;
                       const auto lo = amrex::lbound(bx);
                       const auto hi = amrex::ubound(bx);
                       Real tmp_pos = 1.e37_rt;
                       const Real* prob_lo = geomdata.ProbLo();
                       const Real* dx = geomdata.CellSize();
                       for         (int k = lo.z; k <= hi.z; ++k) {
                           for     (int j = lo.y; j <= hi.y; ++j) {
                               for (int i = lo.x; i <= hi.x; ++i) {
                                   Real lcl_pos = 1.e37_rt;
                                   if ( T_arr(i,j,k,TEMP) > AC_Tcross) {
                                       int idx[3] = {i,j,k};
                                       idx[AC_FlameDir] -= 1;
                                       if ( T_arr(idx[0],idx[1],idx[2],TEMP) < AC_Tcross ) {
                                           Real coor[3] = {0.0};
                                           AMREX_D_TERM(coor[0] = prob_lo[0] + (i+0.5)*dx[0];,
                                                        coor[1] = prob_lo[1] + (j+0.5)*dx[1];,
                                                        coor[2] = prob_lo[2] + (k+0.5)*dx[2];);
                                           Real slope = ((T_arr(i,j,k,TEMP) ) - T_arr(idx[0],idx[1],idx[2],TEMP))/dx[AC_FlameDir];
                                           lcl_pos = coor[AC_FlameDir] - dx[AC_FlameDir] + ( AC_Tcross - T_arr(idx[0],idx[1],idx[2],TEMP) ) / slope;
                                       }
                                   }
                                   tmp_pos = amrex::min(tmp_pos,lcl_pos);
                               }
                           }
                       }
                       return tmp_pos;
                   });
        }
        a_coft = amrex::min(a_coft,lowT);
    }
    ParallelDescriptor::ReduceRealMin(a_coft);
}

void
PeleLM::loadActiveControlHistory()
{
    m_ctrl_nfilled = m_ctrl_NavgPts+1;
    struct stat buffer;
    bool have_history = (stat (m_ctrl_AChistory.c_str(), &buffer) == 0);
    if ( have_history ) {
       if ( m_ctrl_verbose ) {
           Print() << " Setting AC from history from " << m_ctrl_AChistory << "\n";
       }
       std::fstream ACfile (m_ctrl_AChistory.c_str(), std::fstream::in);
       if (ACfile.is_open()) {
          while ( ACfile.good() ) {
             int step_io;
             Real time_io, vel_io, slocal_io, dV_io, s_est_io, coft_old_io;
             ACfile >> step_io
                    >> time_io
                    >> vel_io
                    >> slocal_io
                    >> dV_io
                    >> s_est_io
                    >> coft_old_io;
             if ( ( ( m_nstep - step_io ) >= 0 ) &&
                  ( ( m_nstep - step_io ) <= m_ctrl_NavgPts ) ) {   // Fill the previous step data
                int pos = ( m_nstep - step_io );
                m_ctrl_nfilled -= 1;
                m_ctrl_time_pts[pos] = time_io;
                m_ctrl_velo_pts[pos] = vel_io;
                m_ctrl_cntl_pts[pos] = coft_old_io;
             }
             if ( step_io == m_nstep ) {                          // Found the right step
                m_ctrl_V_in = vel_io;
                m_ctrl_tBase = time_io;
                m_ctrl_dV = dV_io;
                m_ctrl_sest = s_est_io;
                m_ctrl_coftOld = coft_old_io;
             }
          }
          ACfile.close();
       }
       if ( m_ctrl_verbose ) {
          Print() << " AC history arrays: \n";
          for (long int n = 0; n < m_ctrl_time_pts.size(); n++) {
             Print() << "  ["<<n<<"] time: "<< m_ctrl_time_pts[n]
                                <<", velo: "<< m_ctrl_velo_pts[n]
                                <<", coft: "<< m_ctrl_cntl_pts[n] << "\n";
          }
       }
    } else {
       if ( m_ctrl_verbose ) {
           Print() << " AC history file " << m_ctrl_AChistory << " not found, restarting from scratch \n";
       }
    }
}
