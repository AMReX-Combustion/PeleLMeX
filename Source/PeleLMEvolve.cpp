#include <PeleLM.H>

using namespace amrex;

void PeleLM::Evolve() {
   BL_PROFILE("PeleLM::Evolve()");

   bool do_not_evolve = ( (m_max_step == 0) ||
                          ((m_stop_time >= 0.) && (m_cur_time > m_stop_time)) );

   int plt_justDidIt = 0;
   int chk_justDidIt = 0;

   while(!do_not_evolve) {

      plt_justDidIt = 0;
      chk_justDidIt = 0;

      if (m_verbose > 0) {
         amrex::Print() << "\n ====================   NEW TIME STEP   ==================== \n";
      }

      bool regridded = false;
      if ( (m_regrid_int > 0) && (m_nstep > 0) && (m_nstep%m_regrid_int == 0) ) {
         if (m_verbose > 0) amrex::Print() << " Regridding...\n";
         regrid(0, m_cur_time);
         resetMacProjector();
         resetCoveredMask();
         regridded = true;
         updateDiagnostics();
      }
#ifdef PELELM_USE_SPRAY
      // Inject and redistribute spray particles
      if (do_spray_particles && regridded) {
        sprayPostRegrid();
      }
#endif
      int is_init = 0;
      Advance(is_init);
      m_nstep++;
      m_cur_time += m_dt;

#ifdef PELELM_USE_SPRAY
      if (do_spray_particles) {
        sprayInjectRedist();
      }
#endif

      // Temporals
      if (doTemporalsNow()) {
         writeTemporals();
      }

      // Diagnostics
      doDiagnostics();

      // Check message
      bool dump_and_stop = checkMessage();

      // Check for plot file
      if (writePlotNow() || dump_and_stop) {
         WritePlotFile();
         plt_justDidIt = 1;
      }

      if (writeCheckNow() || dump_and_stop) {
         WriteCheckPointFile();
         chk_justDidIt = 1;
      }

      // Check for the end of the simulation
      do_not_evolve = ( (m_max_step >= 0 && m_nstep >= m_max_step) ||
                        (m_stop_time >= 0.0 && m_cur_time >= m_stop_time - 1.0e-12 * m_dt) ||
                        dump_and_stop );

   }

   if (m_verbose > 0) {
      amrex::Print() << "\n >> Final simulation time: " << m_cur_time << "\n";
   }
   if ( (m_plot_int > 0 || m_plot_per_approx > 0. || m_plot_per_exact > 0.) && !plt_justDidIt ) {
      WritePlotFile();
   }
   if ( m_check_int > 0 && !chk_justDidIt ) {
      WriteCheckPointFile();
   }
}

bool
PeleLM::writePlotNow()
{
   bool write_now = false;

   if ( m_plot_int > 0 && (m_nstep % m_plot_int == 0) ) {
      write_now = true;

   } else if ( m_plot_per_exact > 0.0 && (std::abs(std::remainder(m_cur_time, m_plot_per_exact)) < 1.e-12) ) {
      write_now = true;

   } else if (m_plot_per_approx > 0.0) {
      // Check to see if we've crossed a plot_per interval by comparing
      // the number of intervals that have elapsed for both the current
      // time and the time at the beginning of this timestep.
      int num_per_old = static_cast<int>((m_cur_time-m_dt) / m_plot_per_approx);
      int num_per_new = static_cast<int>((m_cur_time     ) / m_plot_per_approx);
      // Before using these, however, we must test for the case where we're
      // within machine epsilon of the next interval. In that case, increment
      // the counter, because we have indeed reached the next plot_per interval
      // at this point.
      const Real eps = std::numeric_limits<Real>::epsilon() * 10.0_rt * std::abs(m_cur_time);
      const Real next_plot_time = (num_per_old + 1) * m_plot_per_approx;
      if ((num_per_new == num_per_old) && std::abs(m_cur_time - next_plot_time) <= eps)
      {
         num_per_new += 1;
      }
      // Similarly, we have to account for the case where the old time is within
      // machine epsilon of the beginning of this interval, so that we don't double
      // count that time threshold -- we already plotted at that time on the last timestep.
      if ((num_per_new != num_per_old) && std::abs((m_cur_time - m_dt) - next_plot_time) <= eps)
      {
         num_per_old += 1;
      }
      if (num_per_old != num_per_new)
      {
         write_now = true;
      }
   }

   return write_now;
}

bool
PeleLM::writeCheckNow()
{
   bool write_now = false;

   if ( m_check_int > 0 && (m_nstep % m_check_int == 0) ) {
      write_now = true;
   }
   // TODO : time controled ?

   return write_now;
}

bool
PeleLM::doTemporalsNow()
{
   bool write_now = false;

   if ( m_do_temporals &&
        (m_nstep % m_temp_int == 0) ) {
      write_now = true;
   }

   return write_now;
}

bool
PeleLM::checkMessage()
{
    bool dump_and_stop = false;
    if (m_nstep % m_message_int == 0) {
        int dumpclose = 0;
        if (ParallelDescriptor::IOProcessor()) {
            FILE *fp;
            if ( (fp=fopen("dump_and_stop","r")) != 0 ) {
                remove("dump_and_stop");
                dumpclose = 1;
                fclose(fp);
            }
        }
        int packed_data[1];
        packed_data[0] = dumpclose;
        ParallelDescriptor::Bcast(packed_data, 1, ParallelDescriptor::IOProcessorNumber());
        dump_and_stop = packed_data[0];
    }
    return dump_and_stop;
}
