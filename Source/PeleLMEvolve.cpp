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

      if ( (m_regrid_int > 0) && (m_nstep > 0) && (m_nstep%m_regrid_int == 0) ) {
         if (m_verbose > 0) amrex::Print() << " Regridding...\n";
         regrid(0, m_cur_time);
         resetMacProjector();
         resetCoveredMask();
      }

      int is_init = 0;
      Advance(is_init);
      m_nstep++;
      m_cur_time += m_dt;

      // Temporals
      if (doTemporalsNow()) {
         writeTemporals();
      }

      // Check for plot file
      if (writePlotNow()) {
         WritePlotFile();
         plt_justDidIt = 1;
      }

      if (writeCheckNow()) {
         WriteCheckPointFile();
         chk_justDidIt = 1;
      }

      // Check for the end of the simulation
      do_not_evolve = ( (m_max_step >= 0 && m_nstep >= m_max_step) ||
                        (m_stop_time >= 0.0 && m_cur_time >= m_stop_time - 1.0e-12 * m_dt) );

   }

   if (m_verbose > 0) {
      amrex::Print() << "\n >> Final simulation time: " << m_cur_time << "\n";
   }
   if ( m_plot_int > 0 && !plt_justDidIt ) {
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
   }
   // TODO : time controled

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
