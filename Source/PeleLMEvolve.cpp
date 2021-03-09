#include <PeleLM.H>

using namespace amrex;

void PeleLM::Evolve() {
   BL_PROFILE_VAR("PeleLM::Evolve()", Evolve);

   bool do_not_evolve = ( (m_max_step == 0) ||
                          ((m_stop_time >= 0.) && (m_cur_time > m_stop_time)) );

   int plt_justDidIt = 0;
   
   while(!do_not_evolve) {

      plt_justDidIt = 0;

      if (m_verbose > 0) {
         amrex::Print() << "\n ====================   NEW TIME STEP   ==================== \n";
      }

      int is_init = 0;
      Advance(is_init);
      m_nstep++;
      m_cur_time += m_dt;

      // Check for plot file
      if (writePlotNow()) {
         WritePlotFile();
         plt_justDidIt = 1;
      }

      // Check for checkpoint file
      // TODO

      // Check for the end of the simulation
      do_not_evolve = ( (m_max_step >= 0 && m_nstep >= m_max_step) ||
                        (m_stop_time >= 0.0 && m_cur_time >= m_stop_time - 1.0e-12 * m_dt) );

      // Move new -> old
      if (!do_not_evolve) {
         for (int lev = 0; lev <= finest_level; ++lev) {
            std::swap(m_leveldata_old[lev],m_leveldata_new[lev]);
         }
      }

   }

   if (m_verbose > 0) {
      amrex::Print() << "\n >> Final simulation time: " << m_cur_time << "\n";
   }
   if ( m_plot_int > 0 && !plt_justDidIt ) {
      WritePlotFile();
   }
   
}

bool
PeleLM::writePlotNow()
{
   bool write_now = false;

   if ( m_plot_int > 0 && (m_nstep % m_plot_int == 0) ) {
      write_now = true;
   }

   return write_now;
}
