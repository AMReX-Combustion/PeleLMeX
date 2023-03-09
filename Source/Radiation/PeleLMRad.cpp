#include <PeleLM.H>
#include <PeleLMDerive.H>

#include <PeleLMRad.hpp>

void PeleLM::computeRadSource(const PeleLM::TimeStamp& a_timestamp,
			      const amrex::Real a_dt) {
  int const co2Indx = rad_model->readRadIndices().co2Indx;
  int const h2oIndx = rad_model->readRadIndices().h2oIndx;
  int const coIndx = rad_model->readRadIndices().coIndx;

  for (int lev = 0; lev <= finest_level; lev++) {
    auto ldata_p = PeleLM::getLevelDataPtr(lev, a_timestamp);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*(m_extSource[lev]), amrex::TilingIfNotGPU());
	 mfi.isValid(); ++mfi) {
      auto const& q_yin_co2 = ldata_p->state.const_array(mfi, FIRSTSPEC + co2Indx);
      auto const& q_yin_h2o = ldata_p->state.const_array(mfi, FIRSTSPEC + h2oIndx);
      auto const& q_yin_co = ldata_p->state.const_array(mfi, FIRSTSPEC + coIndx);
      auto const& q_Tin = ldata_p->state.const_array(mfi, TEMP);
      auto const& q_Pin = ldata_p->state.const_array(mfi, RHORT);
#ifdef PELELM_USE_SOOT
      auto const& q_fvin = ldata_p->state.const_array(mfi, FIRSTSOOT + 1);
      rad_model->updateSpecProp(mfi, q_yin_co2, q_yin_h2o,
				  q_yin_co, q_Tin, q_Pin,
				  q_fvin, lev);
#else
      rad_model->updateSpecProp(mfi, q_yin_co2, q_yin_h2o,
				  q_yin_co, q_Tin, q_Pin, lev);
#endif
    }
  }

  rad_model->evaluateRad();

  for (int lev = 0; lev <= finest_level; lev++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*(m_extSource[lev]), amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
      auto const& source_arr = m_extSource[lev]->array(mfi, RHOH);
      rad_model->calcRadSource(mfi, source_arr, lev);
    }
  }
}
