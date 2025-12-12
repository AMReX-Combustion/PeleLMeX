#ifndef PELELM_UTILS_H
#define PELELM_UTILS_H

template <class T>
amrex::Gpu::DeviceVector<T>
convertToDeviceVector(const amrex::Vector<T>& v)
{
  static_assert(
    std::is_trivially_copyable_v<T>, "T must be trivially copyable");
  int ncomp = v.size();
  amrex::Gpu::DeviceVector<T> v_d(ncomp);
#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(v_d.data(), v.data(), sizeof(T) * ncomp);
#else
  std::memcpy(v_d.data(), v.data(), sizeof(T) * ncomp);
#endif
  return v_d;
}

#endif
