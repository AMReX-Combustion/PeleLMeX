#ifndef PELELM_UTILS_H
#define PELELM_UTILS_H

template <class T>
amrex::Gpu::DeviceVector<T>
convertToDeviceVector(amrex::Vector<T> v)
{
  int ncomp = v.size();
  amrex::Gpu::DeviceVector<T> v_d(ncomp);
#ifdef AMREX_USE_GPU
  amrex::Gpu::htod_memcpy(v_d.data(), v.data(), sizeof(T) * ncomp);
#else
  std::memcpy(v_d.data(), v.data(), sizeof(T) * ncomp);
#endif
  return v_d;
}

template <typename T = amrex::Real, unsigned int N = NUM_SPECIES>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::GpuArray<T, N>
array4_to_array(
  int i, int j, int k, const amrex::Array4<T>& array, int start_comp = 0)
{
  amrex::GpuArray<T, N> celldata;
  for (unsigned int n = 0; n < N; n++) {
    celldata[n] = array(i, j, k, start_comp + n);
  }
  return celldata;
}

#endif
