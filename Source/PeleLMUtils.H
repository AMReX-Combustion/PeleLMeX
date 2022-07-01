#ifndef PELELM_UTILS_H
#define PELELM_UTILS_H
template <class T>
amrex::Gpu::DeviceVector<T>
convertToDeviceVector ( amrex::Vector<T> v)
{
   int ncomp = v.size();
   amrex::Gpu::DeviceVector<T> v_d(ncomp);
#ifdef AMREX_USE_GPU
   amrex::Gpu::htod_memcpy(v_d.data(), v.data(), sizeof(T)*ncomp);
#else
   std::memcpy(v_d.data(), v.data(), sizeof(T)*ncomp);
#endif
   return v_d;
}
#endif
