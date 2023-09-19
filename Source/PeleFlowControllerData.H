#ifndef PELE_FLOWCONTROLDATA_H
#define PELE_FLOWCONTROLDATA_H
#include <type_traits>
#include <AMReX_REAL.H>

struct FlowControllerData
{
  amrex::Real ctrl_dV{0.0};
  amrex::Real ctrl_V_in{0.0};
  amrex::Real ctrl_tBase{0.0};
  unsigned int ctrl_active{0};
};

// SFINAE to detect if a FCData is in the ProbParm or not
template <typename T, typename = void>
struct hasFlowControllerData : std::false_type
{
};

template <typename T>
struct hasFlowControllerData<T, decltype((void)T::FCData, void())>
  : std::true_type
{
};

template <typename T>
FlowControllerData*
getFCDataPtr(const T& a_t, std::false_type /*unused*/)
{
  amrex::ignore_unused(a_t);
  return nullptr;
}

template <typename T>
FlowControllerData*
getFCDataPtr(T& a_t, std::true_type /*unused*/)
{
  return &(a_t.FCData);
}
#endif
