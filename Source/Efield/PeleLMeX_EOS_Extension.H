#ifndef EOS_EXT_H
#define EOS_EXT_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include "mechanism.H"

namespace pele {
namespace physics {
namespace eos {

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
charge_mass(amrex::Real zk[])
{
  CKCHRGMASS(zk);
}

} // namespace eos
} // namespace physics
} // namespace pele
#endif
