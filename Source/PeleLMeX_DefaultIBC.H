// Problem-specific Default Initial Boundary Conditions (IBC)
#ifndef PROBLEM_IBC_H
#define PROBLEM_IBC_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

#include <PMFData.H>

struct ProbParmDefault
{
  amrex::Real P_mean = 101325.0;
};

struct DefaultProbIBC
{
  /**
   * \brief Initialize data at the beginning of the simulation
   *
   * \param i                 i index.
   * \param j                 j index.
   * \param k                 k index.
   * \param is_incompressible incompressible flag.
   * \param state             state Array4.
   * \param aux               auxiliary Array4.
   * \param geomdata          geometry data.
   * \param prob_parm         problem specific parameters.
   * \param pmf_data          PMF (Premixed Flame) data
   */
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void
  initdata(
    int /*i*/,
    int /*j*/,
    int /*k*/,
    int /*is_incompressible*/,
    amrex::Array4<amrex::Real> const& /*state*/,
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDefault const& /*prob_parm*/,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    /** Notes:
     *   1. This function is called at the beginning of the simulation
     *   2. The default implementation does nothing
     *   3. Specific problem implementations should override this function
     */
  }

  /**
   * \brief Define the domain boundary conditions
   *
   * \param x                 the face center coordinates.
   * \param m_nAux            number of auxiliary variables.
   * \param s_in              incoming state in first cell(as a standard C++ array)
   * \param s_ext             outgoing state (as a standard C++ array)
   * \param idir              direction of the face.
   * \param sgn               sign of the face normal.
   * \param time              current time
   * \param geomdata          geometry data.
   * \param prob_parm         problem specific parameters.
   * \param pmf_data          PMF (Premixed Flame) data
   */
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void
  bcnormal(
    const amrex::Real* /*x[AMREX_SPACEDIM]*/,
    const int /*m_nAux*/,
    amrex::Real* /*s_ext[NVAR]*/,
    const int /*idir*/,
    const int /*sgn*/,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDefault const& /*prob_parm*/,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    /** Notes:
     *   1. This function handles the boundary conditions on the domain faces
     *   2. The default implementation does nothing
     *   3. Specific problem implementations should override this function
     */
  }


  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void
  zero_visc(
    int /*i*/,
    int /*j*/,
    int /*k*/,
    amrex::Array4<amrex::Real> const& /*beta*/,
    amrex::GeometryData const& /*geomdata*/,
    amrex::Box const& /*domainBox*/,
    const int /*dir*/,
    const int /*beta_comp*/,
    const int /*nComp*/)
  {
    /** Notes:
     *   1. This function can modify the diffusion coefficient
     *      on the domain faces in order to represent walls on a domain
     *      face otherwise marked as inflow
     *   2. To know for which state component the function is called, the
     *      logic is as follows:
     *       - if beta_comp == 0 and nComp == NUM_SPECIES => call for species
     *       - if beta_comp == 0 and nComp == 1 => call for temperature
     *   3. The default implementation does nothing
     */
  }

}; // namespace DefaultProbIBC

#endif
