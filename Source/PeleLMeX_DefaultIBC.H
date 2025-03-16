// Problem-specific Default Initial Boundary Conditions struct (IBC)
#ifndef PROBLEM_IBC_H
#define PROBLEM_IBC_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

#include <PMFData.H>

struct ProbParmDefault
{
  /**
   * PeleLMeX always needs a mean thermodynamic pressure
   */
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
  static void initdata(
    int /*i*/,
    int /*j*/,
    int /*k*/,
    int /*is_incompressible*/,
    amrex::Array4<amrex::Real> const& /*state*/,
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDefault const& /*prob_parm*/,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/){
    amrex::Abort(
      "Default empty initdata called, please override in your problem specific "
      "struct !")
    /** Notes:
     *   1. This function is called at the beginning of the simulation
     *   2. The default implementation does nothing
     *   3. Specific problem implementations should override this function
     */
  }

  /**
   * \brief Define the domain face boundary conditions
   *
   * \param x                 the face center coordinates.
   * \param s_in              state in first cell (as ptr into Array4)
   * \param s_ext             outgoing state (as a standard C++ array)
   * \param idir              direction of the face.
   * \param sgn               sign of the face normal.
   * \param time              current time
   * \param geomdata          geometry data.
   * \param prob_parm         problem specific parameters.
   * \param pmf_data          PMF (Premixed Flame) data
   */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void bcnormal(
    const amrex::Real* /*x[AMREX_SPACEDIM]*/,
    const amrex::CellData<amrex::Real>& /*s_in*/,
    amrex::Real* /*s_ext[NVAR]*/,
    const int /*idir*/,
    const int /*sgn*/,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDefault const& /*prob_parm*/,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/){
    amrex::Abort(
      "Default empty bcnormal called, please override in your problem specific "
      "struct !")
    /** Notes:
     *   1. This function handles the boundary conditions on the domain faces
     *   2. The default implementation does nothing
     *   3. Specific problem implementations should override this function
     */
  }

  /**
   * \brief Alter the scalars face-centered diffusion coefficients
   *
   * \param i                 i index.
   * \param j                 j index.
   * \param k                 k index.
   * \param beta              face-centered Beta Array4.
   * \param geomdata          geometry data.
   * \param domainBox         the computational domain box.
   * \param dir               direction of the face.
   * \param beta_comp         first components in Beta MultiFab.
   * \param nComp             number of components to change.
   * \param prob_parm         problem specific parameters.
   */
  AMREX_GPU_DEVICE AMREX_FORCE_INLINE static void zero_visc(
    int /*i*/,
    int /*j*/,
    int /*k*/,
    amrex::Array4<amrex::Real> const& /*beta*/,
    amrex::GeometryData const& /*geomdata*/,
    amrex::Box const& /*domainBox*/,
    const int /*dir*/,
    const int /*beta_comp*/,
    const int /*nComp*/,
    ProbParmDefault const& /*prob_parm*/)
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
     *   4. To only modify at the domain faces, one needs to compare the i, j or
     * k index matching the idir direction against the domainBox
     */
  }

  /**
   * \brief Set the embedded boundary type
   *
   * \param x                 the EB face center coordinates.
   * \param EBflagType        the outgoing embedded boundary type.
   * \param geomdata          geometry data.
   * \param prob_parm         problem specific parameters.
   */
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void setEBType(
    const amrex::Real* /*xEBface[AMREX_SPACEDIM]*/,
    int& /*EBflagType*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDefault const& /*prob_parm*/)
  {
    /** Notes:
     *   This function is under development and will be updated in the
     *   near future.
     */
  }

  /**
   * \brief Set the embedded boundary state value
   *
   * \param xEBface           the EB face center coordinates.
   * \param xEBfaceNormal     the EB face normal vector components.
   * \param s_in              the state in the cell the face belong to.
   * \param s_ext             outgoing state (as a standard C++ array).
   * \param time              current time.
   * \param geomdata          geometry data.
   * \param prob_parm         problem specific parameters.
   */
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void setEBState(
    const amrex::Real* /*xEBface[AMREX_SPACEDIM]*/,
    const amrex::Real* /*xEBfaceNormal[AMREX_SPACEDIM]*/,
    const amrex::CellData<const amrex::Real>& /*s_in*/,
    amrex::Real* /*s_ext[NVAR]*/,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDefault const& /*prob_parm*/)
  {
    /** Notes:
     *   This function is under development and will be updated in the
     *   near future.
     */
  }
}; // namespace DefaultProbIBC

#endif
