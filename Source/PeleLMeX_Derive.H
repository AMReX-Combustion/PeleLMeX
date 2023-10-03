#ifndef PELELM_Derive_H
#define PELELM_Derive_H
#include <list>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>

// Mimic AMReX Derive behavior but without StateDescriptor
// Forward declaration of PeleLM
class PeleLM;

using PeleLMDeriveFunc = void (*)(
  PeleLM*,
  const amrex::Box&,
  amrex::FArrayBox&,
  int,
  int,
  const amrex::FArrayBox&,
  const amrex::FArrayBox&,
  const amrex::FArrayBox&,
  const amrex::Geometry&,
  amrex::Real,
  const amrex::Vector<amrex::BCRec>&,
  int);

using PeleLMMFVecDeriveFunc = void (*)(
  amrex::Vector<amrex::MultiFab*>,
  int,
  int,
  amrex::Vector<const amrex::MultiFab*>,
  amrex::Vector<const amrex::MultiFab*>,
  amrex::Vector<const amrex::MultiFab*>,
  const amrex::Geometry&,
  amrex::Real,
  const amrex::Vector<amrex::BCRec>&);
class PeleLMDeriveRec
{
  friend class PeleLMDeriveList;

public:
  using DeriveBoxMap = amrex::Box (*)(const amrex::Box&);
  static amrex::Box TheSameBox(const amrex::Box& box) noexcept;
  static amrex::Box GrowBoxByOne(const amrex::Box& box) noexcept;

  /**
   * \brief The destructor.
   */
  ~PeleLMDeriveRec();

  /**
   * \brief The name of the derived type.
   */
  const std::string& name() const noexcept;

  /**
   * \brief The names of components
   *
   * \param comp
   */
  const std::string& variableName(int comp) const noexcept;

  /**
   * \brief The index of a variable
   *
   * \param name
   */
  int variableComp(const std::string& a_name) const noexcept;

  /**
   * \brief The IndexType of the derived type.
   */
  amrex::IndexType deriveType() const noexcept;

  /**
   * \brief The DeriveFunc used to calculate the derived type.
   */
  PeleLMDeriveFunc derFunc() const noexcept;

  DeriveBoxMap boxMap() const noexcept;

  /**
   * \brief Type of interpolater to use in computing derived type.
   */
  amrex::Interpolater* interp() const noexcept;

  /**
   * \brief Number of components in the derived type.
   */
  int numDerive() const noexcept;

  /**
   * \brief Constructor.
   *
   * \param name
   * \param result_type
   * \param nvar_derive
   * \param der_func
   * \param box_map
   * \param interp
   */
  PeleLMDeriveRec(
    std::string name,
    amrex::IndexType result_type,
    int nvar_derive,
    PeleLMDeriveFunc der_func,
    DeriveBoxMap box_map,
    amrex::Interpolater* interp = &amrex::pc_interp);

  /**
   * \brief Constructor.
   *
   * \param name
   * \param result_type
   * \param nvar_derive
   * \param var_names
   * \param der_func
   * \param box_map
   * \param interp
   */
  PeleLMDeriveRec(
    std::string name,
    amrex::IndexType result_type,
    int nvar_derive,
    amrex::Vector<std::string>& var_names,
    PeleLMDeriveFunc der_func,
    DeriveBoxMap box_map,
    amrex::Interpolater* interp = &amrex::pc_interp);

  /**
   * \brief Constructor.
   *
   * \param name
   * \param result_type
   * \param nvar_derive
   * \param der_func
   * \param box_map
   * \param interp
   */
  PeleLMDeriveRec(
    std::string name,
    amrex::IndexType result_type,
    int nvar_derive,
    DeriveBoxMap box_map,
    amrex::Interpolater* interp = &amrex::pc_interp);

  /**
   * \brief Constructor.
   *
   * \param name
   * \param result_type
   * \param nvar_derive
   * \param var_names
   * \param der_func
   * \param box_map
   * \param interp
   */
  PeleLMDeriveRec(
    std::string name,
    amrex::IndexType result_type,
    int nvar_derive,
    amrex::Vector<std::string>& var_names,
    DeriveBoxMap box_map,
    amrex::Interpolater* interp = &amrex::pc_interp);

private:
  //! Name of derived quantity.
  std::string derive_name;

  //! Names of derived variables
  amrex::Vector<std::string> variable_names;

  //! Type of derived quantity.
  amrex::IndexType der_type;

  //! Number of components in derived quantity.
  int n_derive = 0;

  //! Function that computes derived quantity from state variables.
  PeleLMDeriveFunc func = nullptr;

  //! Interpolater for mapping crse grid derived data to finer levels.
  amrex::Interpolater* mapper = nullptr;

  //! Box mapper that specifies constituent region given derived region.
  DeriveBoxMap bx_map;
};

/**
 * \brief A list of PeleLMDeriveRecs.
 *
 * PeleLMDeriveList manages and provides access to the list of PeleLMDeriveRecs.
 */
class PeleLMDeriveList
{
public:
  /**
   * \brief The default constructor.
   */
  PeleLMDeriveList();

  PeleLMDeriveList(const PeleLMDeriveList&) = delete;
  PeleLMDeriveList& operator=(const PeleLMDeriveList&) = delete;

  /**
   * \brief Determines whether quantity identified by \<name\> is in the
   * registry.
   *
   * \param name
   */
  bool canDerive(const std::string& name) const;

  /**
   * \brief Access the particular record in registry.
   *
   * \param name
   */
  const PeleLMDeriveRec* get(const std::string& name) const;

  /**
   * \brief Adds another entry to the registry.
   *
   * \param name
   * \param result_type
   * \param nvar_derive
   * \param der_func
   * \param box_map
   * \param interp
   */
  void add(
    const std::string& name,
    amrex::IndexType result_type,
    int nvar_der,
    PeleLMDeriveFunc der_func,
    PeleLMDeriveRec::DeriveBoxMap bx_map,
    amrex::Interpolater* interp = &amrex::pc_interp);

  /**
   * \brief Adds another entry to the registry.
   *
   * \param name
   * \param result_type
   * \param nvar_derive
   * \param var_names
   * \param der_func
   * \param box_map
   * \param interp
   */
  void add(
    const std::string& name,
    amrex::IndexType result_type,
    int nvar_der,
    amrex::Vector<std::string>& vars,
    PeleLMDeriveFunc der_func,
    PeleLMDeriveRec::DeriveBoxMap bx_map,
    amrex::Interpolater* interp = &amrex::pc_interp);

  /**
   * \brief Adds another entry to the registry.
   *
   * \param name
   * \param result_type
   * \param nvar_derive
   * \param der_func
   * \param box_map
   * \param interp
   */
  void add(
    const std::string& name,
    amrex::IndexType result_type,
    int nvar_der,
    PeleLMDeriveRec::DeriveBoxMap bx_map,
    amrex::Interpolater* interp = &amrex::pc_interp);

  /**
   * \brief Adds another entry to the registry.
   *
   * \param name
   * \param result_type
   * \param nvar_derive
   * \param var_names
   * \param der_func
   * \param box_map
   * \param interp
   */
  void add(
    const std::string& name,
    amrex::IndexType result_type,
    int nvar_der,
    amrex::Vector<std::string>& vars,
    PeleLMDeriveRec::DeriveBoxMap bx_map,
    amrex::Interpolater* interp = &amrex::pc_interp);

  std::list<PeleLMDeriveRec>& dlist();

  void clear() { lst.clear(); }

private:
  std::list<PeleLMDeriveRec> lst;
};
#endif
