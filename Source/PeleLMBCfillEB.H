#ifndef PELE_BCFILL_EB
#define PELE_BCFILL_EB
#include "PeleLM_Index.H"
#include "pelelm_prob_EB.H"

#include <AMReX_PhysBCFunct.H>
#include <AMReX_Geometry.H>

using namespace amrex;

template <typename, typename T>
struct matchBCNormEB_sign {
    static_assert(std::integral_constant<T, false>::value,
                  "Second template parameter needs to be of function type.");
};

template <typename T, typename Ret, typename... Args>
struct matchBCNormEB_sign <T, Ret(Args...)> {
  private:
    template<typename C>
    static constexpr auto check(T*)
      -> typename std::is_same<decltype( std::declval<C>().bcnormal_eb( std::declval<Args>()... )), Ret>::type;
    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<T>(0)) type;

  public:
    static constexpr bool value = type::value;
};

template <typename T>
struct hasEBBCnormal {
  private:
    template<typename C>
    static constexpr auto check(T*)
    -> matchBCNormEB_sign<T, AMREX_GPU_DEVICE void(const amrex::Real *,
                                                   amrex::Real *,
                                                   const amrex::Real,
                                                   amrex::GeometryData const&,
                                                   ProbParm const&)>;
    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<T>(0)) type;

  public:
    static constexpr bool value = type::value;
};

// A default empty struct when the we don't provide a valid EBHandler
template <typename T, bool B = false>
struct PeleLMUserFillStateEB
{
    ProbParm const* lprobparm;
    T const ebh;

    AMREX_GPU_HOST
    constexpr PeleLMUserFillStateEB(ProbParm const* a_prob_parm, T const &a_ebh)
                                   : lprobparm(a_prob_parm), ebh(a_ebh) {}

    AMREX_GPU_DEVICE
    void operator()(
       int i, int j, int k,
       amrex::Array4<amrex::Real> const& ebState,
       const int firstcomp,
       const int numcomp,
       AMREX_D_DECL(amrex::Array4<const amrex::Real> const &ebfc_x,
                    amrex::Array4<const amrex::Real> const &ebfc_y,
                    amrex::Array4<const amrex::Real> const &ebfc_z),
       amrex::GeometryData const& geom,
       const amrex::Real time) const
    {
    }
};

template <typename T>
struct PeleLMUserFillStateEB <T, true>
{
    ProbParm const* lprobparm;
    T const ebh;

    AMREX_GPU_HOST
    constexpr PeleLMUserFillStateEB(ProbParm const* a_prob_parm, T const &a_ebh)
                                   : lprobparm(a_prob_parm), ebh(a_ebh) {}

    AMREX_GPU_DEVICE
    void operator()(
       int i, int j, int k,
       amrex::Array4<amrex::Real> const& ebState,
       const int firstcomp,
       const int numcomp,
       AMREX_D_DECL(amrex::Array4<const amrex::Real> const &ebfc_x,
                    amrex::Array4<const amrex::Real> const &ebfc_y,
                    amrex::Array4<const amrex::Real> const &ebfc_z),
       amrex::GeometryData const& geom,
       const amrex::Real time) const
    {
        // Get the EBface centroid coordinates
        const amrex::Real* dx = geom.CellSize();
        const amrex::Real* prob_lo = geom.ProbLo();
        const amrex::Real xcell[AMREX_SPACEDIM] = {AMREX_D_DECL(prob_lo[0] + (i + 0.5) * dx[0],
                                                                prob_lo[1] + (j + 0.5) * dx[1],
                                                                prob_lo[2] + (k + 0.5) * dx[2])};
        const amrex::Real xface[AMREX_SPACEDIM] = {AMREX_D_DECL(xcell[0] + ebfc_x(i,j,k) * dx[0],
                                                                xcell[1] + ebfc_y(i,j,k) * dx[1],
                                                                xcell[2] + ebfc_z(i,j,k) * dx[2])};

        //TODO : would be practical to have the current state at the EBface ...
        amrex::Real stateExt[NVAR] = {0.0};

        // User-defined fill function
        ebh.bcnormal_eb(xface, stateExt, time, geom, *lprobparm);

        // Extract requested entries
        for (int n = 0; n < numcomp; n++) {
           ebState(i,j,k,n) = stateExt[firstcomp+n];
        }
    }
};

template <typename, typename T>
struct matchEBBCSignature {
    static_assert(std::integral_constant<T, false>::value,
                  "Second template parameter needs to be of function type.");
};

template <typename T, typename Ret, typename... Args>
struct matchEBBCSignature <T, Ret(Args...)> {
  private:
    template<typename C>
    static constexpr auto check(T*)
      -> typename std::is_same<decltype( std::declval<C>().bctype_eb( std::declval<Args>()... ) ), Ret>::type;
    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<T>(0)) type;

  public:
    static constexpr bool value = type::value;
};

template <typename T>
struct hasEBBCType {
  private:
    template<typename C>
    static constexpr auto check(T*)
    -> matchEBBCSignature<T, AMREX_GPU_DEVICE void(int, amrex::GeometryData const&, ProbParm const&)>;
    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<T>(0)) type;

  public:
    static constexpr bool value = type::value;
};
#endif
