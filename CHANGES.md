## v1.1.0

EB-Inflow capability was added in PR #515. To enable this capability, the optional
`setEBState` and `setEBType` problem-specific functions are being renamed and having the
interfaces updated, becoming `bcnormal_eb` and `bctype_eb`, respectively. See the
`pelelmex_prob.H` file in the EB_BackwardStepFlame case for examples and the new required
interface. Errors will be raised if isothermal EBs or EB inflows are used without the
new functions being properly defined. Note this capability is still under development and
may be subject to further change/enhancements.

## v1.0.0

Version 1.0.0 is meant to be the first stable release of PeleLMeX. This version
includes changes relative to prior versions (the latest being v0.25.4/v25.04) that
require modification of case files for compatibility. The major change is that
many user-defined functions have been gathered into a single `ProblemSpecificFunctions`
struct, which allows for defaults to be used without redefining these functions
for each new case. To update a case files from v0.25.04 to v1.0.0, the following
steps are required:

1. Define the `PeleLM::freeProbParm()` function in `pelelmex_prob.cpp`. This function
can be a no-op unless memory is explicitly allocated in the `PeleLM::readProbParm()`
function:

    ```
    void
    PeleLM::freeProbParm()
    {
    }
    ```

2. Add `#include <PeleLMeX_ProblemSpecificFunctions.H>` at the top of `pelelmex_prob.H`

3. Move the contents of the the `ProbParm` struct defined in the `pelelmex_prob_parm.H` to
a new `MyProbParm` struct defined in `pelelmex_prob.H` that inherits from `DefaultProbParm`:

    ```
    struct MyProbParm : public ProbParmDefault
    {
        // contents of your ProbParm
    };
    ```

    Delete the `pelelmex_prob_parm.H` file and remove the `#include <pelelmex_prob_parm.H>`
    from the top of `pelelmex_prob.H`.

4. Define a struct in `pelelmex_prob.H` that inherits from `DefaultProblemSpecificFunctions`.
Move your `pelelmex_initdata` function definition to a `static` function member of this
struct named `initdata`. Similarly, if you have nontrivial `bcnormal` or `zero_visc`
functions, move these function definitions to be `static` members of the new struct (though
these function definitions may now be optionally omitted if they do not do anything). Note,
if you have a zero_visc function, `MyProbParm` must be added to the interface as indicated
below and `MyProbParm` should replace `ProbParm` in the interfaces for the other functions.

    ```
    struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
    {
        AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE
        static void initdata(
            int i,
            int j,
            int k,
            int /*is_incompressible*/,
            amrex::Array4<amrex::Real> const& state,
            amrex::Array4<amrex::Real> const& /*aux*/,
            amrex::GeometryData const& geomdata,
            MyProbParm const& prob_parm,
            pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
        {
            // Contents of your pelelmex_initdata function
        }

        AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE
        static void bcnormal(
            const amrex::Real* /*x[AMREX_SPACEDIM]*/,
            const amrex::CellData<amrex::Real>& /*s_in*/,
            amrex::Real s_ext[NVAR],
            const int /*idir*/,
            const int sgn,
            const amrex::Real /*time*/,
            amrex::GeometryData const& /*geomdata*/,
            MyProbParm const& prob_parm,
            pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
        {
            // Contents of your bcnormal function
        }

        AMREX_GPU_DEVICE
        AMREX_FORCE_INLINE
        static void zero_visc(
            int i,
            int j,
            int k,
            amrex::Array4<amrex::Real> const& beta,
            amrex::GeometryData const& /*geomdata*/,
            amrex::Box const& domainBox,
            const int dir,
            const int beta_comp,
            const int nComp,
            MyProbParm const& /*prob_parm*/)
        {
            // Contents of your zero_visc function
        }
    };
    ```

    If you have a `bcnormal_aux` function, treat it similarly to `bcnormal`.

5. Adding the following lines to the end of `pelelmex_prob.H`:

    ```
    using ProblemSpecificFunctions = MyProblemSpecificFunctions;
    using ProbParm = MyProbParm;
    ```

6. If your case has a `PeleLMeX_EBUserDefined.H` file, move the `EBUserDefined`,
`setEBState`, and `setEBType` function definitions to also be `static` members of the
`MyProblemSpecificFunctions` struct and delete the `PeleLMeX_EBUserDefined.H` file.
The interface for these functions has been adjusted, so match the interfaces from
`PeleLMeX/Source/PeleLMeX_ProblemSpecificFunctions.H` but replace `ProbParmDefault`
with `MyProbParm`. Note that `setEBType` now sets an integer `EBflagType` (indicating
BC type) and a real `EBfaceFrac`, indicating the fraction of the EB face covered by
that type.

7. If your case has a `PeleLMeX_ProblemSpecificFunctions.cpp` file that defines
non-trivial `modify_ext_sources` or `set_ode_names` functions, move these function
definitions to also be `static` members of the `MyProblemSpecificFunctions` struct and
delete the `PeleLMeX_ProblemSpecificFunctions.cpp` file.

8. If your case has a `PeleLMeX_PatchFlowVariables.cpp` file that defines the
`patchFlowVaraibvles` function, move this function to be a static member of the
`MyProblemSpecificFunctions` struct and delete the `PeleLMeX_PatchFlowVariables.cpp` file.

Feel free to open a discussion on the PeleLMeX GitHub page if you run into issues. You can
also view changes in PeleLMeX Pull Request #477 to see how these changes were implemented
for the standard cases in PeleLMeX.
