.. highlight:: rst

.. _sec:model:

The `PeleLMeX` Model
====================

In this section, we present the actual model that is evolved numerically by `PeleLMeX`, and the numerical algorithms
to do it.  There are many control parameters to customize the solution strategy and process, and in order to actually
set up and run specific problems with `PeleLMeX`, the user must specific the chemical model, and provide routines
that implement initial and boundary data and refinement criteria for the adaptive mesh refinement.  We discuss
setup and control of `PeleLMeX` in later sections. `PeleLMeX` is a non-subcycling version of `PeleLM` and as such, much of the following is identical to the model in `PeleLM` with a few key differences.

Overview of `PeleLMeX`
----------------------

`PeleLMeX` evolves chemically reacting low Mach number flows with block-structured adaptive mesh refinement (AMR). The code depends upon the `AMReX <https://github.com/AMReX-Codes/amrex>` library to provide the underlying data structures, and tools to manage and operate on them across massively parallel computing architectures. `PeleLMeX` also utilizes the source code and algorithmic infrastructure of `AMReX-Hydro <https://github.com/AMReX-Codes/AMReX-Hydro>`. `PeleLMeX` borrows heavily from `PeleLM <https://github.com/AMReX-Combustion/PeleLM>`.  The core algorithms in `PeleLM` are described in the following papers:

* *A conservative, thermodynamically consistent numerical approach for low Mach number combustion. I. Single-level integration*, A. Nonaka, J. B. Bell, and M. S. Day, *Combust. Theor. Model.*, **22** (1) 156-184 (2018)

* *A Deferred Correction Coupling Strategy for Low Mach Number Flow with Complex Chemistry*, A. Nonaka, J. B. Bell, M. S. Day, C. Gilet, A. S. Almgren, and M. L. Minion, *Combust. Theory and Model*, **16** (6) 1053-1088 (2012)

* *Numerical Simulation of Laminar Reacting Flows with Complex Chemistry*, M. S. Day and J. B. Bell, *Combust. Theory Model* **4** (4) 535-556 (2000)

* *An Adaptive Projection Method for Unsteady, Low-Mach Number Combustion*, R. B. Pember, L. H. Howell, J. B. Bell, P. Colella, W. Y. Crutchfield, W. A. Fiveland, and J. P. Jessee, *Comb. Sci. Tech.*, **140** 123-168 (1998)

* *A Conservative Adaptive Projection Method for the Variable Density Incompressible Navier-Stokes Equations,* A. S. Almgren, J. B. Bell, P. Colella, L. H. Howell, and M. L. Welcome, *J. Comp. Phys.*, **142** 1-46 (1998)

The low Mach number flow equations
----------------------------------

`PeleLMeX` solves the reacting Navier-Stokes flow equations in the *low Mach number* regime, where the characteristic fluid velocity is small compared to the sound speed, and the effect of acoustic wave propagation is unimportant to the overall dynamics of the system. Accordingly, acoustic wave propagation can be mathematically removed from the equations of motion, allowing for a numerical time step based on an advective CFL condition, and this leads to an increase in the allowable time step of order :math:`1/M` over an explicit, fully compressible method (:math:`M` is the Mach number).  In this mathematical framework, the total pressure is decomposed into the sum of a spatially constant (ambient) thermodynamic pressure :math:`P_0` and a perturbational pressure, :math:`\pi({\vec x})` that drives the flow.  Under suitable conditions, :math:`\pi/P_0 = \mathcal{O} (M^2)`. 

The set of conservation equations specialized to the low Mach number regime is a system of PDEs with advection, diffusion and reaction (ADR) processes that are constrained to evolve on the manifold of a spatially constant :math:`P_0`:

.. math::

    &\frac{\partial (\rho \boldsymbol{u})}{\partial t} + 
    \nabla \cdot \left(\rho  \boldsymbol{u} \boldsymbol{u} + \tau \right)
    = -\nabla \pi + \rho \boldsymbol{F},\\
    &\frac{\partial (\rho Y_m)}{\partial t} +
    \nabla \cdot \left( \rho Y_m \boldsymbol{u}
    + \boldsymbol{\mathcal{F}}_{m} \right)
    = \rho \dot{\omega}_m,\\
    &\frac{ \partial (\rho h)}{ \partial t} +
    \nabla \cdot \left( \rho h \boldsymbol{u}
    + \boldsymbol{\mathcal{Q}} \right) = 0 ,

where :math:`\rho` is the density, :math:`\boldsymbol{u}` is the velocity, :math:`h` is the mass-weighted enthalpy, :math:`T` is temperature and :math:`Y_m` is the mass fraction of species :math:`m`. :math:`\dot{\omega}_m` is the molar production rate for species :math:`m`, the modeling of which will be described later in this section. :math:`\tau` is the stress tensor, :math:`\boldsymbol{\mathcal{Q}}` is the heat flux and :math:`\boldsymbol{\mathcal{F}}_m` are the species diffusion fluxes. These transport fluxes require the evaluation of transport coefficients (e.g., the viscosity :math:`\mu`, the conductivity :math:`\lambda` and the diffusivity matrix :math:`D`) which are computed using the library EGLIB, as will be described in more depth in the diffusion section. The momentum source, :math:`\boldsymbol{F}`, is an external forcing term.  For example, we have used :math:`\boldsymbol{F}` to implement a long-wavelength time-dependent force to establish and maintain quasi-stationary turbulence.

These evolution equations are supplemented by an equation of state for the thermodynamic pressure.  For example, the ideal gas law,

.. math::

    P_0(\rho,Y_m,T)=\frac{\rho \mathcal{R} T}{W}=\rho \mathcal{R} T
    \sum_m \frac{Y_m}{W_m} .  

In the above, :math:`W_m` and :math:`W` are the species :math:`m`, and mean molecular weights, respectively.  To close the system we also require a relationship between enthalpy, species and temperature.  We adopt the definition used in the CHEMKIN standard,

.. math::

    h=\sum_m Y_m h_m(T)

where :math:`h_m` is the species :math:`m` enthalpy.  Note that expressions for :math:`h_m(T)` see <section on thermo properties> incorporate the heat of formation for each species.


Neither species diffusion nor reactions redistribute the total mass, hence we have :math:`\sum_m \boldsymbol{\mathcal{F}}_m = 0` and :math:`\sum_m \dot{\omega}_m = 0`. Thus, summing the species equations and using the definition :math:`\sum_m Y_m = 1` we obtain the continuity equation:

.. math::

    \frac{\partial \rho}{\partial t} + \nabla \cdot \rho \boldsymbol{u} = 0

This, together with the conservation equations form a differential-algebraic equation (DAE) system that describes an evolution subject to a constraint.  A standard approach to attacking such a system computationally is to differentiate the constraint until it can be recast as an initial value problem.  Following this procedure, we set the thermodynamic pressure constant in the frame of the fluid,

.. math::

    \frac{DP_0}{Dt} = 0

and observe that if the initial conditions satisfy the constraint, an evolution satisfying the above will continue to satisfy the constraint over all time.  Expanding this expression via the chain rule and continuity:

.. math::

    \nabla \cdot \boldsymbol{u} = \frac{1}{T}\frac{DT}{Dt}
    + W \sum_m \frac{1}{W_m} \frac{DY_m}{Dt} = S

The constraint here take the form of a condition on the divergence of the flow.  Note that the actual expressions to use here will depend upon the chosen models for evaluating the transport fluxes.

For the standard ideal gas EOS, the divergence constraint on velocity becomes:

.. math::

    \nabla \cdot \boldsymbol{u} &= \frac{1}{\rho c_p T} \left(\nabla \cdot \lambda\nabla T - \sum_m \boldsymbol{\Gamma_m} \cdot \nabla h_m \right) \\
    &- \frac{1}{\rho} \sum_m \frac{W}{W_m} \nabla \cdot \boldsymbol{\Gamma_m} + \frac{1}{\rho}\sum_m \left(\frac{W}{W_m} - \frac{h_m}{c_p T} \right) \dot \omega \equiv S .



Confined domain ambient pressure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In unconfined domains, the ambient pressure will remain constant in time, but for confined domains, this is not the case. Above, we assumed that :math:`p_0` was constant. If :math:`p_0` is a function of time, the pressure derivative term must be restored in the velocity divergence constraint as:

.. math::

    \nabla \cdot \boldsymbol{u} + \theta \frac{dp_0}{dt} = S ,

where :math:`\theta \equiv 1/(\Gamma_1 p_0)`, with :math:`\Gamma_1 = \partial ln(p)/\partial ln(\rho)|_s` being the first adiabatic exponent. :math:`\Gamma_1` depends on the composition and is not a constant. Both :math:`\theta` and :math:`\S` can be decomposed into mean and fluctuating components and the above equation can be rewritten as:

.. math::

    \nabla \cdot \boldsymbol{u} + (\overline \theta + \delta \theta)\frac{dp_0}{dt} = \overline S + \delta S,

where :math:`\overline \theta` and :math:`\overline S` are the mean values of :math:`\theta` and :math:`S` over the domain, and :math:`\delta \theta` and :math:`\delta S` are the perturbations off their respective means that both integrate to zero over the domain, by definition. This equation can be simplified by integrating over the domain volume:

.. math::
    
    \int_V \nabla \cdot \boldsymbol{u} dV + \int_V (\overline \theta + \delta \theta)\frac{dp_0}{dt} dV = \int_V (\overline S + \delta S) dV

Since the perturbations integrate to zero over the domain volume, the mean values are constants, and :math:`p_0` is only a function of time, the above simplifies to:

.. math::

    \frac{1}{V} \int_V \nabla \cdot \boldsymbol{u} dV + \overline \theta \frac{dp_0}{dt} = \overline S .

Solving for :math:`dp_0/dt` yields an evolution equation of :math:`p_0`:

.. math::

    \frac{dp_0}{dt} = \frac{1}{\overline \theta} \left(\overline S - \frac{1}{V} \int_A \boldsymbol{u} dA \right) ,

where we have used the divergence theorem to convert the volume integral into a surface integral over the domain boundaries: :math:`\int_V \nabla \cdot \boldsymbol{u} dV = \int_A \boldsymbol{u} dA`. The above pressure evolution is accomponied by a modified velocity constraint:

.. math::

    \nabla \cdot \boldsymbol{u} = \delta S - \delta \theta \frac{\overline S}{\overline \theta} - \left(1 + \frac{\theta}{\overline \theta} \right) \frac{1}{V} \int_A \boldsymbol{u} dA

The above equations hold for any fully enclosed or partially enclosed domain where there can be mass flowing into or out of the domain, but the next flowrate is non-zero and therefore the pressure should be adjusted in time. In a perfectly enclosed domain, where there is no mass in or out of the system, :math:`\int_A \boldsymbol{u} dA = 0` and the pressure evolution is simplified to:

.. math::

    \frac{dp_0}{dt} = \frac{\overline S}{\overline \theta} ,

and simplified velocity constraint,

.. math::

     \nabla \cdot \boldsymbol{u} = \delta S - \delta \theta \frac{\overline S}{\overline \theta} .

    
