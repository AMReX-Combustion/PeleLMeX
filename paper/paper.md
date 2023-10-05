---
title: 'PeleLMeX: an AMR Low Mach Number Reactive Flow Simulation Code without level sub-cycling'
tags:
  - C++
  - adaptive mesh refinement
  - hydrodynamics
  - combustion
  - reactions
  - CFD
  - low-Mach number
authors:
  - name: Lucas Esclapez
    orcid: 0000-0002-2438-7292
    affiliation: 1
  - name: Marc Day
    orcid: 0000-0002-1711-3963
    affiliation: 1
  - name: John Bell
    orcid: 0000-0002-5749-334X
    affiliation: 2
  - name: Anne Felden
    affiliation: 2
  - name: Candace Gilet
    affiliation: 4
  - name: Ray Grout
    affiliation: 1
  - name: Marc Henry de Frahan
    orcid: 0000-0001-7742-1565
    affiliation: 1
  - name: Emmanuel Motheau
    affiliation: 2
  - name: Andy Nonaka
    affiliation: 2
  - name: Landon Owen
    affiliation: 3
  - name: Bruce Perry
    orcid: 0000-0002-9150-8103
    affiliation: 1
  - name: Jon Rood
    affiliation: 1
    orcid: 0000-0002-7513-3225
  - name: Nicolas Wimer
    orcid: 0000-0001-5083-0799
    affiliation: 1
  - name: Weiqun Zhang
    affiliation: 2
affiliations:
  - name: High Performance Algorithms and Complex Fluids, National Renewable Energy Laboratory, USA
    index: 1
  - name: Center for Computational Sciences and Engineering, Lawrence Berkeley National Laboratory, USA
    index: 2
  - name: Sandia National Laboratory, USA
    index: 3
  - name: Independent Researcher, USA
    index: 4
date: 28 December 2022
bibliography: paper.bib
---

# Summary

PeleLMeX simulates chemically reacting low Mach number flows with block-structured adaptive mesh refinement (AMR).
The code is built upon the AMReX [@AMReX] library, which provides the underlying data structures and tools to manage
and operate on them across massively parallel computing architectures. PeleLMeX algorithmic features are inherited from its
predecessor PeleLM [@PeleLM] but key improvements allow representation of more complex physical processes. Together with its compressible
flow counterpart PeleC [@PeleC], the thermo-chemistry library PelePhysics and the multi-physics library PeleMP, it forms
the Pele suite of open-source reactive flow simulation codes.

PeleLMeX uses a finite volume approach to solve the multi-species reacting Navier-Stokes equations in
their low Mach number limit [@Day:2000], where the characteristic fluid velocity is small compared to the speed of sound,
and the effect of acoustic wave propagation is unimportant to the overall dynamics of the system. Accordingly,
acoustic wave propagation can be mathematically removed from the equations of motion, allowing for a numerical time
step based on an advective CFL condition.
This low Mach number limit mathematically translates into a constraint on the divergence of the velocity field [@Majda:1986]. The
momentum equation is then solved for using a predictor/corrector method initially developed for incompressible flows [@Almgren1998]
and later extended to reactive, variable-density flows [@Pember:1998]. In the low Mach framework, the thermodynamic pressure is
uniform in space but can evolve in time when simulating closed domains with chemical reactions and additional mass injections [@Nonaka18].
PeleLMeX uses an iterative Spectral Deferred Correction (SDC) time advancement scheme [@Nonaka12;@Nonaka18] to ensure a tight coupling
of the fast diffusion/reaction and the comparatively slow advection, while iteratively enforcing
the low Mach number constraint.
Advection terms are treated explicitly using second-order Godunov schemes [@AMReX-Hydro], diffusion terms are treated
semi-implicitly with a Crank-Nicholson scheme and the often stiffer reaction term is obtained using a fully implicit
Backward Differentiation Formula scheme (specifically, the CVODE integrator [@balos2021enabling] of the Sundials
suite [@SUNDIALS;@gardner2022sundials]). The solution of the linear systems arising in the implicit diffusion and velocity projections are
handled using AMReX's native geometric multigrid (GMG) solver, but can also be transferred to HYPRE [@Hypre2002] if GMG fails.
In contrast with PeleLM, PeleLMeX relies on a non-subcycling approach to advance the numerical solution on an AMR hierarchy,
where all the levels are advanced together using the same time step, the size of which is prescribed by a CFL condition across all the levels.
This distinctive feature drove the development of PeleLMeX as it enables extending the closed chamber algorithm described in
[@Nonaka18] to an AMR hierarchy and incorporating more complex physical processes such as flame/plasma interactions [@Esclapez:2020].

In addition, PeleLMeX uses an Embedded Boundary (EB) approach to represent complex geometries: an arbitrary surface can
be intersected with the Cartesian matrix of uniform cells, and the numerical stencils are modified near cells that are cut
by the EB. Redistribution schemes [@Giuliani:2022] are then used for the explicit advection and diffusion updates in order to alleviate the
constraint associated with small cut cells. Through its dependency to the multi-physics library PeleMP, PeleLMeX also inherits
the ability to include Lagrangian sprays as well as soot and radiation models.

PeleLMeX is written in C++ and is built upon the AMReX [@AMReX] library from which it inherits its parallel paradigm.
It uses a MPI+X approach where MPI is used to distribute AMR grid patches across CPU ranks and each grid can be further divided into
logical tiles spread across threads using OpenMP for multi-core CPU machines, or spread across GPU threads using CUDA/HIP/SYCL
on GPU-accelerated machines.

# Statement of Need

Several software tools for reactive flow simulations can found online (often with limited access), including unstructured body-fitted
solvers based on OpenFOAM [@Hassanaly:2018], the structured solver NGA2 [@NGA2], and the Sierra/Fuego solver [@Domingo:2003].
In contrast with the aforementioned solvers, PeleLMeX is fully publicly available and documented.
Its unique features consist in combining an AMR approach with a low
Mach number formulation to achieve high performances from a small desktop station to the world's largest supercomputer.
Recent code developments focused on enabling massively parallel simulations at scale on high-performance accelerated computer
architectures to tackle the challenging requirements of fundamental and applied combustion research, as well as extending
the solver modeling capabilities by including Large Eddy Simulation (LES) closure models and support for data-driven
combustion models [@Perry:2022].

PeleLMeX is intended for students, researchers and engineers interested in understanding complex combustion processes
by performing high fidelity simulations. Although it can be used to study laminar flames, its distinctive features make it
particularly attractive for studying the fine scale flame/turbulence interactions in combustion applications where AMR is necessary
to tackle the large scale separation and the computational resources available on the latest heterogeneous exascale
platform can be leveraged. In order to achieve energy, transport and industrial decarbonization, fuel-flexible combustion devices
must be designed and deployed to accommodate hydrogen, ammonia and a wide range of biofuels. In this context, PeleLMeX can prove a
valuable tool to study alternative fuels combustion characteristics, flame dynamics or pollutant formation mechanisms,
both in academic idealized cases [@Howarth:2023] as well as in device scale simulations [@Appukuttan:2023].

# Acknowledgments

This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE
organizations -- the Office of Science and the National Nuclear Security Administration -- responsible for the planning and
preparation of a capable exascale ecosystem -- including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing imperative. This research used resources of the Oak Ridge
Leadership Computing Facility, which is a DOE Office of Science User Facility supported under Contract DE-AC05-00OR22725 and computing
resources sponsored by the Department of Energy's Office of Energy Efficiency and Renewable Energy and located at the National Renewable Energy
Laboratory. This work was authored in part by the National Renewable Energy Laboratory, operated by Alliance for Sustainable Energy,
LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308. The views expressed in the article do not
necessarily represent the views of the DOE or the U.S. Government.

# References
