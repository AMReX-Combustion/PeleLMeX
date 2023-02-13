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
  - name: Ray Grout
    affiliation: 1
  - name: Emmanuel Motheau
    affiliation: 2
  - name: Andy Nonaka
    affiliation: 2
  - name: Landon Owen
    affiliation: 3
  - name: Bruce Perry
    affiliation: 1
  - name: Jon Rood
    affiliation: 1
  - name: Nicolas Wimer
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
date: 28 December 2022
bibliography: paper.bib
---

# Summary

PeleLMeX evolves chemically reacting low Mach number flows with block-structured adaptive mesh refinement (AMR). 
The code is built upon the AMReX [@AMReX] library which provides the underlying data structures and tools to manage 
and operate on them across massively parallel computing architectures. Together with its compressible flow counterpart 
PeleC [@PeleC], the thermo-chemistry library PelePhysics and the multi-physics library PeleMP, it forms the Pele suite of 
open-source reactive flow simulation codes.

PeleLMeX uses a finite volume approach to solve the multi-species reacting Navier-Stokes equations in 
their low Mach number limit [@Day:2000], where the characteristic fluid velocity is small compared to the speed of sound, 
and the effect of acoustic wave propagation is unimportant to the overall dynamics of the system. Accordingly, 
acoustic wave propagation can be mathematically removed from the equations of motion, allowing for a numerical time 
step based on an advective CFL condition.
This low Mach number limit mathematically translates into a constraint on the divergence of the velocity field [@Majda:1986]. The 
momemtum equation is then solved for using a predictor/corrector method initially developed for incompressible flows [@Almgren1998]
and later extended to reactive, variable-density flows [@Pember:1998]. In the low Mach framework, the thermodynamic pressure is 
uniform in space but can evolve in time when simulating closed domains with chemical reactions and additional mass injections [@Nonaka18].
PeleLMeX uses an iterative Spectral Defered Correction (SDC) time advancement scheme [@Nonaka12;@Nonaka18] to ensure a tight coupling
of the fast diffusion/reaction and the comparatively slow advection, while iteratively enforcing 
the low Mach number contraint.
Advection terms are treated explicitly using second-order Godunov schemes [@AMReX-Hydro], diffusion terms are treated
semi-implicitly with a Crank-Nicholson scheme and the often stiffer reaction term is obtained using a fully implicit 
Backward Differentiation Formula schemes (specifically, the CVODE integrator [@balos2021enabling] of the Sundials
suite [@SUNDIALS;@gardner2022sundials]). The resolution of the linear systems arising in the implicit diffusion and velocity projections are
handled using AMReX's native geometric multigrid (GMG) solver, but can also be transfered to HYPRE [@Hypre2002] if GMG fails.
PeleLMeX relies on a non-subcycling approach to advance the numerical solution on an AMR hierarchy, where all the levels
are advanced together using the same time step, the size of which is precribed by a CFL condition accross all the levels. The consistency of
the numerical fluxes at coarse/fine interfaces is ensured by averaging down fluxes from fine to coarse levels.

In addition, PeleLMeX uses an Embedded Boundary (EB) approach to represent complex geometries: an arbitrary surface can 
be intersected with the Cartesian matrix of uniform cells, and the numerical stencils are modified near cells that are cut 
by the EB. Redistribution schemes [@Berger:2021] are then used for the explicit advection and diffusion updates in order to alleviate the 
constraint associated with small cut cells. Through its dependency to the multi-physics library PeleMP, PeleLMeX also inherits 
the ability to include Lagrangian sprays as well as soot and radiation models. 

PeleLMeX is written in C++ and is built upon the AMReX [@AMReX] library from which it inherits its parallel paradigm.
It uses a MPI+X approach where MPI is used to distribute AMR grid patches across CPU ranks and each grid can be further divided into 
logical tiles spread across threads using OpenMP for multi-core CPU machines, or spread across GPU threads using CUDA/HIP/SYCL 
on GPU-accelerated machines.

# Statement of Need

While there exist several reactive flow Direct Numerical Simulation codes, PeleLMeX presents a unique set of features. 
From its inception, under the name LMC in the early 2000, the motivation was to combine an AMR approach with a low Mach number 
formulation to achieve high performances from a small desktop stations to the world largest supercomputer, and to this day
it remains the only publicly available code to offers these features. Recent code developments focused on enabling
massively parallel simulations at scale on high-performance accelerated computer architectures to tackle the challenging
requirements of fundamental and applied combustion research, as well as extending the solver modeling capabilities by including
Large Eddy Simulation (LES) closure models and support for data-driven combustion models [@Perry:2021].

PeleLMeX is predominantly used to study the fine scale interactions between turbulence and chemical reactions occuring in many
combustion applications. A better understanding of these interactions is the basis for developing accurate modeling approaches
that can be used to design the next generation of low-emission combustion devices.

# Acknowledgments

This research was supported by the Exascale Computing Project (ECP), Project Number: 17-SC-20-SC, a collaborative effort of two DOE 
organizations -- the Office of Science and the National Nuclear Security Administration -- responsible for the planning and 
preparation of a capable exascale ecosystem -- including software, applications, hardware, advanced system engineering, and 
early testbed platforms -- to support the nation's exascale computing imperative.

# References

