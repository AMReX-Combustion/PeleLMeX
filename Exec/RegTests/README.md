# PeleLMeX testing

This folder contains a set of small-scale cases testing the various functionalities of PeleLMeX.
A number of these cases are automatically compiled or compiled & run upon pull request into the development branch,
The following provide a brief summary of each case and its intended purpose.

## EB\_EnclosedFlame
A 2D (3D) cylindrical flame enclosed in a cylindrical EB chamber. Testing the closed chamber algorithm with EB.

## EB\_EnclosedVortex
A 2D (3D) vortex enclosed in a cylindrical EB chamber. Testing EB for pure incompressible flows.

## EB\_FlowPastCylinder
A 2D (3D) inflow/outflow setup with either an EB cylinder in the middle of the flow or an EB
bump on a WallNoSlipAdiab domain boundary. Testing EB for incompressible flows and pure mixing cases (no reactions).

## EB\_PipeFlow
2D and 3D setup of both channel and pipe flow, either in an inflow/outflow setup or a 
periodc setup driven by a background pressure gradient.

## EnclosedFlame
A 2D (3D) cylindrical flame enclosed in a square (rectangular) domain. Testing the closed chamber algorithm.

## EnclosedInjection
A simple jet injection into an EB cylindrical chamber. Testing the closed chamber algorithm with mass injection.

## FlameSheet
A 2D (3D) harmonically perturbed flame sheet, initial solution from a Cantera simulation provided for 3 mechanisms
(drm19, dodecane\_lu and dodecane\_lu\_qss). This is the basis for weak scaling studies in PeleLMeX and test all the 
reactive pieces of the algorithm as well as transport options (Unity Lewis number, Soret effect, ...)

## HITDecay
A 3D decaying HIT case, where the initial solution is generated with a Passot-Pouquet spectrum. Testing basic incompressible
integration algorithm and Large Eddy Simulation implementation.

## HotBubble
Case of a 2D/2D-RZ/3D bubble of light gases (either hotter or lighter mixture composition) lifted under the effect of 
a gravity field. Testing of RZ algorthim and gravity forces.

## PeriodicCases
A set of 2D/3D periodic cases: convected vortex in different directions, convected temperature/mixture Gaussian bump,
pure diffusion of species/temperature. This is the base case to test PeleLMeX accuracy and order of convergence of
the spatial and temporal schemes.

## SprayTest
Basic testing of the coupling between PeleLMeX and PeleMP Lagrangian particle module.

## TaylorGreen
3D Taylor-Green vortex case. Enable testing of the different advection schemes available in PeleLMeX.

## TurbInflow
Injection of a 3D precomputed turbulent velocity field at the boundary of a PeleLMeX simulation. Testing the 
turbulence injection functionality.
