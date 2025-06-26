## TurbForce FlameSheet

A 3D flame sheet, starting with an initial solution from a Cantera simulation provided. The domain is latterly periodic with an outflow at the top and a wall at the bottom. The flame starts high up in the domain and burns towards the wall. HIT is generated and maintained via an additional forcing function found in PelePhysics/Source/Utility/TurbForce.

# Workflow
For best results it is advised to generate an initial turbulence plotfile to get to generate the inital turbulence profile. This can be done cheaply using an incompressible simulation (see Exec/RegTests/TurbForce_incompressible). Then provide the resulting plotfile as a velocity_plotfile. Make sure then to set turbforce.time_offset to the time of that plotfile to keep the forcing term in sync with the turbulence.

# Citation
To cite the forcing scheme, please cite the [CAMCOS article](http://dx.doi.org/10.2140/camcos.2008.3.103)
```
@article{aspden2008analysis,
  title={{Analysis of implicit LES methods}},
  author={Aspden, Andrew and Nikiforakis, Nikos and Dalziel, Stuart and Bell, John},
  journal={Communications in Applied Mathematics and Computational Science},
  volume={3},
  number={1},
  pages={103--126},
  year={2008},
  publisher={Mathematical Sciences Publishers}
}
```