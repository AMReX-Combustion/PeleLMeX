## HIT
Two different turbulent options are possible (HIT decay and maintained HIT).

## HITDecay
A 3D decaying Homogeneous Isotropic Turbulence (HIT) case, where the initial solution is generated with a Passot-Pouquet spectrum. Testing the basic incompressible
integration algorithm and Large Eddy Simulation implementation.

## TurbForce Incompressible
A periodic domain that can be used to generate and maintain HIT. This is typically used to generate the initial turbulence profile used in maintained HIT flame sheets
(see Exec/RegTests/FlameSheet input file flamesheet-drm19-HITForced-3d.inp).

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