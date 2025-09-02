## Single Droplet Evaporation Validation
This case compares results from PeleLMeX against experimental literature data. Additional details and case descriptions are provided in the PelePhysics documentation at [https://amrex-combustion.github.io/PelePhysics/Spray.html#single-droplet-tests](https://amrex-combustion.github.io/PelePhysics/Spray.html#single-droplet-tests)

The input file `mp-single-drop-evap-heptane-decane.inp` is set up for the `WongLin` case with the original PeleMP liquid properties.  This can be run using:
~~~
mpirun -np 4 ./<PeleLMeX_EXE> mp-single-drop-evap-heptane-decane.inp
~~~

Similarly, the input file `gcm-single-drop-evap-heptane-decane.inp` is set up for the WongLin case with the group contribution method (GCM) liquid properties. This requires setting ``SPRAY_GCM=TRUE`` in the ``GNUmakefile``.

All cases can be run by opening ``Validate.py`` and setting the case name from the table in the PelePhysics documentation listed above
~~~
LiqPropsType = "mp" # "mp" or "gcm"
case = TestCaseName(LiqPropsType)
~~~
then running 
~~~
python Validate.py
~~~
Be sure to use the correct compile-time flag for ``SPRAY_GCM``. 

Case options include:
* `Nomura`
* `WongLin`
* `Daif`
* `RungeHep`, `RungeDec`, `RungeMix`, and `RungeJP8`

Note that multicomponent evaporation is a work in progress as illustrated by the `RungeJP8` test case. 