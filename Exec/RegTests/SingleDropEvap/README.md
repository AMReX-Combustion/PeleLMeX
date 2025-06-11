## Single Droplet Evaporation Validation
This case compares results from PeleLMeX against experimental literature data. Additional details and case descriptions are provided in the PelePhysics documentation at [https://amrex-combustion.github.io/PelePhysics/Spray.html#single-droplet-tests](https://amrex-combustion.github.io/PelePhysics/Spray.html#single-droplet-tests)

The input file `single-drop-evap.inp` is set up for the `WongLin` case.  This can be run using:
~~~
mpirun -np 4 ./<PeleLMeX_EXE> single-drop-evap.inp
~~~

All cases can be run by opening ``Validate.py`` and setting the case name from the table in the PelePhysics documentation listed above
~~~
case = TestCaseName()
~~~
then running 
~~~
python Validate.py
~~~

Case options include:
* `Nomura(471)` and `Nomura(741)`
* `WongLin()`
* `Daif()`
* `RungeHep()`, `RungeDec()`, `RungeMix()`, and `RungeJP8()`

Note that multicomponent evaporation is a work in progress as illustrated by the `RungeJP8` test case. 