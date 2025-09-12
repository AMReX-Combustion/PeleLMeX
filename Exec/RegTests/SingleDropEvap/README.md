## Single Droplet Evaporation Validation
This case compares results from PeleLMeX against experimental literature data. Additional details and case descriptions are provided in the PelePhysics documentation at [https://amrex-combustion.github.io/PelePhysics/Spray.html#single-droplet-tests](https://amrex-combustion.github.io/PelePhysics/Spray.html#single-droplet-tests)

There are two general input files:
* `single-drop-evap-mp.inp`
* `single-drop-evap-gcm.inp`
Both are set up for the `WongLin` case with either the original PeleMP liquid properties or the GCM liquid properties.  The case with the PeleMP liquid properties can be run by compiling with `SPRAY_FUEL_NUM = 2` and 
`SPRAY_GCM = FALSE`, then running:
~~~
mpirun -np 4 ./<PeleLMeX_EXE> single-drop-evap-mp.inp
~~~

Similarly, the case with the GCM liquid properties can be run by compiling with `SPRAY_FUEL_NUM = 2` and `SPRAY_GCM = TRUE`, then running:
~~~
mpirun -np 4 ./<PeleLMeX_EXE> single-drop-evap-gcm.inp
~~~

All cases provided in the PelePhysics documentation can be run by opening `Validate.py` and setting the `case_name` from the table in the PelePhysics documentation listed above, the model for the liquid properties, and the model for estimating saturated vapor pressure for the PeleMP model, which defaults to the Antoine fit. For example:
~~~
# Case to run
case_name = "Daif"

# Liquid properties model: "mp" or "gcm"
LiqPropsType = "mp"

# Psat model for PeleMP: "Antoine" or "Clasius-Clapeyron"
PeleMP_PsatModel = "Antoine"
~~~
then run
~~~
python Validate.py
~~~
Be sure to use the correct compile-time flag for `SPRAY_GCM` or the python script will thrown an error. 

Case options include:
* `Nomura`
* `WongLin`
* `Daif`
* `RungeHep`, `RungeDec`, `RungeMix`, and `RungeJP8`

Users can compare the results from the various tests/configurations by running the `CompareLiqPropsType.py` script with the desired `case_name`. The script will search the current directory for all available data for each case.

Note that multicomponent evaporation is a work in progress as illustrated by the `RungeJP8` test case. 

